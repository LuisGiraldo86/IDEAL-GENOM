"""Shared utility functions for genomic analysis pipelines."""

import gzip
import os
import psutil
import logging
import requests
import shutil
import zipfile
from pathlib import Path
from typing import Optional, List

logger = logging.getLogger(__name__)

def get_optimal_threads(reserve: int = 2, default: int = 10, max_threads: Optional[int] = None) -> int:
    """
    Calculate optimal thread count for genomic analysis operations.
    
    Determines the number of threads to use based on the CPU cores actually
    schedulable by this process (respecting affinity restrictions from
    taskset/cpuset cgroups/SLURM allocations where available), reserving
    some cores for system operations with robust fallback handling.
    
    Parameters
    ----------
    reserve : int, default=2
        Number of cores to reserve for system operations
    default : int, default=10
        Default thread count if CPU detection fails
    max_threads : int, optional
        Maximum number of threads to use (caps the result)
        
    Returns
    -------
    int
        Optimal number of threads to use (always >= 1)
        
    Examples
    --------
    >>> threads = get_optimal_threads()  # On 16-core system, returns 14
    >>> threads = get_optimal_threads(reserve=4)  # Returns 12
    >>> threads = get_optimal_threads(max_threads=8)  # Never exceeds 8
    """
    if reserve < 0:
        raise ValueError(f"reserve must be >= 0, got {reserve}")
    if default < 1:
        raise ValueError(f"default must be >= 1, got {default}")
    if max_threads is not None and max_threads < 1:
        raise ValueError(f"max_threads must be >= 1, got {max_threads}")
    
    # Try multiple methods to get CPU count
    cpu_count = None

    # Method 1: os.sched_getaffinity() - preferred, respects CPU affinity
    # restrictions (taskset, cpuset cgroups, SLURM/PBS allocations on a
    # shared node). os.cpu_count() ignores these and reports the total
    # cores on the machine, which can lead to requesting far more PLINK
    # threads than the process can actually schedule on.
    try:
        cpu_count = len(os.sched_getaffinity(0))
        if cpu_count and cpu_count > 0:
            logger.debug(f"CPU count from os.sched_getaffinity(0): {cpu_count}")
        else:
            cpu_count = None
    except (AttributeError, NotImplementedError, OSError) as e:
        # Not available on non-Linux platforms (e.g. macOS, Windows)
        logger.debug(f"os.sched_getaffinity() unavailable: {e}")
        cpu_count = None

    # Method 2: os.cpu_count() fallback
    if cpu_count is None:
        try:
            cpu_count = os.cpu_count()
            if cpu_count and cpu_count > 0:
                logger.debug(f"CPU count from os.cpu_count(): {cpu_count}")
            else:
                cpu_count = None
        except Exception as e:
            logger.warning(f"Failed to get CPU count from os.cpu_count(): {e}")
            cpu_count = None

    # Method 3: psutil fallback
    if cpu_count is None:
        try:
            cpu_count = psutil.cpu_count(logical=True)
            if cpu_count and cpu_count > 0:
                logger.debug(f"CPU count from psutil.cpu_count(): {cpu_count}")
            else:
                cpu_count = None
        except Exception as e:
            logger.warning(f"Failed to get CPU count from psutil: {e}")
            cpu_count = None
    
    # Calculate optimal threads
    if cpu_count is not None and cpu_count > 0:
        optimal = max(1, cpu_count - reserve)
    else:
        logger.warning(f"Could not detect CPU count, using default: {default}")
        optimal = default
    
    # Apply max_threads cap if specified
    if max_threads is not None:
        optimal = min(optimal, max_threads)
    
    logger.debug(f"Calculated optimal threads: {optimal} (CPU: {cpu_count}, reserve: {reserve})")
    return optimal

def _get_cgroup_available_memory_mb() -> Optional[float]:
    """
    Get memory actually available under a cgroup limit (Docker --memory,
    Kubernetes memory limits, SLURM --mem allocations), in MB.

    psutil.virtual_memory() reads host-wide /proc/meminfo and has no notion
    of cgroup memory limits, so inside a container/job restricted to less
    memory than the host it silently over-reports what's usable. Returns
    None if no cgroup limit is set (or cgroup files aren't readable), in
    which case the host-wide psutil figures are the correct source.
    """
    # cgroup v2
    try:
        max_path = Path('/sys/fs/cgroup/memory.max')
        current_path = Path('/sys/fs/cgroup/memory.current')
        if max_path.is_file() and current_path.is_file():
            raw_limit = max_path.read_text().strip()
            if raw_limit != 'max':
                limit_bytes = int(raw_limit)
                current_bytes = int(current_path.read_text().strip())
                return max(0.0, (limit_bytes - current_bytes) / (1024 * 1024))
    except (OSError, ValueError) as e:
        logger.debug(f"Could not read cgroup v2 memory limit: {e}")

    # cgroup v1
    try:
        limit_path = Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')
        usage_path = Path('/sys/fs/cgroup/memory/memory.usage_in_bytes')
        if limit_path.is_file() and usage_path.is_file():
            limit_bytes = int(limit_path.read_text().strip())
            # cgroup v1 uses a very large sentinel value for "no limit"
            if limit_bytes < (1 << 62):
                usage_bytes = int(usage_path.read_text().strip())
                return max(0.0, (limit_bytes - usage_bytes) / (1024 * 1024))
    except (OSError, ValueError) as e:
        logger.debug(f"Could not read cgroup v1 memory limit: {e}")

    return None

def get_available_memory(
    fraction: float = 2/3, 
    min_mb: int = 512, 
    max_mb: Optional[int] = None,
    safety_buffer_mb: int = 1024
) -> int:
    """
    Calculate available memory for genomic analysis operations with safety checks.
    
    Determines the amount of memory to allocate based on currently available
    system memory, using a configurable fraction with minimum and maximum limits
    to avoid system instability.
    
    Parameters
    ----------
    fraction : float, default=2/3
        Fraction of available memory to use (should be between 0 and 1)
    min_mb : int, default=512
        Minimum memory to allocate in MB
    max_mb : int, optional
        Maximum memory to allocate in MB (None for no limit)
    safety_buffer_mb : int, default=1024
        Safety buffer to always leave available for system (MB)
        
    Returns
    -------
    int
        Memory in MB to allocate
        
    Raises
    ------
    ValueError
        If fraction is not between 0 and 1, or other parameters are invalid
    RuntimeError
        If insufficient memory is available
        
    Examples
    --------
    >>> memory_mb = get_available_memory()  # Uses 2/3 of available memory
    >>> memory_mb = get_available_memory(fraction=0.5, max_mb=8192)  # Uses half, max 8GB
    """
    if not 0 < fraction <= 1:
        raise ValueError(f"fraction must be between 0 and 1, got {fraction}")
    if min_mb < 0:
        raise ValueError(f"min_mb must be >= 0, got {min_mb}")
    if max_mb is not None and max_mb < min_mb:
        raise ValueError(f"max_mb ({max_mb}) must be >= min_mb ({min_mb})")
    if safety_buffer_mb < 0:
        raise ValueError(f"safety_buffer_mb must be >= 0, got {safety_buffer_mb}")
    
    try:
        memory_info = psutil.virtual_memory()
        
        # Log memory information for debugging
        total_mb = memory_info.total / (1024 * 1024)
        available_mb = memory_info.available / (1024 * 1024)
        used_mb = memory_info.used / (1024 * 1024)
        
        logger.debug(f"System memory - Total: {total_mb:.0f}MB, Available: {available_mb:.0f}MB, Used: {used_mb:.0f}MB")

        # Cap by cgroup limit if this process is running under one (Docker
        # --memory, Kubernetes limits, SLURM --mem) — psutil only sees the
        # host-wide figures above.
        cgroup_available_mb = _get_cgroup_available_memory_mb()
        if cgroup_available_mb is not None and cgroup_available_mb < available_mb:
            logger.debug(f"Capping available memory to cgroup limit: {cgroup_available_mb:.0f}MB (host-wide: {available_mb:.0f}MB)")
            available_mb = cgroup_available_mb

        # Calculate usable memory after safety buffer
        usable_memory_mb = max(0, available_mb - safety_buffer_mb)
        
        # Apply fraction
        target_memory_mb = usable_memory_mb * fraction
        
        # Apply minimum constraint
        if target_memory_mb < min_mb:
            if usable_memory_mb < min_mb:
                raise RuntimeError(
                    f"Insufficient memory available. Need at least {min_mb}MB, "
                    f"but only {usable_memory_mb:.0f}MB usable after safety buffer."
                )
            target_memory_mb = min_mb
            logger.warning(f"Using minimum memory {min_mb}MB (would be {target_memory_mb:.0f}MB with fraction)")
        
        # Apply maximum constraint if specified
        if max_mb is not None and target_memory_mb > max_mb:
            target_memory_mb = max_mb
            logger.debug(f"Capped memory allocation to maximum {max_mb}MB")
        
        result = int(round(target_memory_mb, 0))
        logger.debug(f"Allocated memory: {result}MB (fraction: {fraction}, available: {available_mb:.0f}MB)")
        
        return result
        
    except Exception as e:
        logger.error(f"Failed to get memory information: {e}")
        raise RuntimeError(f"Could not determine available memory: {e}")

def count_file_lines(file_path: Path) -> int:
    """
    Count lines in a file efficiently.
    
    Uses a generator expression for memory-efficient line counting,
    suitable for large genomic data files.
    
    Parameters
    ----------
    file_path : Path
        Path to the file to count
        
    Returns
    -------
    int
        Number of lines in the file
        
    Raises
    ------
    FileNotFoundError
        If the file does not exist
    IOError
        If the file cannot be read
        
    Examples
    --------
    >>> from pathlib import Path
    >>> count = count_file_lines(Path('variants.bim'))
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)

def validate_input_file(file_path: Path, extensions: Optional[List[str]] = None) -> Path:
    """
    Validate that a file exists and optionally check its extension.
    
    Validates file existence and optionally ensures the file has one of the 
    specified extensions. Useful for genomic data files that must have specific
    formats (e.g., .vcf, .bim, .fam).
    
    Parameters
    ----------
    file_path : Path
        Path to the file to validate
    extensions : List[str], optional
        List of valid file extensions (including the dot, e.g., ['.vcf', '.vcf.gz']).
        If None, no extension validation is performed.
        
    Returns
    -------
    Path
        Validated file path
        
    Raises
    ------
    TypeError
        If file_path is not a Path object
    FileNotFoundError
        If the file does not exist
    IsADirectoryError
        If the path points to a directory instead of a file
    ValueError
        If the file extension is not in the allowed extensions list
        
    Examples
    --------
    >>> from pathlib import Path
    >>> vcf_file = validate_input_file(Path('data.vcf'), ['.vcf', '.vcf.gz'])
    >>> any_file = validate_input_file(Path('output.txt'))  # No extension check
    """
    if not isinstance(file_path, Path):
        raise TypeError(f"file_path must be a Path object, got {type(file_path)}")
    
    if not file_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {file_path}")
    
    if not file_path.is_file():
        raise IsADirectoryError(f"Path is not a file: {file_path}")
    
    # Validate extension if specified
    if extensions is not None:
        if not isinstance(extensions, (list, tuple)):
            raise TypeError("extensions must be a list or tuple of strings")
        
        # Check if any of the allowed extensions match
        file_extension = ''.join(file_path.suffixes)  # Handles .vcf.gz, .tar.gz, etc.
        
        # Also check just the last suffix for single extensions
        single_extension = file_path.suffix
        
        valid_extension = (file_extension in extensions or 
                          single_extension in extensions or
                          str(file_path).endswith(tuple(extensions)))
        
        if not valid_extension:
            raise ValueError(
                f"File '{file_path}' has invalid extension. "
                f"Expected one of {extensions}, got '{file_extension}' or '{single_extension}'"
            )
    
    return file_path

def validate_file_path(file_path: Path, must_exist: bool = True, must_be_file: bool = True) -> Path:
    """
    Generic file path validation with flexible requirements.
    
    Provides flexible validation for file paths with configurable requirements.
    Useful when you need different validation rules for different scenarios
    (e.g., input files that must exist vs. output files that may not exist yet).
    
    Parameters
    ----------
    file_path : Path
        Path to validate
    must_exist : bool, default=True
        If True, the path must already exist
    must_be_file : bool, default=True
        If True, the path must be a file (not a directory).
        Only checked if must_exist=True and the path exists.
        
    Returns
    -------
    Path
        Validated file path
        
    Raises
    ------
    TypeError
        If file_path is not a Path object
    FileNotFoundError
        If must_exist=True and the path does not exist
    IsADirectoryError
        If must_be_file=True and the path is a directory
        
    Examples
    --------
    >>> from pathlib import Path
    >>> # Validate existing input file
    >>> input_file = validate_file_path(Path('input.txt'))
    >>> 
    >>> # Validate output file path (may not exist yet)
    >>> output_file = validate_file_path(Path('output.txt'), must_exist=False)
    >>> 
    >>> # Validate path that could be file or directory
    >>> path = validate_file_path(Path('data'), must_be_file=False)
    """
    if not isinstance(file_path, Path):
        raise TypeError(f"file_path must be a Path object, got {type(file_path)}")
    
    # Check existence requirement
    if must_exist and not file_path.exists():
        raise FileNotFoundError(f"Path does not exist: {file_path}")
    
    # Check file type requirement (only if path exists)
    if file_path.exists() and must_be_file and not file_path.is_file():
        if file_path.is_dir():
            raise IsADirectoryError(f"Path is a directory, expected a file: {file_path}")
        else:
            raise ValueError(f"Path exists but is neither a file nor directory: {file_path}")
    
    return file_path

def validate_output_dir(output_dir: Path, create: bool = True) -> Path:
    """
    Validate and optionally create output directory.
    
    Parameters
    ----------
    output_dir : Path
        Path to the output directory
    create : bool, default=True
        If True, create the directory if it doesn't exist
        
    Returns
    -------
    Path
        Validated output directory path
        
    Raises
    ------
    FileNotFoundError
        If directory doesn't exist and create=False
    PermissionError
        If directory cannot be created due to permissions
        
    Examples
    --------
    >>> from pathlib import Path
    >>> output = validate_output_dir(Path('/data/results'))
    """
    if not output_dir.exists():
        if create:
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
    
    if not output_dir.is_dir():
        raise NotADirectoryError(f"Path exists but is not a directory: {output_dir}")
    
    return output_dir

def format_memory_size(bytes_size: int) -> str:
    """
    Format byte size into human-readable string.
    
    Parameters
    ----------
    bytes_size : int
        Size in bytes
        
    Returns
    -------
    str
        Formatted size string (e.g., '1.5 GB', '256 MB')
        
    Examples
    --------
    >>> format_memory_size(1536 * 1024 * 1024)
    '1.50 GB'
    >>> format_memory_size(512 * 1024)
    '512.00 KB'
    """
    size = float(bytes_size)
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} PB"

def get_system_resource_info() -> dict:
    """
    Get comprehensive system resource information.
    
    Returns
    -------
    dict
        Dictionary containing CPU, memory, and disk information
        
    Examples
    --------
    >>> info = get_system_resource_info()
    >>> print(f"Available memory: {info['memory']['available_mb']:.0f} MB")
    """
    try:
        # CPU information
        cpu_count_physical = psutil.cpu_count(logical=False)
        cpu_count_logical = psutil.cpu_count(logical=True)
        cpu_freq = psutil.cpu_freq()
        
        # Memory information
        memory = psutil.virtual_memory()
        swap = psutil.swap_memory()
        
        # Disk information for current working directory
        disk = psutil.disk_usage('.')
        
        return {
            'cpu': {
                'physical_cores': cpu_count_physical,
                'logical_cores': cpu_count_logical,
                'frequency_mhz': cpu_freq.current if cpu_freq else None,
            },
            'memory': {
                'total_mb': memory.total / (1024 * 1024),
                'available_mb': memory.available / (1024 * 1024),
                'used_mb': memory.used / (1024 * 1024),
                'percent_used': memory.percent,
            },
            'swap': {
                'total_mb': swap.total / (1024 * 1024),
                'used_mb': swap.used / (1024 * 1024),
                'percent_used': swap.percent,
            },
            'disk': {
                'total_gb': disk.total / (1024 * 1024 * 1024),
                'free_gb': disk.free / (1024 * 1024 * 1024),
                'used_gb': disk.used / (1024 * 1024 * 1024),
                'percent_used': (disk.used / disk.total) * 100,
            }
        }
    except Exception as e:
        logger.error(f"Failed to get system resource information: {e}")
        return {}

def download_file(url:str, local_filename: Path) -> None:

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    return

def unzip_file_flat(in_file: Path, target_file: str, out_dir: Path, remove_zip: bool = False) -> Path:
    """Extracts a specific file from a ZIP archive, decompresses it if it's a .gz file, and optionally deletes original files.

    Args:
        in_file (str): Path to the ZIP file.
        target_file (str): The file inside the ZIP to extract.
        out_dir (str): Directory where the extracted file will be saved.
        remove_zip (bool): If True, delete the original ZIP file after extraction.
        remove_gz (bool): If True, delete the .gz file after decompression.

    Returns:
        Path: Path to the final extracted file.
    """
    in_file = Path(in_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

    extracted_gz_path = out_dir / Path(target_file).name  # Target extracted .gz file

    try:
        with zipfile.ZipFile(in_file, "r") as zip_ref:
            if target_file in zip_ref.namelist():
                # Extract the .gz file from ZIP
                with zip_ref.open(target_file) as source, open(extracted_gz_path, "wb") as dest:
                    dest.write(source.read())
                print(f"Extracted: {extracted_gz_path}")
            else:
                print(f"File {target_file} not found in the archive.")
                return Path()

        # Optionally delete the ZIP file
        if remove_zip:
            in_file.unlink()
            print(f"Deleted ZIP file: {in_file}")

        return extracted_gz_path

    except zipfile.BadZipFile:
        print(f"Error: {in_file} is not a valid ZIP file.")
    except Exception as e:
        print(f"Unexpected error: {e}")

    return Path()  # Return None if extraction fails

def extract_gz_file(gz_file: Path, out_dir: Path, remove_gz: bool = False) -> Path:
    """Extracts a .gz file and saves the decompressed content in the same directory.

    Args:
        gz_file (str): Path to the .gz file.
        out_dir (str): Directory where the decompressed file will be saved.
        remove_gz (bool): If True, delete the .gz file after extraction.

    Returns:
        Path: Path to the extracted file.
    """
    gz_file = Path(gz_file)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

    # Define the decompressed file path (removes .gz extension)
    decompressed_file = out_dir / gz_file.stem

    try:
        with gzip.open(gz_file, "rb") as f_in, open(decompressed_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)  # Copy content from .gz to uncompressed file
        print(f"Decompressed: {decompressed_file}")

        if remove_gz:
            gz_file.unlink()  # Delete the .gz file
            print(f"Removed original .gz file: {gz_file}")

    except Exception as e:
        print(f"Error extracting {gz_file}: {e}")

    return decompressed_file