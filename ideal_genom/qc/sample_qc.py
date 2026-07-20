"""
Module to perform sample quality control
"""

import os
import psutil
import gc
import logging
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colormaps
import seaborn as sns

from ..core.executor import run_plink2, run_plink
from ..core.utils import get_optimal_threads, get_available_memory, count_file_lines
from ..core.get_references import FetcherLDRegions

from pathlib import Path
from typing import Optional, Union

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class SampleQC:

    def __init__(
        self, 
        input_path: Path, 
        input_name: str, 
        output_path: Path, 
        output_name: str, 
        high_ld_regions_file: Path, 
        build: str = '38'
    ) -> None:
        
        """Initialize SampleQC class for quality control of genetic data.
        
        This class handles quality control procedures for genetic data files in PLINK binary format
        (bed, bim, fam). It sets up the directory structure and validates input files.
        
        Parameters
        ----------
        input_path : Path
            Directory path containing the input PLINK files
        input_name : str
            Base name of the input PLINK files (without extension)
        output_path : Path
            Directory path where output files will be saved
        output_name : str
            Base name for output files (without extension)
        high_ld_regions_file : Path
            Path to file containing high LD regions. If not found, will be fetched from package
        built : str, optional
            Genome build version, either '37' or '38' (default='38')
        
        Raises
        ------
        TypeError
            If input types are incorrect
        ValueError
            If genome build version is not '37' or '38'
        FileNotFoundError
            If input paths or required PLINK files are not found
        
        Attributes
        ----------
        renamed_snps : bool
            Flag indicating if SNPs should be renamed
        hh_to_missing : bool
            Flag indicating if heterozygous haploid genotypes should be set to missing
        pruned_file : None
            Placeholder for pruned file path
        results_dir : Path
            Directory for all QC results
        fails_dir : Path
            Directory for failed samples
        clean_dir : Path
            Directory for clean files
        plots_dir : Path
            Directory for QC plots
        """

        if not isinstance(input_path, Path) or not isinstance(output_path, Path):
            raise TypeError("input_path and output_path should be of type Path")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str")
        if not isinstance(high_ld_regions_file, Path):
            raise TypeError("high_ld_regions_file should be of type Path")
        
        if not isinstance(build, str):
            raise TypeError("built should be of type str")
        if build not in ['37', '38']:
            raise ValueError("built should be either 37 or 38")
        
        if not input_path.exists() or not output_path.exists():
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not (input_path / f"{input_name}.bed").exists():
            raise FileNotFoundError(".bed file not found")
        if not (input_path / f"{input_name}.fam").exists():
            raise FileNotFoundError(".fam file not found")
        if not (input_path / f"{input_name}.bim").exists():
            raise FileNotFoundError(".bim file not found")
        
        if not high_ld_regions_file.is_file():
            logger.info(f"High LD file not found at {high_ld_regions_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions(build=build)
            ld_fetcher.get_ld_regions()

            ld_regions = ld_fetcher.ld_regions
            if ld_regions is None:
                raise ValueError("Failed to fetch high LD regions file")
            logger.info(f"High LD file fetched from the package and saved at {ld_regions}")
        else:
            logger.info(f"High LD file found at {high_ld_regions_file}")
            ld_regions = high_ld_regions_file
        
        self.input_path  = Path(input_path)
        self.output_path = Path(output_path)
        self.input_name  = input_name
        self.output_name = output_name
        self.high_ld_regions_file = ld_regions

        self.processed_files = None
        self.pruned_file = None

        # create results folder
        self.results_dir = self.output_path / 'sample_qc_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # create fails folder
        self.fails_dir = self.results_dir / 'fail_samples'
        self.fails_dir.mkdir(parents=True, exist_ok=True)

        # create clean files folder
        self.clean_dir = self.results_dir / 'clean_files'
        self.clean_dir.mkdir(parents=True, exist_ok=True)
        
        # create figures folder
        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

    def execute_preprocessing(
        self, 
        rename: bool = True, 
        hh_to_missing: bool = True
    ) -> None:
        """Executes the SNP ID renaming process and/or Convert haploid genotypes to missing valuesusing PLINK2.
        
        This method renames SNP IDs in the PLINK binary files to a standardized format of 'chr:pos:a1:a2'.
        The renaming is performed using PLINK2's --set-all-var-ids parameter.
        This method uses PLINK's --set-hh-missing flag to convert haploid genotypes to 
        missing values in the genotype data. This is often useful for quality control 
        of genetic data, particularly for variants on sex chromosomes.
        
        Parameters
        ----------
        rename (bool, optional): Flag to control whether SNP renaming should be performed. 
            Defaults to True.
        hh_to_missing (bool, optional): Flag to control whether haploid genotypes should be converted to missing values.
            Defaults to True.

        Returns
        -------
            None

        Raises
        ------
            TypeError: If rename parameter is not a boolean.
            TypeError: If hh_to_missing parameter is not a boolean.

        Notes
        -----
            - The renamed files will be saved with '-renamed' suffix
            - Thread count is optimized based on available CPU cores
            - The new SNP ID format will be: chromosome:position:allele1:allele2
            - Sets self.renamed_snps to True if renaming is performed
        """

        if not isinstance(rename, bool):
            raise TypeError("rename must be a boolean")
        if not isinstance(hh_to_missing, bool):
            raise TypeError("hh_to_missing must be a boolean")
        
        max_threads = get_optimal_threads()
        memory = get_available_memory()
        
        if not rename and not hh_to_missing:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Skipping renaming of SNPs in the study data")
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Skipping conversion of haploid genotypes to missing values")
            return
        elif rename and not hh_to_missing:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Renaming SNPs in the study data to the format chr_pos_a1_a2")
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Converting haploid genotypes to missing values in the study data")
            
            plink2_cmd = [
            '--bfile', str(self.input_path / self.input_name),
            '--set-all-var-ids', '@:#:$r:$a',
            '--threads', str(max_threads),
            '--memory', str(memory),
            '--make-bed',
            '--out', str(self.input_path / (self.input_name + '-processed'))
            ]

        elif not rename and hh_to_missing:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Skipping renaming of SNPs in the study data")
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Converting haploid genotypes to missing values in the study data")
            
            plink2_cmd = [
            '--bfile', str(self.input_path / self.input_name),
            '--set-invalid-haploid-missing',
            '--threads', str(max_threads),
            '--memory', str(memory),
            '--make-bed',
            '--out', str(self.input_path / (self.input_name + '-processed'))
            ]

        else:
            logger.info(f"STEP: Rename SNPs. `rename` set to {rename}. Renaming SNPs in the study data to the format chr_pos_a1_a2")
            logger.info(f"STEP: Convert haploid genotypes to missing values. `hh_to_missing` set to {hh_to_missing}. Converting haploid genotypes to missing values in the study data")
            
            plink2_cmd = [
            '--bfile', str(self.input_path / self.input_name),
            '--set-all-var-ids', '@:#:$r:$a',
            '--set-invalid-haploid-missing',
            '--threads', str(max_threads),
            '--memory', str(memory),
            '--make-bed',
            '--out', str(self.input_path / (self.input_name + '-processed'))
            ]
        

        # Execute PLINK2 command
        run_plink2(plink2_cmd)

        self.processed_files = self.input_path / (self.input_name + '-processed')

        return
    

    def execute_ld_pruning(
        self, 
        ind_pair: list = [50, 5, 0.2]
    ) -> None:
        """Execute LD (Linkage Disequilibrium) pruning on genetic data using PLINK.
        
        This method performs LD pruning in three steps:
        1. Excludes complex/high LD regions
        2. Identifies SNPs for pruning using indep-pairwise test
        3. Creates final pruned dataset

        Parameters
        ----------
        ind_pair : list, optional
            List of three elements for LD pruning parameters:
            - Window size (int): Number of SNPs to analyze in each window
            - Step size (int): Number of SNPs to shift window at each step
            - r² threshold (float): Correlation coefficient threshold for pruning
            Default is [50, 5, 0.2]
        
        Raises
        ------
        TypeError
            If ind_pair is not a list
            If first two elements of ind_pair are not integers
            If third element of ind_pair is not float
        ValueError
            If ind_pair does not contain exactly three elements
            If window size or step size is not positive
            If r² threshold is not between 0 and 1
        FileNotFoundError
            If required pruning input file is not found
        
        Notes
        -----
        - Uses available CPU cores (leaving 2 cores free) and 2/3 of available memory
        - Creates intermediate and final files with suffixes:
          * '-LDpruning' (prune.in/prune.out variant lists only, no bed/bim/fam)
          * '-LDpruned'
        - Updates self.pruned_file with path to final pruned dataset
        """
        
        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        if len(ind_pair) != 3:
            raise ValueError("ind_pair must have exactly three elements")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (window size and step size)")
        if ind_pair[0] <= 0 or ind_pair[1] <= 0:
            raise ValueError("Window size and step size must be positive integers")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        if not (0 < ind_pair[2] <= 1):
            raise ValueError("The r^2 threshold must be a float between 0 and 1")

        logger.info("STEP: LD pruning")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        if self.processed_files is not None:
            ld_input = self.processed_files
        else:
            ld_input = self.input_name

        prune_in_file = (self.results_dir / (self.input_name + '-LDpruning')).with_suffix('.prune.in')

        # Compute the LD-pruning variant list, excluding high-LD regions in
        # the same pass instead of pre-materializing a filtered bed/bim/fam
        # just to read it back in immediately after.
        run_plink2([
            '--bfile', str(self.input_path / ld_input),
            '--exclude', str(self.high_ld_regions_file),
            '--indep-pairwise', str(ind_pair[0]), str(ind_pair[1]), str(ind_pair[2]),
            '--memory', str(memory),
            '--threads', str(max_threads),
            '--out', str(self.results_dir / (self.input_name + '-LDpruning'))
        ])

        # Extract the pruned SNPs directly from the original input. The
        # high-LD exclusion must be re-applied here too, since prune_in_file
        # only lists variants that already passed it in the step above.
        run_plink2([
            '--bfile', str(self.input_path / ld_input),
            '--exclude', str(self.high_ld_regions_file),
            '--extract', str(prune_in_file),
            '--make-bed',
            '--out', str(self.results_dir / (self.input_name + '-LDpruned')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])

        self.pruned_file = self.results_dir / (self.input_name + '-LDpruned')

        return
    
    def execute_miss_genotype(self) -> None:
        """Execute missing genotype analysis using PLINK to generate sample missingness statistics.
        
        This method generates genome-wide missingness statistics for all samples in the dataset.
        The statistics are used later to identify samples with high missingness rates during
        the quality control process.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If the output .smiss file is not generated
        
        Notes
        -----
        This function creates one file:
        - {input_name}-missing.smiss: Contains missingness statistics for all samples
        
        The method automatically optimizes thread count and memory usage based on available
        system resources (uses max(CPU cores - 2, 1) threads and 2/3 of available memory).
        """

        logger.info(f"STEP: Missing genotype check.")

        if self.processed_files is not None:
            missingness_input = self.processed_files
        else:
            missingness_input = self.input_name

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        # PLINK2 command: run missingness across file genome-wide
        run_plink2([
            '--bfile', str(self.input_path / missingness_input),
            '--missing',
            '--memory', str(memory),
            '--threads', str(max_threads),
            '--out', str(self.results_dir / (self.input_name + '-missing'))
        ])

        self.call_rate_miss = (self.results_dir / (self.input_name+'-missing')).with_suffix('.smiss')
        if not self.call_rate_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.call_rate_miss}")

        return
    
    def execute_sex_check(
        self, 
        sex_check: list = [0.2, 0.8]
    ) -> None:
        """Execute sex check using PLINK to identify potential sex discrepancies in genetic data.
        
        This method performs sex check analysis by:
        1. Running PLINK's --check-sex command on pruned data
        2. Calculating missingness rates for X chromosome SNPs (filtered inline, no intermediate file)
        
        Parameters
        ----------
        sex_check : list of float, default=[0.2, 0.8]
            List containing two float values that define the F-statistic boundaries for sex determination.
            First value is the lower bound (max-female-xf), second is the upper bound (min-male-xf).
            Samples with F-statistics below the first value are called female, above the second value are called male.
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If sex_check is not a list or if its elements are not floats
        ValueError
            If sex_check doesn't contain exactly 2 elements
        
        Notes
        -----
        The method creates the following output files:
        - {output_name}-sexcheck.sexcheck : Contains sex check results
        - {output_name}-xchr-missing.smiss : X chromosome missingness data
        
        The number of threads used is automatically determined based on available CPU cores,
        using max(available cores - 2, 1) or falling back to half of logical cores if CPU count
        cannot be determined.
        """

        if not isinstance(sex_check, list):
            raise TypeError("sex_check should be a list")
        if len(sex_check) != 2:
            raise ValueError("sex_check must have two elements")
        if not all(isinstance(i, float) for i in sex_check):
            raise TypeError("All elements in sex_check must be floats")
        
        logger.info(f"STEP: Check discordant sex information.")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        # Check sex
        run_plink2([
            '--bfile', str(self.pruned_file),
            '--check-sex', f'max-female-xf={sex_check[0]}', f'min-male-xf={sex_check[1]}',
            '--threads', str(max_threads),
            '--memory', str(memory),
            '--out', str(self.results_dir / (self.output_name + '-sexcheck'))
        ])
        # time.sleep(2)

        # X chromosome missingness, filtering to chr 23 in the same pass
        # instead of materializing a separate -xchr bed/bim/fam first.
        run_plink2([
            '--bfile', str(self.pruned_file),
            '--chr', '23',
            '--threads', str(max_threads),
            '--missing',
            '--out', str(self.results_dir / (self.output_name + '-xchr-missing'))
        ])
        # time.sleep(2)

        self.sexcheck_miss = self.results_dir / (self.output_name + '-sexcheck.sexcheck')
        self.xchr_miss = self.results_dir / (self.output_name + '-xchr-missing.smiss')

        return

    def execute_heterozygosity_rate(
        self, 
        maf: float = 0.01    
    ) -> None:
        """Executes heterozygosity rate analysis on genetic data using PLINK.

        This method performs a series of PLINK commands to analyze heterozygosity rates in genetic data,
        separating SNPs based on minor allele frequency (MAF) threshold and computing heterozygosity
        for both groups.

        Parameters
        ----------
        maf : float, optional
            Minor allele frequency threshold used to split SNPs into two groups.
            Must be between 0 and 0.5. Default is 0.01.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If maf is not a float
        ValueError
            If maf is not between 0 and 0.5
        FileNotFoundError
            If any of the expected output files are not created

        Notes
        -----
        The method:
        1. Extracts autosomal SNPs
        2. Splits SNPs based on MAF threshold
        3. Computes missingness
        4. Converts to PED/MAP format
        5. Computes heterozygosity for both MAF groups

        The computation uses optimized threading based on available CPU cores and memory.
        """

        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf <= 0 or maf >= 0.5:
            raise ValueError("maf should be between 0 and 0.5")

        logger.info(f"STEP: Heterozygosity rate check. `maf` set to {maf}")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        # Extract autosomal SNPs
        run_plink2([
            '--bfile', str(self.pruned_file),
            '--autosome',
            '--memory', str(memory),
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22'))
        ])
        # time.sleep(5)

        # Extract SNPs with MAF greater than threshold
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22')),
            '--maf', str(maf),
            '--make-bed',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafgreater')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        # Extract SNPs with MAF less than threshold
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22')),
            '--exclude', str((self.results_dir / (self.output_name + '-chr1-22-mafgreater')).with_suffix('.bim')),
            '--make-bed',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafless')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        # Get missingness for MAF greater group
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22-mafgreater')),
            '--missing',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafgreater-missing')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        # Get missingness for MAF less group
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22-mafless')),
            '--missing',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafless-missing')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        # Compute heterozygosity for MAF greater group
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22-mafgreater')),
            '--het',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafgreater')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        # Compute heterozygosity for MAF less group
        run_plink2([
            '--bfile', str(self.results_dir / (self.output_name + '-chr1-22-mafless')),
            '--het',
            '--out', str(self.results_dir / (self.output_name + '-chr1-22-mafless')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])
        # time.sleep(2)

        self.maf_greater_het= self.results_dir / (self.output_name+'-chr1-22-mafgreater.het')
        if not self.maf_greater_het.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_greater_het}")
        self.maf_less_het   = self.results_dir / (self.output_name+'-chr1-22-mafless.het')
        if not self.maf_less_het.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_less_het}")
        self.maf_greater_smiss= self.results_dir / (self.output_name+'-chr1-22-mafgreater-missing.smiss')
        if not self.maf_greater_smiss.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_greater_smiss}")
        self.maf_less_smiss   = self.results_dir / (self.output_name+'-chr1-22-mafless-missing.smiss')
        if not self.maf_less_smiss.exists():
            raise FileNotFoundError(f"Missing file: {self.maf_less_smiss}")

        return

    def execute_ibd(self) -> None:
        """
        Execute Identity by Descent (IBD) analysis using PLINK.

        This method performs duplicate and relatedness checks using IBD analysis. It runs two PLINK commands:
        1. Generates genome-wide IBD estimates
        2. Calculates missing genotype rates

        The method uses optimal thread count based on available CPU cores and validates input/output files.

        Returns:
        --------
            None

        Raises:
        -------
            FileNotFoundError: If required input pruned file is missing or if expected output files are not generated

        Required instance attributes:
            pruned_file: Path to pruned PLINK binary file
            results_dir: Directory path for output files
            output_name: Base name for output files
            ibd_miss: Path to missing genotype rate file (set by method)
            genome: Path to IBD estimates file (set by method)
        """

        logger.info("STEP: Duplicates and relatedness check with IBD")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        if self.pruned_file is None or not all([
            (self.pruned_file.with_suffix('.bed')).exists(),
            (self.pruned_file.with_suffix('.bim')).exists(),
            (self.pruned_file.with_suffix('.fam')).exists()
        ]):
            raise FileNotFoundError(f"Missing PLINK files: {self.pruned_file}")

        # Compute IBD
        run_plink([
            '--bfile', str(self.pruned_file),
            '--genome',
            '--out', str(self.results_dir / (self.output_name + '-ibd')),
            '--threads', str(max_threads),
            '--memory', str(memory)
        ])

        # Calculate missing genotype rates
        run_plink2([
            '--bfile', str(self.pruned_file),
            '--allow-no-sex',
            '--missing',
            '--out', str(self.results_dir / (self.output_name + '-ibd-missing')),
            '--threads', str(max_threads),
            '--memory', str(memory)
        ])

        self.ibd_miss = self.results_dir / (self.output_name+'-ibd-missing.smiss')
        if not self.ibd_miss.exists():
            raise FileNotFoundError(f"Missing file: {self.ibd_miss}")
        self.genome = self.results_dir / (self.output_name+'-ibd.genome')
        if not self.genome.exists():
            raise FileNotFoundError(f"Missing file: {self.genome}")

        return

    def execute_kinship(self, kinship: float = 0.354) -> None:
        """Execute kinship analysis to identify and handle sample relatedness.

        This method performs kinship analysis using PLINK2 to identify duplicate samples and related individuals.
        It first computes a kinship coefficient matrix for all samples and then prunes samples based on the
        specified kinship threshold.
        
        Parameters
        ----------
        kinship : float, optional
            The kinship coefficient threshold used to identify related samples. Must be between 0 and 1.
            Samples with kinship coefficients above this threshold will be marked for removal.
            Default is 0.354 (equivalent to first-degree relatives).
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If kinship parameter is not a float.
        ValueError
            If kinship parameter is not between 0 and 1.
        FileNotFoundError
            If the expected output file from PLINK2 is not created.
        
        Notes
        -----
        - Uses PLINK2 to compute kinship coefficients and perform sample pruning
        - Automatically determines optimal thread count and memory usage based on system resources
        - Creates output files with kinship coefficient matrix and list of samples to be removed
        - Updates self.kinship_miss with path to file containing samples to be removed
        """

        if not isinstance(kinship, float):
            raise TypeError("kinship should be a float")
        if kinship < 0 or kinship >1:
            raise ValueError("kinship should be between 0 and 1")

        logger.info(f"STEP: Duplicates and relatedness check with Kinship. `kinship` set to {kinship}")

        if self.processed_files is not None:
            kinship_input = self.processed_files
        else:
            kinship_input = self.input_name

        max_threads = get_optimal_threads()
        memory = get_available_memory()
        
        # Compute kinship-coefficient matrix for all samples
        run_plink2([
            '--bfile', str(self.input_path / kinship_input),
            '--make-king', 'triangle', 'bin',
            '--out', str(self.results_dir / (self.output_name + '-kinship-coefficient-matrix')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])

        # Prune for Monozygotic Twins OR Duplicates
        run_plink2([
            '--bfile', str(self.input_path / kinship_input),
            '--king-cutoff', str(self.results_dir / (self.output_name + '-kinship-coefficient-matrix')), str(kinship),
            '--out', str(self.results_dir / (self.output_name + '-kinship-pruned-duplicates')),
            '--memory', str(memory),
            '--threads', str(max_threads)
        ])

        self.kinship_miss = (self.results_dir / (self.output_name+'-kinship-pruned-duplicates')).with_suffix('.king.cutoff.out.id')

        # Check if the file exists
        if not self.kinship_miss.exists():
            raise FileNotFoundError(f"Expected file {self.kinship_miss} was not created. Ensure the PLINK2 command executed successfully.")

        return
    
    def execute_duplicate_relatedness(self, kinship: float = 0.354, use_kinship: bool = True) -> None:
        """
        Execute duplicate and relatedness analysis on the genotype data.
        This method performs either IBD (Identity by Descent) or KING kinship coefficient
        analysis to identify duplicate samples and related individuals in the dataset.
        
        Parameters
        ----------
        kinship : float, optional
            The KING kinship coefficient threshold for identifying related samples.
            Default is 0.354, which corresponds to duplicates/MZ twins.
        use_kinship : bool, optional
            If True, uses KING algorithm for relatedness analysis.
            If False, uses traditional IBD analysis.
            Default is True.
        
        Returns
        -------
        None

        Raises
        ------
        TypeError
            If kinship is not a float or use_kinship is not a boolean.
        
        Notes
        -----
        The method will store the analysis type (KING or IBD) in the use_kinship attribute.
        """

        if not isinstance(use_kinship, bool):
            raise TypeError("use_kinship must be a boolean")
        if not isinstance(kinship, float):
            raise TypeError("kinship must be a float")

        logger.info("STEP: Duplicates and relatedness check")

        if use_kinship:
            self.execute_kinship(kinship)
        else:
            self.execute_ibd()

        self.use_kinship = use_kinship

        return
    
    def _analyze_ibd_failures(self, ibd_threshold: float, chunk_size: int = 100000) -> pd.DataFrame:
        """Helper method to analyze IBD failures in chunks.
        
        Parameters
        ----------
        ibd_threshold : float
            PI_HAT threshold for identifying related samples
        chunk_size : int
            Number of rows to process at a time
            
        Returns
        -------
        pd.DataFrame
            DataFrame with failed samples (FID, IID, Failure columns)
        """
        imiss_path = self.results_dir / (self.output_name + '-ibd-missing.smiss')
        genome_path = self.results_dir / (self.output_name + '-ibd.genome')

        # Load missingness data
        df_imiss = pd.read_csv(imiss_path, sep=r'\s+', engine='python')
        df_imiss.columns = [col.lstrip('#') for col in df_imiss.columns]

        # Process genome file in chunks
        duplicates = []
        for chunk in pd.read_csv(
            genome_path,
            usecols=['FID1', 'IID1', 'FID2', 'IID2', 'PI_HAT'],
            sep=r'\s+',
            engine='python',
            chunksize=chunk_size,
        ):
            filtered_chunk = chunk[chunk['PI_HAT'] > ibd_threshold]
            if not filtered_chunk.empty:
                duplicates.append(filtered_chunk)

        if not duplicates:
            return pd.DataFrame(columns=['FID', 'IID', 'Failure'])

        df_dup = pd.concat(duplicates, ignore_index=True)

        # Merge with missingness
        imiss_related1 = pd.merge(
            df_dup[['FID1', 'IID1']],
            df_imiss[['FID', 'IID', 'F_MISS']],
            left_on=['FID1', 'IID1'],
            right_on=['FID', 'IID'],
        ).rename(columns={'F_MISS': 'F_MISS_1'})

        imiss_related2 = pd.merge(
            df_dup[['FID2', 'IID2']],
            df_imiss[['FID', 'IID', 'F_MISS']],
            left_on=['FID2', 'IID2'],
            right_on=['FID', 'IID'],
        ).rename(columns={'F_MISS': 'F_MISS_2'})

        # Decide which samples to remove (keep one with lower missingness)
        to_remove = pd.concat(
            [
                imiss_related1[['FID1', 'IID1', 'F_MISS_1']],
                imiss_related2[['FID2', 'IID2', 'F_MISS_2']],
            ],
            axis=1,
        )

        to_remove['FID'], to_remove['IID'] = np.where(
            to_remove['F_MISS_1'] > to_remove['F_MISS_2'],
            (to_remove['FID1'], to_remove['IID1']),
            (to_remove['FID2'], to_remove['IID2']),
        )

        fail_duplicates = to_remove[['FID', 'IID']].drop_duplicates().reset_index(drop=True)
        fail_duplicates['Failure'] = 'Duplicates and relatedness (IBD)'

        return fail_duplicates
    
    def _analyze_heterozygosity_failures(self, het_file: Path, std_deviation_het: float, maf_het: float, maf_direction: str) -> pd.DataFrame:
        """Helper method to analyze heterozygosity failures.
        
        Parameters
        ----------
        het_file : Path
            Path to the .het file
        std_deviation_het : float
            Number of standard deviations for outlier detection
        maf_het : float
            MAF threshold used
        maf_direction : str
            Direction of MAF comparison ('>' or '<')
            
        Returns
        -------
        pd.DataFrame
            DataFrame with failed samples (FID, IID, Failure columns)
        """
        # Load heterozygosity data
        df_het = pd.read_csv(het_file, sep=r'\s+', engine='python')
        df_het.columns = [col.lstrip('#') for col in df_het.columns]
        df_het["HET_RATE"] = 1 - (df_het["O(HOM)"] / df_het["OBS_CT"])

        # Calculate mean and standard deviation
        mean_percent = df_het['HET_RATE'].mean()
        sd_percent = df_het['HET_RATE'].std()

        # Identify outliers
        mask_plus = df_het['HET_RATE'] > mean_percent + std_deviation_het * sd_percent
        mask_minus = df_het['HET_RATE'] < mean_percent - std_deviation_het * sd_percent

        # Filter failed samples
        fail_het = df_het[mask_plus | mask_minus][['FID', 'IID']].copy()
        fail_het['Failure'] = f'Heterozygosity rate (MAF {maf_direction} {maf_het})'

        return fail_het
    
    def get_fail_samples(self, call_rate_thres: float, std_deviation_het: float, maf_het: float, ibd_threshold: float) -> pd.DataFrame:
        """Identify and compile all samples that failed quality control checks.
        
        This method aggregates samples that failed various QC checks including call rate,
        sex check, heterozygosity rate, and duplicate/relatedness checks. It uses chunked
        reading for memory efficiency with large datasets and provides a summary of failures.
        
        Parameters
        ----------
        call_rate_thres : float
            Call rate threshold (F_MISS). Samples with missing rate above this value fail.
            Recommended range: 0.02 to 0.1.
        std_deviation_het : float
            Number of standard deviations from mean heterozygosity rate for outlier detection.
            Typical value: 3.
        maf_het : float
            Minor allele frequency threshold used in heterozygosity analysis.
            Typical value: 0.01.
        ibd_threshold : float
            PI_HAT threshold for IBD-based relatedness detection (only used if use_kinship=False).
            Typical values: >0.185 for 2nd degree relatives, >0.5 for 1st degree.
        
        Returns
        -------
        pd.DataFrame
            Summary DataFrame with columns:
            - Failure: Type of QC failure
            - count: Number of samples failing each check
            Includes additional rows for duplicated sample IDs and totals.
        
        Raises
        ------
        FileNotFoundError
            If any required input file from previous QC steps is missing.
        
        Notes
        -----
        The method performs the following:
        
        - Reads large files in 10,000-row chunks to manage memory
        - Identifies samples failing multiple checks (reported as duplicates)
        - Saves two output files:
        
          * fail_samples.txt: List of unique samples to remove (FID, IID)
          * fail_summary.txt: Summary statistics of failures by type
        
        - Uses helper methods _analyze_heterozygosity_failures() and _analyze_ibd_failures()
          for complex analyses
        
        The analysis respects the use_kinship attribute to determine whether to use
        KING-based or IBD-based relatedness detection.
        """

        # Check if required files exist
        required_files = [
            self.call_rate_miss,
            self.sexcheck_miss,
            self.xchr_miss,
            self.maf_greater_het,
            self.maf_less_het        ]

        if not self.use_kinship:
            required_files.append(self.results_dir / (self.output_name + '-ibd-missing.smiss'))
            required_files.append(self.results_dir / (self.output_name + '-ibd.genome'))
        else:
            required_files.append(self.kinship_miss)

        for file in required_files:
            if not file.exists():
                raise FileNotFoundError(f"Required file not found: {file}")
            
        # Load and filter call rate failures in chunks
        fail_call_rate_chunks = []
        for chunk in pd.read_csv(
            self.call_rate_miss,
            sep=r'\s+',
            engine='python',
            chunksize=10000
        ):
            chunk.columns = [col.lstrip('#') for col in chunk.columns]
            failed_chunk = chunk[chunk['F_MISS'] > call_rate_thres][['FID', 'IID']].copy()
            if not failed_chunk.empty:
                fail_call_rate_chunks.append(failed_chunk)
        
        if fail_call_rate_chunks:
            fail_call_rate = pd.concat(fail_call_rate_chunks, ignore_index=True)
            fail_call_rate['Failure'] = 'Call rate'
        else:
            fail_call_rate = pd.DataFrame(columns=['FID', 'IID', 'Failure'])

        # Load and filter sex check failures in chunks
        fail_sexcheck_chunks = []
        for chunk in pd.read_csv(
            self.results_dir / (self.output_name + '-sexcheck.sexcheck'),
            sep=r'\s+',
            engine='python',
            chunksize=10000
        ):
            chunk.columns = [col.lstrip('#') for col in chunk.columns]
            failed_chunk = chunk[chunk['STATUS'] == 'PROBLEM'][['FID', 'IID']].copy()
            if not failed_chunk.empty:
                fail_sexcheck_chunks.append(failed_chunk)
        
        if fail_sexcheck_chunks:
            fail_sexcheck = pd.concat(fail_sexcheck_chunks, ignore_index=True)
            fail_sexcheck['Failure'] = 'Sex check'
        else:
            fail_sexcheck = pd.DataFrame(columns=['FID', 'IID', 'Failure'])

        # Heterozygosity failures for MAF greater than threshold
        fail_het_greater = self._analyze_heterozygosity_failures(
            het_file=self.maf_greater_het,
            std_deviation_het=std_deviation_het,
            maf_het=maf_het,
            maf_direction='>'
        )

        # Heterozygosity failures for MAF less than threshold
        fail_het_less = self._analyze_heterozygosity_failures(
            het_file=self.maf_less_het,
            std_deviation_het=std_deviation_het,
            maf_het=maf_het,
            maf_direction='<'
        )

        # Duplicates and relatedness check
        if self.use_kinship:
            # Load kinship-based duplicates/related samples
            df_duplicates = pd.read_csv(
                self.kinship_miss,
                sep=r'\s+',
                engine='python'
            )
            df_duplicates.columns = ['FID', 'IID']
            fail_duplicates = df_duplicates[['FID', 'IID']].copy()
            fail_duplicates['Failure'] = 'Duplicates and relatedness (Kinship)'
            
            del df_duplicates
        else:
            # Use IBD analysis (complex chunked processing)
            fail_duplicates = self._analyze_ibd_failures(ibd_threshold=ibd_threshold)

        # Merge all failure DataFrames
        fails = [fail_call_rate, fail_sexcheck, fail_het_greater, fail_het_less, fail_duplicates]
        df_all_fails = pd.concat(fails, axis=0, ignore_index=True)

        # Create summary statistics
        summary = df_all_fails['Failure'].value_counts().reset_index()
        num_dup = df_all_fails.duplicated(subset=['FID', 'IID']).sum()

        # Remove duplicate sample IDs (samples failing multiple checks)
        df_unique_fails = df_all_fails.drop_duplicates(subset=['FID', 'IID'])

        # Save to file
        df_unique_fails[['FID', 'IID']].to_csv(self.fails_dir / 'fail_samples.txt', index=False, sep='\t')

        # Add summary rows
        totals = summary.select_dtypes(include="number").sum() - num_dup
        dups_row = pd.DataFrame({'Failure': ['Duplicated Sample IDs'], 'count': [-num_dup]})
        total_row = pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})
        summary = pd.concat([summary, dups_row, total_row], ignore_index=True)

        logger.info(f"Total samples failing QC: {len(df_unique_fails)}")
        logger.info(f"Samples failing multiple checks: {num_dup}")

        summary.to_csv(self.results_dir / 'fail_summary.txt', index=False, sep='\t')

        return summary
    
    def execute_drop_samples(self) -> None:
        """
        Execute the removal of samples that failed quality control checks using PLINK.
        This method performs the following steps:
        1. Determines the appropriate binary file name based on previous processing steps
        2. Reads the fail_samples.txt file containing samples to be removed
        3. Executes PLINK command to create new binary files excluding failed samples
        
        Raises:
        -------
            FileNotFoundError: If the fail_samples.txt file is not found in the fails directory

        Returns:
        --------
            None

        Notes:
        ------
            - The output files will be created with suffix '-clean-samples'
            - The method preserves allele order during the operation
            - Input files must be in PLINK binary format (.bed, .bim, .fam)
        """
        
        logger.info("STEP: Drop samples that failed quality control checks")

        if self.processed_files is not None:
            binary_name = self.processed_files
        else:
            binary_name = self.input_name

        logger.info(f"Binary file name: {binary_name}")

        # Drop samples
        fail_samples_file = self.fails_dir / 'fail_samples.txt'
        if not fail_samples_file.exists():
            raise FileNotFoundError(f"Required file {fail_samples_file} not found. Ensure the fail_samples.txt file is generated before executing this step.")

        # Execute PLINK2 command
        run_plink2([
            '--bfile', str(self.input_path / binary_name),
            '--remove', str(fail_samples_file),
            '--make-bed',
            '--out', str(self.clean_dir / (self.output_name))
        ])

        return
    
    def execute_sample_qc_pipeline(self, sample_params: dict) -> None:
        """
        Execute the complete sample quality control pipeline.
        
        This method runs all sample QC steps in the correct order with proper
        memory management and logging. It encapsulates the entire workflow
        that was previously handled in the main script.
        
        Parameters
        ----------
        sample_params : dict
            Dictionary containing all sample QC parameters with keys:
            - 'rename_snp': bool, whether to rename SNPs
            - 'hh_to_missing': bool, whether to convert haploid to missing
            - 'ind_pair': list, LD pruning parameters [window, step, r2]
            - 'mind': float, missing genotype rate threshold
            - 'sex_check': list, sex check F-statistic thresholds
            - 'maf': float, minor allele frequency for heterozygosity
            - 'kinship': float, kinship coefficient threshold
            - 'use_kinship': bool, whether to use KING vs IBD
            - 'het_deviation': float, standard deviations for het filtering
            - 'ibd_threshold': float, IBD threshold for relatedness
        
        Returns
        -------
        None
        
        Notes
        -----
        The pipeline executes steps in this order:
        1. Rename SNPs (optional)
        2. Convert haploid to missing (optional)
        3. LD pruning
        4. Missing genotype analysis
        5. Sex check
        6. Heterozygosity rate analysis
        7. Duplicate/relatedness analysis
        8. Identify failed samples
        9. Remove failed samples
        10. Clean up temporary files
        
        Memory usage is monitored after each step and garbage collection
        is performed to prevent memory issues with large datasets.
        """
        
        sample_qc_steps = {
            'pre-processing'        : (self.execute_preprocessing, {"rename": sample_params['rename_snp'], "hh_to_missing": sample_params['hh_to_missing']}),
            'ld_pruning'            : (self.execute_ld_pruning, {"ind_pair": sample_params['ind_pair']}),
            'miss_genotype'         : (self.execute_miss_genotype, {}),
            'sex_check'             : (self.execute_sex_check, {"sex_check": sample_params['sex_check']}),
            'heterozygosity'        : (self.execute_heterozygosity_rate, {"maf": sample_params['maf']}),
            'duplicates_relatedness': (self.execute_duplicate_relatedness, {"kinship": sample_params['kinship'], "use_kinship": sample_params['use_kinship']}),
            'get_fail_samples'      : (self.get_fail_samples, {"call_rate_thres": sample_params['mind'], "std_deviation_het": sample_params['het_deviation'], "maf_het": sample_params['maf'], "ibd_threshold": sample_params['ibd_threshold']}),
            'drop_fail_samples'     : (self.execute_drop_samples, {}),
        }

        step_description = {
            'pre-processing'        : 'Rename SNPs to chr:pos:ref:alt and solve hh warnings by setting to missing',
            'ld_pruning'            : 'Perform LD pruning',
            'miss_genotype'         : 'Get samples with high missing rate',
            'sex_check'             : 'Get samples with discordant sex information',
            'heterozygosity'        : 'Get samples with high heterozygosity rate',
            'duplicates_relatedness': 'Get samples with high relatedness rate or duplicates',
            'get_fail_samples'      : 'Get samples that failed quality control',
            'drop_fail_samples'     : 'Drop samples that failed quality control',
            'clean_input_files'     : 'Clean input folder',
            'clean_results_files'   : 'Clean results folder',
        }

        logger.info("Starting Sample Quality Control Pipeline")
        
        for name, (func, params) in sample_qc_steps.items():
            print(f"\033[1m{step_description[name]}.\033[0m")
            
            try:
                func(**params)
                logger.info(f"Successfully completed step: {name}")
            except Exception as e:
                logger.error(f"Error in step '{name}': {str(e)}")
                raise

            # Memory management
            # time.sleep(3)  # to avoid overwhelming the system with too many operations at once
            gc.collect()  # clear memory after each step

            mem = psutil.virtual_memory()
            logger.info(f"Memory usage after {name}: {mem.percent}%")
            
        logger.info("Sample Quality Control Pipeline completed successfully")

        return

class SampleQCReport:

    def __init__(self, output_path: Path) -> None:
        self.output_path = output_path
        pass

    def report_sample_qc(self,
            call_rate_smiss: Path,
            sexcheck_miss: Path,
            xchr_miss: Path,
            maf_greater_het: Path,
            maf_less_het: Path,
            maf_greater_smiss: Path,
            maf_less_smiss: Path,
            genome: Optional[Path]= None,
            generate_ibd_report: bool = False,
            f_coeff_thresholds: list = [0.2, 0.8],
            call_rate_thres: float = 0.2, 
            std_deviation_het: float=3, 
            maf_het: float=0.01, 
            ibd_threshold: Optional[float] = None
    ) -> None:
        """Generate comprehensive quality control visualization reports for sample-level QC.
        
        This method orchestrates the generation of all sample QC visualization reports by
        calling individual reporting methods for call rate, sex check, heterozygosity, and
        optionally IBD analysis. It generates plots and visualizations without performing
        any QC execution or sample filtering.
        
        Parameters
        ----------
        call_rate_smiss : Path
            Path to the sample missingness file (.smiss) from PLINK --missing command
        sexcheck_miss : Path
            Path to the sex check results file (.sexcheck) from PLINK --check-sex
        xchr_miss : Path
            Path to the X chromosome missingness file (.smiss)
        maf_greater_het : Path
            Path to heterozygosity file (.het) for SNPs with MAF > threshold
        maf_less_het : Path
            Path to heterozygosity file (.het) for SNPs with MAF < threshold
        maf_greater_smiss : Path
            Path to missingness file (.smiss) for SNPs with MAF > threshold
        maf_less_smiss : Path
            Path to missingness file (.smiss) for SNPs with MAF < threshold
        genome : Optional[Path], default=None
            Path to the PLINK .genome file for IBD analysis. Required if generate_ibd_report=True
        generate_ibd_report : bool, default=False
            Whether to generate IBD analysis visualization
        call_rate_thres : float, default=0.2
            Call rate threshold for visualization reference line
        std_deviation_het : float, default=3
            Number of standard deviations for heterozygosity outlier visualization
        maf_het : float, default=0.01
            Minor allele frequency threshold used in the analysis
        ibd_threshold : Optional[float], default=None
            PI_HAT threshold for IBD visualization reference line. Required if generate_ibd_report=True
        
        Returns
        -------
        None
            This method generates plots as side effects and does not return data.
        
        Raises
        ------
        ValueError
            If generate_ibd_report=True but ibd_threshold or genome is None
        
        Notes
        -----
        The method generates the following visualizations in sequence:
        1. Call rate distribution plots (histogram and scatterplots)
        2. Sex check scatter plot showing F-statistics vs X chr missingness
        3. Heterozygosity rate plots for MAF > threshold
        4. Heterozygosity rate plots for MAF < threshold
        5. IBD PI_HAT distribution histogram (if generate_ibd_report=True)
        
        All plots are saved to the output_path directory specified during initialization.
        This method is intended to be used after execute_sample_qc_pipeline() has completed
        and all required QC result files have been generated.
        """
       


        # ==========================================================================================================
        #                                             CALL RATE CHECK
        # ==========================================================================================================

        # load samples who failed call rate check
        self.report_call_rate(
            smiss_file   =call_rate_smiss,
            threshold    =call_rate_thres, 
            plots_dir    =self.output_path,
            y_axis_cap   =10
        )

        logger.info('Call rate report done')

        # ==========================================================================================================
        #                                             SEX CHECK
        # ==========================================================================================================

        self.report_sex_check(
            sex_check_filename =sexcheck_miss, 
            xchr_imiss_filename=xchr_miss,
            plots_dir          =self.output_path,
            f_coeff_thresholds=f_coeff_thresholds
        )

        logger.info('Sex check report done')

        # ==========================================================================================================
        #                                       HETETROZYGOSITY RATE CHECK
        # ==========================================================================================================

        self.report_heterozygosity_rate(
            het_filename        = maf_greater_het, 
            autosomal_filename  = maf_greater_smiss, 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '>',
            plots_dir           = self.output_path
        )

        logger.info(f'Heterozygosity rate check done for MAF > {maf_het}')

        self.report_heterozygosity_rate(
            het_filename        = maf_less_het, 
            autosomal_filename  = maf_less_smiss, 
            std_deviation_het   = std_deviation_het,
            maf                 = maf_het,
            split               = '<',
            plots_dir           = self.output_path
        )

        logger.info(f'Heterozygosity rate check done for MAF < {maf_het}')

        # ==========================================================================================================
        #                                       DUPLICATES-RELATEDNESS CHECK
        # ==========================================================================================================

        if generate_ibd_report:
            
            if ibd_threshold is not None and genome is not None:

                self.report_ibd_analysis(
                    genome       =genome, 
                    ibd_threshold=ibd_threshold
                )
                logger.info('Duplicates and relatedness check done with IBD')

            else:
                raise ValueError("To generate IBD report, please provide ibd_threshold, ibd_smiss, and genome files.")

        return
    
    def report_call_rate(self, smiss_file: Path, threshold: float, plots_dir: Path, y_axis_cap: Union[int, float] = 10, color: str = '#1B9E77', line_color: str = '#D95F02', format: str = 'png') -> None:
        """
        Generate sample call rate analysis plots.
        This method reads a PLINK-format missing rate file and creates visualization plots
        showing the distribution of missing SNPs across samples. It generates two sets of plots:
        1. Histograms showing the distribution of missing SNPs (F_MISS)
        2. Scatterplots showing different views of the call rate data
        
        Parameters
        ----------
        smiss_file : Path
            Path to the PLINK format missing rate file (.smiss or .imiss)
        threshold : float
            Call rate threshold for visualization reference line (in terms of F_MISS)
        plots_dir : Path
            Directory where plots will be saved
        y_axis_cap : Union[int, float], optional
            Maximum value for y-axis in capped histogram plots. Default is 10
        color : str, optional
            Color for the main plot elements. Default is '#1B9E77'
        line_color : str, optional
            Color for threshold lines in plots. Default is '#D95F02'
        format : str, optional
            Format for saving plots (e.g., 'png', 'pdf', 'svg'). Default is 'png'
        
        Returns
        -------
        None
            This method generates plots as side effects and does not return data.
        
        Notes
        -----
        The method generates two image files:
        - call_rate_{threshold}_histogram.<format>: Contains histogram plots
        - call_rate_{threshold}_scatterplot.<format>: Contains scatter plots
        """
        
        if not isinstance(smiss_file, Path):
            raise TypeError("smiss_file should be a Path object")
        if  not smiss_file.exists():
            raise FileNotFoundError(f"File {smiss_file} does not exist")
        if not smiss_file.is_file():
            raise FileNotFoundError(f"Path {smiss_file} is not a file")
        if not isinstance(threshold, float):
            raise TypeError("threshold should be a float object")
        if not (0 <= threshold <= 1):
            raise ValueError("threshold should be between 0 and 1")
        if not isinstance(plots_dir, Path):
            raise TypeError("plots_dir should be a Path object")
        if not isinstance(y_axis_cap, (int, float)):
            raise TypeError("y_axis_cap should be an int or float")
        if not isinstance(color, str):
            raise TypeError("color should be a string representing a color")
        if not isinstance(line_color, str):
            raise TypeError("line_color should be a string representing a color")
        if not isinstance(format, str):
            raise TypeError("format should be a string representing the file format (e.g., 'png', 'pdf')")

        # load call rate data
        df_call_rate = pd.read_csv(
            smiss_file,
            sep=r'\s+',
            engine='python'
        )
        df_call_rate.columns = [col.lstrip('#') for col in df_call_rate.columns]

        # filter samples that fail call rate
        fail_call_rate = df_call_rate[df_call_rate['F_MISS'] > threshold][['FID', 'IID']].reset_index(drop=True)
        fail_call_rate['Failure'] = 'Call rate'

        # Create the figure and subplots
        fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        # First subplot: Full histogram
        axes1[0] = sns.histplot(df_call_rate['F_MISS'], bins=30, color=color, alpha=0.7, ax=axes1[0]) # type: ignore
        axes1[0].set_title("Sample Call Rate Distribution")
        axes1[0].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes1[0].set_ylabel("Frequency")

        # Second subplot: Histogram with capped y-axis
        axes1[1] = sns.histplot(df_call_rate['F_MISS'], bins=30, color=color, alpha=0.7, ax=axes1[1]) # type: ignore
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_title("Sample Call Rate Distribution (Capped)")
        axes1[1].set_xlabel("Proportion of missing SNPs (F_MISS)")

        plt.tight_layout()
        plt.savefig(plots_dir / f"call_rate_{threshold}_histogram.{format}", dpi=400)
        plt.show(block=False)
        plt.close()

        fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

        # First subplot: capped y-axis
        axes2[0] = sns.histplot(df_call_rate['F_MISS'], bins=30, color=color, alpha=0.7, ax=axes2[0]) # type: ignore
        axes2[0].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes2[0].set_title("Sample Call Rate Distribution (Capped)")
        axes2[0].set_xlabel("Proportion of missing SNPs (F_MISS)")

        # Add a vertical line at the threshold
        axes2[0].axvline(threshold, linewidth=2, color=line_color, linestyle='dashed')

        # Second subplot: Number of samples vs F_MISS
        df_call_rate_sorted = pd.DataFrame({
            'Index': range(len(df_call_rate['F_MISS'])),
            'F_MISS': sorted(df_call_rate['F_MISS'])
        })

        axes2[1] = sns.scatterplot(
            data  =df_call_rate_sorted,
            x     ='Index',
            y     ='F_MISS',
            marker='o',
            edgecolor='none',
            color =color,
            ax    =axes2[1]
        ) 
        axes2[1].set_title("Sample Call Rate")
        axes2[1].set_xlabel(f"Number of samples")
        axes2[1].set_ylabel("F_MISS")

        # Add a vertical line at the threshold
        axes2[1].axhline(threshold, linewidth=2, color=line_color, linestyle='dashed')

        # third subplot: Number of samples vs F_MISS
        axes2[2] = sns.scatterplot(
            x      =df_call_rate['F_MISS'],
            y      =np.random.normal(size=len(df_call_rate['F_MISS'])),
            markers='o',
            s      =20,
            color =color,
        )
        axes2[2].set_title("Sample Call Rate")
        axes2[2].set_xlabel("Proportion of missing SNPs (F_MISS)")
        axes2[2].set_ylabel(f"Samples")
        axes2[2].set_yticks([])
    

        # Add a vertical line at the threshold
        axes2[2].axvline(threshold, linewidth=2, color=line_color, linestyle='dashed')

        plt.tight_layout()
        plt.savefig(plots_dir / f"call_rate_{threshold}_scatterplot.{format}", dpi=400)
        plt.show(block=False)
        plt.close()

        return
    
    def report_sex_check(self, sex_check_filename: Path, xchr_imiss_filename: Path, plots_dir: Path, f_coeff_thresholds: list = [0.2, 0.8], format: str = 'png', fig_size: tuple = (8,6)) -> None:
        """
        Creates a sex check visualization based on PLINK's sex check results.
        This function reads sex check data and X chromosome missingness data, merges them,
        and generates a scatter plot to visualize potential sex discrepancies.
        
        Parameters
        ----------
        sex_check_filename : Path
            Path to PLINK's sex check results file (typically .sexcheck file)
        xchr_imiss_filename : Path
            Path to X chromosome missingness data file (.smiss or .imiss)
        plots_dir : Path
            Directory where the plot will be saved
        f_coeff_thresholds : list, optional
            List of two F coefficient thresholds [lower, upper] for reference lines.
            Default is [0.2, 0.8]
        format : str, optional
            Format for saving the plot (e.g., 'png', 'pdf', 'svg'). Default is 'png'
        fig_size : tuple, optional
            Figure size as (width, height). Default is (8, 6)
        
        Returns
        -------
        None
            This method generates plots as side effects and does not return data.
        
        Notes
        -----
        The function creates a scatter plot with:
        - Blue hollow circles for samples with Male PEDSEX
        - Green hollow circles for samples with Female PEDSEX
        - Red filled circles for problematic samples
        - Dotted red vertical lines at the specified F coefficient thresholds
        The plot is saved as 'sex_check.{format}' in the specified plots directory.
        """

        if not isinstance(sex_check_filename, Path):
            raise TypeError("sex_check_filename should be a Path object")
        if not isinstance(xchr_imiss_filename, Path):
            raise TypeError("xchr_imiss_filename should be a Path object")
        if not isinstance(plots_dir, Path):
            raise TypeError("plots_dir should be a Path object")
        if not isinstance(format, str):
            raise TypeError("format should be a string")
        if format not in ['png', 'jpeg', 'jpg', 'svg', 'pdf', 'ps']:
            raise ValueError("format should be one of ['png', 'jpeg', 'jpg', 'svg', 'pdf', 'ps]")
        

        df_sexcheck = pd.read_csv(
            sex_check_filename,
            sep   =r'\s+',
            engine='python'
        )
        df_sexcheck.columns = [col.lstrip('#') for col in df_sexcheck.columns]

        df_xchr_smiss = pd.read_csv(
            xchr_imiss_filename,
            sep   =r'\s+',
            engine='python'
        )
        df_xchr_smiss.columns = [col.lstrip('#') for col in df_xchr_smiss.columns]

        df = pd.merge(df_sexcheck, df_xchr_smiss, on=['FID', 'IID'], how='inner')

        fail_sexcheck = df[df['STATUS'] == 'PROBLEM'][['FID', 'IID']].reset_index(drop=True)
        fail_sexcheck['Failure'] = 'Sex check'

        df['Category'] = 'General'
        df.loc[df['PEDSEX'] == 1, 'Category'] = 'Male PEDSEX'
        df.loc[df['PEDSEX'] == 2, 'Category'] = 'Female PEDSEX'

        df_problem = df[df['STATUS'] == 'PROBLEM'].reset_index(drop=True)
        df = df[df['STATUS'] != 'PROBLEM'].reset_index(drop=True)

        # Define the palette (color mapping)
        palette = {
            "Male PEDSEX"  : "blue",
            "Female PEDSEX": "green"
        }

        # Define the size mapping
        size_mapping = {
            "Male PEDSEX"  : 40,
            "Female PEDSEX": 40
        }

        # Create the Matplotlib scatter plot
        fig, ax = plt.subplots(figsize=fig_size)

        # Iterate through categories to plot each group separately
        for category, group in df.groupby("Category"):
            ax.scatter(
                group["F"], 
                group["F_MISS"], 
                edgecolors=palette[category],     # Map color # type: ignore
                facecolors='none',                # Hollow circles
                s         =size_mapping[category],# Map size # type: ignore
                label     =category               # Add label for legend
            )

        ax.scatter(
            df_problem["F"], 
            df_problem["F_MISS"], 
            color     ='red',
            s         =25,
            marker    ='o',
            label     ='Problem Status',
            edgecolors=palette['Female PEDSEX'],
        )

        ax.set_ylim(bottom=-0.01)  

        # Add vertical lines
        plt.axvline(x=f_coeff_thresholds[0], color='red', linestyle='dotted')
        plt.axvline(x=f_coeff_thresholds[1], color='red', linestyle='dotted')

        # Customize labels and legend
        plt.title("Sex Check")
        plt.xlabel("X chr inbreeding (homozygosity) estimate F")
        plt.ylabel("Proportion of missing SNPs for the X chr")
        plt.legend(title='', loc='best')
        
        plt.tight_layout()
        plt.savefig(plots_dir / f'sex_check.{format}', dpi=400)
        plt.show(block=False)
        plt.close()

        return
    
    def report_heterozygosity_rate(self, het_filename: Path, autosomal_filename: Path, std_deviation_het: Union[float, int], maf: float, split: str, plots_dir: Path, y_axis_cap: Union[float, int] = 80, format: str = 'png', scatter_fig_size: tuple = (10, 6)) -> None:
        """
        Generate heterozygosity rate visualization plots for quality control analysis.
        This function loads heterozygosity and autosomal call rate data, merges them,
        and generates visualization plots showing samples with deviant heterozygosity rates.
        
        Parameters
        ----------
        het_filename : Path
            Path to the file containing heterozygosity information (.het file)
        autosomal_filename : Path
            Path to the autosomal missingness file (.smiss or .imiss)
        std_deviation_het : Union[float, int]
            Number of standard deviations to use as threshold for identifying deviant samples
        maf : float
            Minor allele frequency threshold used in the analysis
        split : str
            Direction of MAF comparison ('>' or '<') to indicate which MAF subset is being analyzed
        plots_dir : Path
            Directory where plot files will be saved
        y_axis_cap : Union[float, int], optional
            Maximum value for y-axis in capped histogram plot. Default is 80
        format : str, optional
            Format for saving plots (e.g., 'png', 'pdf', 'svg'). Default is 'png'
        scatter_fig_size : tuple, optional
            Figure size for scatter plot as (width, height). Default is (10, 6)
        
        Returns
        -------
        None
            This method generates plots as side effects and does not return data.
        
        Notes
        -----
        The function generates two types of plots:
        1. Histograms of heterozygosity rates (both uncapped and capped)
        2. Scatter plot of heterozygosity rate vs missing SNP proportion
        Files are saved with names:
        - heterozygosity_rate_greater_{maf}_histogram.{format} (if split='>')
        - heterozygosity_rate_less_{maf}_histogram.{format} (if split='<')
        - heterozygosity_rate_greater_{maf}_scatterplot.{format} (if split='>')
        - heterozygosity_rate_less_{maf}_scatterplot.{format} (if split='<')
        """

        if not isinstance(het_filename, Path):
            raise TypeError("het_filename should be a Path object")
        if not isinstance(autosomal_filename, Path):
            raise TypeError("autosomal_filename should be a Path object")
        if not isinstance(std_deviation_het, (float, int)):
            raise TypeError("std_deviation_het should be a float or int")
        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf<0 or maf>0.5:
            raise ValueError("maf should be between 0 and 0.5")
        if not isinstance(split, str):
            raise TypeError("split should be a string")
        if not isinstance(plots_dir, Path):
            raise TypeError("plots_dir should be a Path object")
        if not isinstance(y_axis_cap, (float, int)):
            raise TypeError("y_axis_cap should be a float or int")
        if not isinstance(format, str):
            raise TypeError("format should be a string")
        if format not in ['png', 'jpeg', 'jpg', 'svg', 'pdf', 'ps']:
            raise ValueError("format should be one of ['png', 'jpeg', 'jpg', 'svg', 'pdf', 'ps]")
        
        # load samples that failed heterozygosity rate check with MAF > threshold
        df_het = pd.read_csv(
            het_filename,
            sep   =r'\s+',
            engine='python'
        )
        df_het.columns = [col.lstrip('#') for col in df_het.columns]
        df_het["HET_RATE"] = 1 - (df_het["O(HOM)"] / df_het["OBS_CT"])


        # autosomal call rate per individual
        df_autosomal = pd.read_csv(
            autosomal_filename,
            sep   =r'\s+',
            engine='python'
        )
        df_autosomal.columns = [col.lstrip('#') for col in df_autosomal.columns]

        # merge both dataframes
        df_het = pd.merge(
            df_het[['FID', 'IID', 'HET_RATE']],
            df_autosomal[['FID', 'IID', 'F_MISS']],
            on =['FID', 'IID'],
            how='inner'
        )

        mean_percent= df_het['HET_RATE'].mean()
        sd_percent  = df_het['HET_RATE'].std()

        mask_plus = df_het['HET_RATE'] > mean_percent + std_deviation_het*sd_percent
        mask_minus= df_het['HET_RATE'] < mean_percent - std_deviation_het*sd_percent

        fail_het = df_het[mask_plus | mask_minus][['FID', 'IID']].reset_index(drop=True)

        if split == '>':
            fail_het['Failure'] = 'Heterozygosity rate greater'
        else:
            fail_het['Failure'] = 'Heterozygosity rate less'

        # plots

        fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

        axes1[0] = sns.histplot(df_het['HET_RATE'], bins=30, color='green', alpha=0.7, ax=axes1[0]) # type: ignore
        axes1[0].set_title("Autosomal heterozygosity")
        axes1[0].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[0].set_ylabel("Frequency")

        axes1[1] = sns.histplot(df_het['HET_RATE'], bins=30, color='green', alpha=0.7, ax=axes1[1]) # type: ignore
        axes1[1].set_title("Autosomal heterozygosity (capped)")
        axes1[1].set_xlabel(f"% Heterozygosity MAF {split} {maf}")
        axes1[1].set_ylim(0, y_axis_cap)  # Cap y-axis
        axes1[1].set_ylabel("Frequency")

        plt.tight_layout()
        
        if split == '>':
            plt.savefig(plots_dir / f"heterozygosity_rate_greater_{maf}_histogram.{format}", dpi=400)
        else:
            plt.savefig(plots_dir / f"heterozygosity_rate_less_{maf}_histogram.{format}", dpi=400)
        
        plt.show(block=False)
        plt.close()

        df_het['Deviated'] = 'Not Excluded'
        df_het.loc[mask_plus, 'Deviated'] = f'{std_deviation_het}xSD Excluded'
        df_het.loc[mask_minus, 'Deviated']= f'{std_deviation_het}xSD Excluded'

        # Create the scatter plot
        fig= plt.figure(figsize=scatter_fig_size)
        ax = fig.add_subplot(111)

        sns.scatterplot(
            data   =df_het,
            x      ='HET_RATE',
            y      ='F_MISS',
            hue    ='Deviated',
            palette={'Not Excluded': 'blue', f'{std_deviation_het}xSD Excluded': 'red'},
            markers={'Not Excluded': 'o', f'{std_deviation_het}xSD Excluded': 'o'},
            size   ='Deviated',
            sizes  ={'Not Excluded': 20, f'{std_deviation_het}xSD Excluded': 30}
        )

        ax.set_ylim(bottom=-0.01) 

        plt.title("Autosomal heterozygosity and call rate")
        plt.xlabel(f"% Heterozygosity MAF {split} {maf}")
        plt.ylabel("Proportion of missing SNPs")
        plt.legend(title='Exclusion', loc='best')

        plt.tight_layout()
        if split == '>':
            plt.savefig(plots_dir / f"heterozygosity_rate_greater_{maf}_scatterplot.{format}", dpi=400)
        else:
            plt.savefig(plots_dir / f"heterozygosity_rate_less_{maf}_scatterplot.{format}", dpi=400)
        plt.show(block=False)
        plt.close()

        return

    def report_ibd_analysis(self, genome: Path, ibd_threshold: float = 0.185, chunk_size: int = 100000) -> None:
        """Generate visualization of IBD (Identity By Descent) analysis results.
        
        This method processes IBD analysis results and creates a histogram showing the
        distribution of PI_HAT values for related sample pairs. The visualization includes
        a reference line at the specified threshold to help assess relatedness in the dataset.
        
        Parameters
        ----------
        genome : Path
            Path to the PLINK .genome file containing pairwise IBD estimates
        ibd_threshold : float, default=0.185
            The PI_HAT threshold for the reference line in the plot.
            Typical values: >0.98 for duplicates, >0.5 for first-degree relatives,
            >0.185 for second-degree relatives.
        chunk_size : int, default=100000
            Number of rows to process at a time when reading the genome file.
        
        Returns
        -------
        None
            This method generates a plot as a side effect and does not return data.
        
        Raises
        ------
        TypeError
            If ibd_threshold is not a float or chunk_size is not an integer.
        FileNotFoundError
            If the genome file is not found.
        
        Notes
        -----
        The method creates a histogram visualization showing:
        - Distribution of PI_HAT values for sample pairs with PI_HAT > 0.1
        - A vertical line indicating the ibd_threshold
        The plot is saved as 'ibd_pihat_distribution.png' in the output directory.
        """
        
        if not isinstance(ibd_threshold, float):
            raise TypeError("ibd_threshold should be a float")
        if not isinstance(chunk_size, int):
            raise TypeError("chunk_size should be an integer")

        # File paths
        if not genome.exists():
            raise FileNotFoundError(f"Missing file: {genome}")

        # Initialize dataframe for duplicates
        filtered_chunks = []

        # Process the .genome file in chunks
        for chunk in pd.read_csv(
            genome,
            usecols  =['FID1', 'IID1', 'FID2', 'IID2', 'PI_HAT'],
            sep      =r'\s+',
            engine   ='python',
            chunksize=chunk_size,
        ):
            # Filter rows with PI_HAT > 0.1 to visualize the full distribution of related samples
            # (lower than ibd_threshold to show context around the cutoff)
            filtered_chunk = chunk[chunk['PI_HAT'] > 0.1]
            if not filtered_chunk.empty:
                filtered_chunks.append(filtered_chunk)


        # Concatenate all filtered chunks
        df_chunks = pd.concat(filtered_chunks, ignore_index=True)

        # Generate histogram of PI_HAT values
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.histplot(data=df_chunks, x='PI_HAT', bins=50, color='#1B9E77', alpha=0.7, ax=ax)
        ax.axvline(ibd_threshold, color='#D95F02', linestyle='--', linewidth=2, label=f'Threshold ({ibd_threshold})')
        ax.set_title('Distribution of PI_HAT Values for Related Samples')
        ax.set_xlabel('PI_HAT (Proportion IBD)')
        ax.set_ylabel('Frequency')
        ax.legend()
        plt.tight_layout()
        plt.savefig(self.output_path / 'ibd_pihat_distribution.png', dpi=600)
        plt.show(block=False)
        plt.close()
       
        return 

class SampleQCCleanUp:

    def __init__(self, output_path: Path, input_path: Path) -> None:
        self.output_path = output_path
        self.input_path = input_path

    def clean_input_files(self) -> None:

        """Remove intermediate files from input directory.

        This method deletes temporary files created during preprocessing steps:
        - Files ending with 'processed' (.bed, .bim, .fam)

        Returns
        -------
        None

        Notes
        -----
        Only removes files if they exist. No error is raised if files are not found.
        """

        logger.info("Cleaning up processed files from input directory")

        extensions = ['.bed', '.bim', '.fam']
        
        # Remove processed files
        for ext in extensions:
            for file in self.input_path.glob(f'*-processed{ext}'):
                if file.exists():
                    file.unlink()
                    logger.info(f"Deleted: {file}")
        
        return
    
    def clean_results_files(self) -> None:

        """Remove intermediate files from output directory.

        This method deletes temporary files created during sample QC steps:
        - Files ending with '.bed', '.bim', '.fam', '.vmiss', '.smiss', '.nosex', '.sexcheck', '.het', '.genome'
        - Files ending with 'prune.in', 'prune.out'
        - Files ending with 'king.cutoff.out.id', 'king.cutoff.in.id', 'king.id', 'king.bin'

        Returns
        -------
        None

        Notes
        -----
        Only removes files if they exist. No error is raised if files are not found.
        """

        logger.info("Cleaning up intermediate files from output directory")

        # Remove intermediate files
        extensions = ['.bed', '.bim', '.fam', '.vmiss', '.smiss', '.nosex', '.sexcheck', '.het', '.genome']
        for ext in extensions:
            for file in self.output_path.glob(f'*{ext}'):
                if file.exists():
                    file.unlink()
                    logger.info(f"Deleted: {file}")
        
        # Remove prune files
        prune_patterns = ['*.prune.in', '*.prune.out']
        for pattern in prune_patterns:
            for file in self.output_path.glob(pattern):
                if file.exists():
                    file.unlink()
                    logger.info(f"Deleted: {file}")
        
        # Remove kinship files
        kinship_patterns = ['*.king.cutoff.out.id', '*.king.cutoff.in.id', '*.king.id', '*.king.bin']
        for pattern in kinship_patterns:
            for file in self.output_path.glob(pattern):
                if file.exists():
                    file.unlink()
                    logger.info(f"Deleted: {file}")

        return
    
    def clean_all(self) -> None:
        """Remove all intermediate files from input and output directories.

        This method calls clean_input_files() and clean_results_files()
        to remove temporary files created during preprocessing and sample QC steps.

        Returns
        -------
        None
        """

        self.clean_input_files()
        self.clean_results_files()

        return