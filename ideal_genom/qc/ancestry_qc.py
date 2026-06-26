import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from typing import Union, Literal, Optional
from scipy.spatial import distance as dist

from ..core.executor import run_plink, run_plink2
from ..core.utils import get_optimal_threads, get_available_memory, count_file_lines
from ..core.get_references import Fetcher1000Genome, FetcherLDRegions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class ReferenceGenomicMerger:

    def __init__(
        self,
        input_path: Path,
        input_name: str,
        output_path: Path,
        output_name: str,
        high_ld_regions_file: Path,
        reference_files: dict,
        build: str = '38'
    ) -> None:
        """
        Initialize ReferenceGenomicMerger class.

        This class performs reference-study data merging for genetic data by harmonizing
        study samples with reference populations.

        Parameters
        ----------
        input_path : Path
            Path to directory containing input files
        input_name : str
            Name of input file without extension
        output_path : Path
            Path to directory for output files
        output_name : str
            Name for output files without extension
        high_ld_regions_file : Path
            Path to file containing high LD regions to exclude
        reference_files : dict
            Dictionary containing paths to reference population files
        build : str, default='38'
            Genome build version ('37' or '38')

        Raises
        ------
        TypeError
            If input arguments are not of correct type
        ValueError
            If genome build version is not '37' or '38'
        FileNotFoundError
            If required input files/directories do not exist

        Attributes
        ----------
        reference_AC_GT_filtered : Path or None
            Filtered reference allele counts and genotypes
        study_AC_GT_filtered : Path or None
            Filtered study allele counts and genotypes
        pruned_reference : Path or None
            LD-pruned reference data
        pruned_study : Path or None
            LD-pruned study data
        reference_fixed_chr : Path or None
            Reference data with fixed chromosomes
        reference_fixed_pos : Path or None
            Reference data with fixed positions
        reference_flipped : Path or None
            Reference data with flipped alleles
        reference_cleaned : Path or None
            Final cleaned reference data
        """


        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_regions_file, Path):
            raise TypeError("high_ld_regions_file should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str):
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(build, str):
            raise TypeError("build should be a string")
        if build not in ['37', '38']:
            raise ValueError("build should be either '37' or '38'")
        
        if not input_path.exists():
            raise FileNotFoundError(f"input_path does not exist: {input_path}")
        if not output_path.exists():
            raise FileNotFoundError(f"output_path does not exist: {output_path}")
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

        self.input_path: Path = input_path
        self.input_name: str = input_name
        self.output_path: Path = output_path
        self.output_name: str = output_name
        self.build: str = build
        self.high_ld_regions_file: Path = ld_regions
        # Convert reference_files values to Path objects if they are strings
        self.reference_files: dict = {k: Path(v) if isinstance(v, str) else v for k, v in reference_files.items()}
        
        # Validate that all reference files exist
        self._validate_reference_files()

        self.reference_AC_GT_filtered: Optional[Path] = None
        self.study_AC_GT_filtered: Optional[Path] = None
        self.pruned_reference: Optional[Path] = None
        self.pruned_study: Optional[Path] = None
        self.reference_fixed_chr: Optional[Path] = None
        self.reference_fixed_pos: Optional[Path] = None
        self.reference_flipped: Optional[Path] = None
        self.reference_cleaned: Optional[Path] = None

    def execute_rename_snpid(self) -> None:
        """
        Execute the SNP ID renaming process using PLINK2.
        
        This method renames SNP IDs in the PLINK binary files to a standardized format 
        of 'chr:pos:a1:a2' using PLINK2's --set-all-var-ids parameter.
        
        Returns
        -------
        None

        Notes
        -----
        - The renamed files will be saved with '-renamed' suffix
        - Thread count is optimized based on available CPU cores
        - The new SNP ID format will be: chromosome:position:allele1:allele2
        """

        logger.info("STEP: Renaming SNP IDs in the study data using PLINK2")

        max_threads = get_optimal_threads()

        # Execute PLINK2 command
        run_plink2([
            '--bfile', str(self.input_path / self.input_name),
            '--set-all-var-ids', '@:#:$r:$a',
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.output_path / (self.input_name + '-renamed'))
        ])

        return

    def execute_filter_prob_snps(self) -> None:
        """
        Filter problematic SNPs (A/T and C/G) from both study and reference data.
        
        This method performs the following operations:
        1. Identifies and filters A/T and C/G SNPs from study data
        2. Identifies and filters A/T and C/G SNPs from reference data
        3. Creates new PLINK binary files excluding the identified problematic SNPs
        4. Uses optimal available CPU threads and memory
        
        Returns
        -------
        None

        Notes
        -----
        - Creates filtered SNP list files in the output directory
        - Creates new PLINK binary files (.bed, .bim, .fam) in the output directory
        - Sets self.reference_AC_GT_filtered and self.study_AC_GT_filtered paths
        - Logs progress and statistics of filtering operations
        - Requires valid PLINK binary files for both study and reference data
        """

        logger.info("STEP: Filtering A->T and C->G SNPs from study and reference data.")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        # Find and filter problematic SNPs in both datasets
        filtered_study = self._filter_non_AT_or_GC_snps(
            target_bim=self.output_path / f"{self.input_name}-renamed.bim",
            output_filename=self.input_name
        )
        filtered_reference = self._filter_non_AT_or_GC_snps(
            target_bim=self.reference_files['bim'],
            output_filename=self.reference_files['bim'].stem
        )

        # Set output paths
        self.reference_AC_GT_filtered = self.output_path / f"{self.reference_files['bim'].stem}-no_ac_gt_snps"
        self.study_AC_GT_filtered = self.output_path / f"{self.input_name}-no_ac_gt_snps"

        # Log filtering statistics
        study_count = count_file_lines(filtered_study)
        ref_count = count_file_lines(filtered_reference)
        logger.info(f"STEP: Filtering problematic SNPs from the study data: {study_count} SNPs filtered")
        logger.info(f"STEP: Filtering problematic SNPs from the reference data: {ref_count} SNPs filtered")


        # Make sure the reference bim path is valid and extract the base filename
        if not self.reference_files.get('bim') or not isinstance(self.reference_files['bim'], Path):
            raise ValueError("reference_files dictionary must contain a valid 'bim' Path")
        
        reference_base = self.reference_files['bim'].with_suffix('')
        
        # Execute PLINK2 command: generate cleaned study data files
        run_plink2([
            '--bfile', str(self.output_path / (self.input_name + '-renamed')),
            '--chr', '1-22',
            '--exclude', str(filtered_study),
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.study_AC_GT_filtered),
            '--memory', str(int(memory))
        ])

        # Execute PLINK2 command: generate cleaned reference data files
        run_plink2([
            '--bfile', str(reference_base),
            '--max-alleles', '2',
            '--chr', '1-22',
            '--exclude', str(filtered_reference),
            '--allow-extra-chr',
            '--memory', str(int(memory)),
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.reference_AC_GT_filtered)
        ])

        return
    
    def execute_ld_pruning(self, ind_pair:list) -> None:
        """
        Execute linkage disequilibrium (LD) pruning on study and reference data.
        
        This method performs LD-based pruning using PLINK to remove highly correlated SNPs 
        from both study and reference datasets. The pruning is done using a sliding window 
        approach where SNPs are removed based on their pairwise correlation (r²).
        
        Parameters
        ----------
        ind_pair : list
            A list containing three elements:
            
            - ind_pair[0] (int): Window size in SNPs  
            - ind_pair[1] (int): Number of SNPs to shift the window at each step  
            - ind_pair[2] (float): r² threshold for pruning
        
        Raises
        ------
        TypeError
            If ind_pair is not a list.
        TypeError
            If first two elements of ind_pair are not integers.
        TypeError
            If third element of ind_pair is not a float.
        
        Returns
        -------
        None
        
        Notes
        -----
        - Uses PLINK's `--indep-pairwise` command for pruning.
        - Excludes high LD regions specified in `self.high_ld_regions_file`.
        - Creates pruned datasets for both study and reference data.
        - Updates `self.pruned_reference` and `self.pruned_study` with paths to pruned files.
        - Uses all available CPU threads except 2 for processing.
        """


        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not isinstance(ind_pair[0], int) or not isinstance(ind_pair[1], int):
            raise TypeError("The first two elements in ind_pair values should be integers (windows size and step size)")
        
        if not isinstance(ind_pair[2], float):
            raise TypeError("The third element in ind_pair should be a float (r^2 threshold)")
        
        logger.info("STEP: LD-based pruning of study and reference data")

        max_threads = get_optimal_threads()

        # Execute PLINK2 command: generates prune.in and prune.out files from study data
        run_plink2([
            '--bfile', str(self.study_AC_GT_filtered),
            '--exclude', 'range', str(self.high_ld_regions_file),
            '--keep-allele-order',
            '--indep-pairwise', str(ind_pair[0]), str(ind_pair[1]), str(ind_pair[2]),
            '--threads', str(max_threads),
            '--out', str(self.output_path / self.input_name)
        ])

        # Execute PLINK2 command: prune study data and creates a filtered binary file
        run_plink2([
            '--bfile', str(self.study_AC_GT_filtered),
            '--extract', str((self.output_path / self.input_name).with_suffix('.prune.in')),
            '--keep-allele-order',
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.output_path / (self.input_name + '-pruned'))
        ])

        # Execute PLINK2 command: generates a pruned reference data files
        run_plink2([
            '--bfile', str(self.reference_AC_GT_filtered),
            '--extract', str((self.output_path / self.input_name).with_suffix('.prune.in')),
            '--keep-allele-order',
            '--make-bed',
            '--threads', str(max_threads),
            '--out', str(self.output_path / (self.reference_files['bim'].stem + '-pruned'))
        ])

        self.pruned_reference = self.output_path / (self.reference_files['bim'].stem+'-pruned')
        self.pruned_study = self.output_path / (self.input_name+'-pruned')

        return
    
    def execute_fix_chromosome_mismatch(self) -> None:
        """
        Fix chromosome mismatch between study data and reference panel.

        This method executes PLINK commands to correct any chromosome mismatches between the study data
        and reference panel datasets. It identifies mismatches using internal methods and updates
        the chromosome assignments in the reference panel to match the study data.

        The method performs the following steps:
        1. Identifies chromosome mismatches between study and reference BIM files
        2. Creates an update file for chromosome reassignment
        3. Executes PLINK command to update chromosome assignments in reference panel

        Returns
        -------
        None

        Notes
        -----
        - Creates new PLINK binary files with updated chromosome assignments
        - The updated files are saved with '-updateChr' suffix
        
        Raises
        ------
        ValueError
            If pruned_study or pruned_reference is None, meaning execute_ld_pruning() was not called first
        """

        logger.info("STEP: Fixing chromosome mismatch between study data and reference panel")

        # Check if pruned_study and pruned_reference have been set
        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.pruned_reference is None:
            raise ValueError("pruned_reference is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")

        max_threads = get_optimal_threads()

        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.pruned_reference.with_suffix('.bim')

        to_update_chr_file = self._find_chromosome_mismatch(study_bim, reference_bim)

        self.reference_fixed_chr = self.output_path / f"{self.reference_files['bim'].stem}-updateChr"

        line_count = count_file_lines(to_update_chr_file)
        logger.info(f"STEP: Fixing chromosome mismatch between study data and reference panel: {line_count} SNPs to update")

        if line_count == 0:
            # No SNPs to update, copy the pruned reference to fixed chromosome reference
            self.reference_fixed_chr = self.pruned_reference
            logger.info("No chromosome mismatches found. Skipping chromosome update step.")
            return

        # Execute PLINK command
        run_plink([
            '--bfile', str(self.pruned_reference),
            '--update-chr', str(to_update_chr_file),
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.reference_fixed_chr)
        ])

        return
    
    def execute_fix_position_mismatch(self) -> None:
        """
        Fix position mismatches between study data and reference panel.
        
        This method updates the positions of SNPs in the reference panel to match those 
        in the study data, creating new binary PLINK files with corrected positions.
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If pruned_study or reference_fixed_chr is None
        
        Notes
        -----
        - Creates new PLINK binary files (.bed, .bim, .fam) with updated positions
        - Logs the number of SNPs being updated
        - Updates self.reference_fixed_pos with path to corrected files
        - Requires PLINK to be installed and accessible
        - The updated files are saved with '-updatePos' suffix
        """

        logger.info("STEP: Fixing position mismatch between study data and reference panel")

        max_threads = get_optimal_threads()

        # Check if pruned_study and reference_fixed_chr have been properly set
        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.reference_fixed_chr is None:
            raise ValueError("reference_fixed_chr is not set. Make sure execute_fix_chromosome_mismatch() is called before this method.")
            
        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.reference_fixed_chr.with_suffix('.bim')

        to_update_pos_file = self._find_position_mismatch(study_bim, reference_bim)

        self.reference_fixed_pos = self.output_path / f"{self.reference_files['bim'].stem}-updatePos"

        line_count = count_file_lines(to_update_pos_file)
        logger.info(f"STEP: Fixing position mismatch between study data and reference panel: {line_count} SNPs to update")

        if line_count == 0:
            # No SNPs to update, copy the fixed chromosome reference to fixed position reference
            self.reference_fixed_pos = self.reference_fixed_chr
            logger.info("No position mismatches found. Skipping position update step.")
            return

        # Execute PLINK2 command
        run_plink2([
            '--bfile', str(self.reference_fixed_chr),
            '--update-map', str(to_update_pos_file),
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.reference_fixed_pos)
        ])

        return
    
    def execute_fix_allele_flip(self) -> None:
        """
        Execute the allele flipping process between study data and reference panel.
        
        This method identifies and flips SNPs where alleles are reversed between 
        study and reference datasets, generating a new reference panel with corrected alleles.
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If pruned_study or reference_fixed_pos is None
        
        Notes
        -----
        - Creates a .toFlip file containing SNPs requiring allele flipping
        - Generates new PLINK binary files (.bed, .bim, .fam) with flipped alleles
        - Logs the number of SNPs requiring flipping
        - Updates self.reference_flipped with the path to new flipped reference files
        - The flipped files are saved with '-flipped' suffix
        """
        logger.info("STEP: Allele flipping between study data and reference panel")

        max_threads = get_optimal_threads()
            
        # Check if pruned_study and reference_fixed_pos have been properly set
        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.reference_fixed_pos is None:
            raise ValueError("reference_fixed_pos is not set. Make sure execute_fix_position_mismatch() is called before this method.")
            
        # File paths - using with_suffix for consistency and reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.reference_fixed_pos.with_suffix('.bim')

        to_flip_file = self.output_path / f"{self.reference_files['bim'].stem}.toFlip"
        self._find_allele_flip(study_bim, reference_bim, to_flip_file)

        self.reference_flipped = self.output_path / f"{self.reference_files['bim'].stem}-flipped"

        line_count = count_file_lines(to_flip_file)
        logger.info(f"STEP: Allele flipping between study data and reference panel: {line_count} SNPs to flip")

        if line_count == 0:
            # No SNPs to flip, copy the fixed position reference to flipped reference
            self.reference_flipped = self.reference_fixed_pos
            logger.info("No SNPs require allele flipping. Skipping flipping step.")
            return

        # Execute PLINK command
        run_plink([
            '--bfile', str(self.reference_fixed_pos),
            '--flip', str(to_flip_file),
            '--keep-allele-order',
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.reference_flipped)
        ])

    def execute_remove_mismatches(self) -> None:
        """
        Remove mismatched SNPs from the reference data.
        
        This method identifies and removes SNPs where alleles don't match between 
        study and reference datasets, generating a cleaned reference dataset.
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If pruned_study or pruned_reference is None
        
        Notes
        -----
        - Creates a file listing SNPs to be removed (.toRemove)
        - Generates cleaned reference files with '-cleaned' suffix
        - Logs the number of SNPs being removed
        - Updates self.reference_cleaned with path to cleaned files
        """

        logger.info("STEP: Removing mismatched SNPs from reference data")

        max_threads = get_optimal_threads()

        if self.pruned_study is None:
            raise ValueError("pruned_study is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")
        if self.pruned_reference is None:
            raise ValueError("pruned_reference is not set. Make sure execute_ld_pruning() is called before this method and completed successfully.")

        # File paths - using with_suffix instead of with_name for more reliability
        study_bim = self.pruned_study.with_suffix('.bim')
        reference_bim = self.pruned_reference.with_suffix('.bim')

        mismatches_file = self.output_path / f"{self.reference_files['bim'].stem}.toRemove"
        self._find_allele_flip(study_bim, reference_bim, mismatches_file)

        self.reference_cleaned = self.output_path / f"{self.reference_files['bim'].stem}-cleaned"

        line_count = count_file_lines(mismatches_file)
        logger.info(f"STEP: Removing mismatched SNPs from reference data: {line_count} SNPs to remove")

        if line_count == 0:
            # No mismatches to remove, copy the flipped reference to cleaned reference
            self.reference_cleaned = self.reference_flipped
            logger.info("No mismatched SNPs found. Skipping removal step.")
            return
        
        # Execute PLINK2 command
        run_plink2([
            '--bfile', str(self.reference_flipped),
            '--exclude', str(mismatches_file),
            '--keep-allele-order',
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.reference_cleaned)
        ])

        return
    
    def execute_merge_data(self) -> None:
        """
        Merge study and reference data using PLINK.
        
        This method merges the pruned study data with the cleaned reference data 
        using PLINK's --bmerge functionality.
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If reference_cleaned is None
        
        Notes
        -----
        - Creates merged PLINK binary files (.bed, .bim, .fam) in the output directory
        - The merged files are saved with '-merged' suffix
        - Automatically determines optimal number of threads based on available CPU cores
        """

        logger.info("STEP: Merging study and reference data")

        max_threads = get_optimal_threads()

        if self.reference_cleaned is None:
            raise ValueError("reference_cleaned is not set. Make sure execute_remove_mismatches() is called before this method and completed successfully.")

        # Execute PLINK command
        run_plink([
            '--bfile', str(self.pruned_study),
            '--bmerge',
            str(self.reference_cleaned.with_suffix('.bed')),
            str(self.reference_cleaned.with_suffix('.bim')),
            str(self.reference_cleaned.with_suffix('.fam')),
            '--keep-allele-order',
            '--threads', str(max_threads),
            '--make-bed',
            '--out', str(self.output_path / (self.output_name + '-merged'))
        ])

        return

    def _filter_non_AT_or_GC_snps(self, target_bim: Path, output_filename: str) -> Path:
        """
        Filter SNPs that are not A/T or G/C variants from a PLINK BIM file.
        This method reads a BIM file and identifies SNPs that are either A/T or G/C variants.
        These variants are known as strand-ambiguous SNPs because their complementary alleles 
        are the same as their original alleles, making it impossible to determine the correct 
        strand without additional information.
        
        Parameters
        ----------
        target_bim : Path
            Path to the input BIM file containing SNP information
        output_filename : str
            Base name for the output file (without extension)
        
        Returns
        -------
        Path
            Path to the output file containing filtered SNP IDs with .ac_get_snps extension
        
        Notes
        -----
        The input BIM file should be tab-delimited with standard PLINK BIM format.
        Only columns containing SNP ID (column 2) and alleles (columns 5 and 6) are used.
        Memory-efficient implementation that processes the file in chunks to handle large datasets.
        """

        output_file = self.output_path / f"{output_filename}.ac_get_snps"
        
        # Define the problematic allele combinations
        ambiguous_alleles = {"AT", "TA", "GC", "CG"}
        
        # Process file in chunks to reduce memory usage
        chunk_size = 50000  # Process 50K rows at a time
        
        with open(output_file, 'w') as out_f:
            # Use pandas chunking to process large files efficiently
            for chunk in pd.read_csv(
                target_bim, 
                sep="\t", 
                header=None, 
                usecols=[1, 4, 5], 
                names=["SNP", "A1", "A2"],
                chunksize=chunk_size,
                dtype=str  # Ensure alleles are treated as strings
            ):
                # Vectorized operation to check for ambiguous alleles
                # Create sorted allele pairs more efficiently
                allele_pairs = chunk['A1'] + chunk['A2']
                # Sort each pair: convert to list, sort, join back
                sorted_pairs = allele_pairs.apply(lambda x: ''.join(sorted(x)))
                
                # Filter for ambiguous SNPs
                ambiguous_mask = sorted_pairs.isin(ambiguous_alleles)
                ambiguous_snps = chunk.loc[ambiguous_mask, 'SNP']
                
                # Write to file immediately to free memory
                if len(ambiguous_snps) > 0:
                    ambiguous_snps.to_csv(out_f, index=False, header=False)

        return output_file
    
    def _find_chromosome_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:
        """
        Find chromosome mismatches between study and reference BIM files.

        This function identifies SNPs where the chromosome assignment differs between
        the study dataset and the reference panel, despite having the same rsID.
        Sex chromosomes (X, Y) are excluded from the update list.

        Parameters
        ----------
        study_bim : Path
            Path to the study BIM file to check for mismatches
        reference_bim : Path
            Path to the reference BIM file to compare against

        Returns
        -------
        Path
            Path to output file containing SNPs that need chromosome updates.
            File format is tab-separated with columns: chromosome, rsID
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Find mismatches where rsID is the same but chromosome differs
        mismatch_df = reference_df.merge(study_df[["chr", "rsid"]], on="rsid", suffixes=("_ref", "_study"))
        chromosome_mismatch_df = mismatch_df[mismatch_df["chr_ref"] != mismatch_df["chr_study"]]

        # Exclude chromosomes X and Y from updates
        mismatch_df = mismatch_df[~mismatch_df["chr_study"].astype(str).isin(["X", "Y"])]

        to_update_chr_file = self.output_path / "all_phase3.toUpdateChr"

        # Save the mismatch data to a file
        chromosome_mismatch_df[["chr_study", "rsid"]].to_csv(to_update_chr_file, sep="\t", header=False, index=False)

        return to_update_chr_file
    
    def _find_position_mismatch(self, study_bim: Path, reference_bim: Path) -> Path:
        """
        Find SNPs with mismatched positions between study and reference datasets.

        This method compares the base pair positions of SNPs between a study dataset and a 
        reference dataset to identify SNPs that have different positions despite having the 
        same rsID.

        Parameters
        ----------
        study_bim : Path
            Path to the PLINK .bim file of the study dataset.
        reference_bim : Path
            Path to the PLINK .bim file of the reference dataset.

        Returns
        -------
        Path
            Path to the output file containing SNPs that need position updates.
            The output file contains two columns (rsID and new position) without headers.

        Notes
        -----
        The output file format is compatible with PLINK's --update-map command for updating
        SNP positions in the study dataset.
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary from file1 with column 2 as key and column 4 as value
        a = dict(zip(study_df['rsid'], study_df['pos_bp']))

        # Filter rows in reference_df where column 2 exists in 'a' and the values isn column 4 differ
        filtered = reference_df[reference_df['rsid'].map(a).notna() & (reference_df['pos_bp'] != reference_df['rsid'].map(a))]

        # Print the result to a file
        to_update_pos_file = self.output_path / f"{self.reference_files['bim'].stem}.toUpdatePos"
        filtered[['rsid', 'pos_bp']].to_csv(to_update_pos_file, sep="\t", header=False, index=False)

        return to_update_pos_file
    
    def _find_allele_flip(self, study_bim: Path, reference_bim: Path, output_filename: Path) -> None:
        """
        Find SNPs with allele flips between study and reference datasets.

        This method identifies Single Nucleotide Polymorphisms (SNPs) where the alleles are
        flipped between the study and reference datasets. A flip occurs when the allele
        pairs don't match in either order.

        Parameters
        ----------
        study_bim : Path
            Path to the study .bim file containing SNP information
        reference_bim : Path
            Path to the reference .bim file containing SNP information
        output_filename : Path
            Path where the list of flipped SNPs will be saved

        Returns
        -------
        None
            Writes rsids of flipped SNPs to the specified output file

        Notes
        -----
        The .bim files should be tab-separated with columns:
        chromosome, rsid, genetic_distance, base_pair_position, allele1, allele2

        The output file will contain one rsid per line for SNPs where alleles don't match
        between study and reference in either order (A1/A2 or A2/A1).
        """

        col_names = ["chr", "rsid", "pos_cm", "pos_bp", "allele1", "allele2"]
        study_df = pd.read_csv(study_bim, sep='\t', names=col_names)
        reference_df = pd.read_csv(reference_bim, sep='\t', names=col_names)

        # Create a dictionary with the composite key from file1
        a = {f"{row['chr']}{row['rsid']}{row['pos_bp']}": f"{row['allele1']}{row['allele2']}" for _, row in study_df.iterrows()}

        # Filtering the rows in file2 based on the conditions
        filtered = reference_df[
            reference_df.apply(
                lambda row: (
                    f"{row['chr']}{row['rsid']}{row['pos_bp']}" in a and 
                    a[f"{row['chr']}{row['rsid']}{row['pos_bp']}"] not in {f"{row['allele1']}{row['allele2']}", f"{row['allele2']}{row['allele1']}"}
                ), axis=1
            )
        ]

        # Save the second column of filtered rows to a file
        filtered['rsid'].to_csv(output_filename, sep="\t", header=False, index=False)

        return
    
    def execute_merging_pipeline(self, ind_pair: list = [50, 5, 0.2]) -> None:
        """
        Execute complete reference-study merging pipeline.
        
        This method runs all steps required to harmonize and merge study data
        with reference panel data:
        1. Rename SNP IDs to standardized format
        2. Filter problematic SNPs (A/T and C/G variants)
        3. Perform LD-based pruning
        4. Fix chromosome mismatches
        5. Fix position mismatches
        6. Fix allele flips
        7. Remove remaining mismatches
        8. Merge the datasets
        
        Parameters
        ----------
        ind_pair : list, default [50, 5, 0.2]
            LD pruning parameters: [window size, step size, r² threshold]
            
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If ind_pair is not a list
        
        Notes
        -----
        All intermediate files are created in the output_path directory.
        The final merged files will have '-merged' suffix.
        """
        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        logger.info("Starting Reference-Study Merging Pipeline")
        logger.info(f"LD pruning parameters: window={ind_pair[0]}, step={ind_pair[1]}, r²={ind_pair[2]}")
        
        # Execute pipeline steps sequentially
        self.execute_rename_snpid()
        self.execute_filter_prob_snps()
        self.execute_ld_pruning(ind_pair=ind_pair)
        self.execute_fix_chromosome_mismatch()
        self.execute_fix_position_mismatch()
        self.execute_fix_allele_flip()
        self.execute_remove_mismatches()
        self.execute_merge_data()
        
        logger.info("Reference-Study Merging Pipeline completed successfully")
        
        return
    
    def _validate_reference_files(self) -> None:
        """
        Validate that all reference files exist and have correct extensions.
        If validation fails, automatically fetch 1000 Genomes reference files.
        
        Raises
        ------
        ValueError
            If reference_files dictionary is empty or missing required keys after fetching
        FileNotFoundError
            If any of the reference files do not exist after fetching
        ValueError
            If files don't have expected extensions after fetching
        
        Notes
        -----
        Expected reference files: 'bed', 'bim', 'fam' (and optionally 'psam')
        If user-provided files are invalid, automatically downloads and sets up 1000 Genomes reference files
        """
        try:
            # First try to validate user-provided reference files
            self._validate_user_reference_files()
            logger.info("User-provided reference files validation successful")
        except (ValueError, FileNotFoundError, TypeError) as e:
            logger.warning(f"User-provided reference files validation failed: {e}")
            logger.info(f"Automatically fetching 1000 Genomes reference files for build {self.build}")
            
            # Import here to avoid circular imports
            from ..core.get_references import Fetcher1000Genome
            
            # Fetch 1000 Genomes reference files
            fetcher = Fetcher1000Genome(build=self.build)
            fetcher.get_1000genomes()
            fetcher.get_1000genomes_binaries()

            # Update reference_files with fetched files
            self.reference_files = {
                'bim': fetcher.bim_file,
                'bed': fetcher.bed_file,
                'fam': fetcher.fam_file,
                'psam': fetcher.psam_file
            }
            
            logger.info("Successfully fetched and set 1000 Genomes reference files")
            
            # Validate the fetched files to ensure they're properly set up
            try:
                self._validate_user_reference_files()
                logger.info("Fetched reference files validation successful")
            except Exception as fetch_error:
                raise RuntimeError(f"Failed to fetch valid reference files: {fetch_error}")
    
    def _validate_user_reference_files(self) -> None:
        """
        Validate user-provided reference files without automatic fetching.
        
        Raises
        ------
        ValueError
            If reference_files dictionary is empty or missing required keys
        FileNotFoundError
            If any of the reference files do not exist
        TypeError
            If reference file paths are not Path objects
        ValueError
            If files don't have expected extensions
        """
        if not self.reference_files:
            raise ValueError("reference_files dictionary cannot be empty")
        
        # Required reference file types
        required_files = {'bed', 'bim', 'fam'}
        optional_files = {'psam'}
        
        # Check if all required files are present
        missing_files = required_files - set(self.reference_files.keys())
        if missing_files:
            raise ValueError(f"Missing required reference files: {missing_files}")
        
        # Expected file extensions
        expected_extensions = {
            'bed': '.bed',
            'bim': '.bim', 
            'fam': '.fam',
            'psam': '.psam'
        }
        
        # Validate each reference file
        for file_type, file_path in self.reference_files.items():
            if file_type not in (required_files | optional_files):
                logger.warning(f"Unknown reference file type '{file_type}', skipping validation")
                continue
                
            if not isinstance(file_path, Path):
                raise TypeError(f"Reference file '{file_type}' must be a Path object, got {type(file_path)}")
            
            # Check if file exists
            if not file_path.exists():
                raise FileNotFoundError(f"Reference {file_type} file does not exist: {file_path}")
            
            # Check if it's a file (not a directory)
            if not file_path.is_file():
                raise ValueError(f"Reference {file_type} path is not a file: {file_path}")
            
            # Check file extension
            if file_type in expected_extensions:
                expected_ext = expected_extensions[file_type]
                if not file_path.name.endswith(expected_ext):
                    logger.warning(f"Reference {file_type} file doesn't have expected extension '{expected_ext}': {file_path}")
        
        logger.info(f"Successfully validated {len(self.reference_files)} reference files")


class GenomicOutlierAnalyzer:

    def __init__(
        self, 
        input_path: Path, 
        input_name: str, 
        merged_file: Path, 
        reference_tags: Path, 
        output_path: Path, 
        output_name: str
    ) -> None:
        """
        Initialize GenomicOutlierAnalyzer object with input and output parameters.

        Parameters
        ----------
        input_path : Path
            Path to input directory containing files to process
        input_name : str
            Name of input file 
        merged_file : Path
            Path to merged genotype file
        reference_tags : Path
            Path to file containing reference population tags
        output_path : Path
            Path to output directory
        output_name : str
            Name for output files

        Raises
        ------
        TypeError
            If input arguments are not of correct type
        FileNotFoundError
            If required input files/directories do not exist

        Attributes
        ----------
        einvectors : numpy.ndarray, None
            Principal component eigenvectors, initialized as None
        eigenvalues : numpy.ndarray, None
            Principal component eigenvalues, initialized as None
        ancestry_fails : list, None
            List of samples failing ancestry QC, initialized as None
        population_tags : pandas.DataFrame, None
            DataFrame containing population reference tags, initialized as None
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(merged_file, Path):
            raise TypeError("merged_file should be a Path object")
        if not isinstance(reference_tags, Path):
            raise TypeError("reference_tags should be a Path object")
        if not isinstance(input_name, str):
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        
        if not input_path.exists():
            raise FileNotFoundError(f"input_path does not exist: {input_path}")
        if not output_path.exists():
            raise FileNotFoundError(f"output_path does not exist: {output_path}")
        if not merged_file.with_suffix('.bed').exists():
            raise FileNotFoundError(f"merged_file does not exist: {merged_file}")
        if not reference_tags.exists():
            raise FileNotFoundError(f"reference_tags does not exist: {reference_tags}")

        self.merged_file: Path = merged_file
        self.reference_tags: Path = reference_tags
        self.output_path: Path = output_path
        self.output_name: str = output_name
        self.input_path: Path = input_path
        self.input_name: str = input_name

        self.eigenvectors: Optional[Path] = None
        self.eigenvalues: Optional[Path] = None
        self.ancestry_fails: Optional[Path] = None
        self.population_tags: Optional[Path] = None

    def execute_pca(self, pca: int = 10, maf: float = 0.01) -> None:
        """
        Perform Principal Component Analysis (PCA) on the genetic data using PLINK.

        This method executes PCA on the merged genetic data file, calculating the specified
        number of principal components. It automatically determines the optimal number of
        threads and memory allocation based on system resources.

        Parameters
        ----------
        pca : int, default=10
            Number of principal components to calculate.
            Must be a positive integer.
        maf : float, default=0.01
            Minor allele frequency threshold for filtering variants.
            Must be between 0 and 0.5.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If pca is not an integer or maf is not a float
        ValueError
            If pca is not positive or maf is not between 0 and 0.5

        Notes
        -----
        The method creates two output files:
        - {output_name}-pca.eigenvec: Contains the eigenvectors (PC loadings)
        - {output_name}-pca.eigenval: Contains the eigenvalues

        The results are stored in self.einvectors and self.eigenvalues attributes.
        """

        if not isinstance(pca, int):
            raise TypeError("pca should be an integer")
        if pca <= 0:
            raise ValueError("pca should be a positive integer")
        if not isinstance(maf, float):
            raise TypeError("maf should be a float")
        if maf < 0 or maf > 0.5:
            raise ValueError("maf should be a float between 0 and 0.5")

        logger.info("STEP: Performing principal component decomposition")

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        # Execute PLINK2 command: generate PCA for reference data
        run_plink2([
            '--bfile', str(self.merged_file),
            '--keep-allele-order',
            '--maf', str(maf),
            '--out', str(self.output_path / (self.output_name + '-pca')),
            '--pca', str(pca),
            '--memory', str(int(memory)),
            '--threads', str(max_threads)
        ])

        self.eigenvectors = self.output_path / (self.output_name+'-pca.eigenvec')
        self.eigenvalues = self.output_path / (self.output_name+'-pca.eigenval')

        return
    
    def find_ancestry_outliers(
        self, 
        ref_threshold: float, 
        stu_threshold: float, 
        reference_pop: str, 
        num_pcs: int = 2, 
        fails_dir: Path = Path(), 
        distance_metric: Union[str, float] = 'infinity'
    ) -> None:
        """
        Identifies ancestry outliers in the dataset based on PCA analysis.
        This method analyzes population structure using principal component analysis (PCA) and identifies
        samples that are potential ancestry outliers based on their distance from reference populations.
        
        Parameters
        ----------
        ref_threshold : float
            Distance threshold for reference population samples
        stu_threshold : float
            Distance threshold for study population samples
        reference_pop : str
            Name of the reference population to compare against
        num_pcs : int, optional
            Number of principal components to use in the analysis (default is 2)
        fails_dir : Path, optional
            Directory path to save failed samples information (default is empty Path)
        distance_metric : str or float, optional
            Distance metric to use for outlier detection:
            - 'infinity' or 'chebyshev' → Chebyshev distance (L∞ norm)
            - numeric p >= 1 → Minkowski distance with order p (e.g., 2 for Euclidean)
            Default is 'infinity' (Chebyshev distance)
        
        Returns
        -------
        None
            Results are stored in the ancestry_fails attribute

        Raises
        ------
        TypeError
            If parameters are not of the expected type
        ValueError
            If num_pcs is not a positive integer
            If distance_metric has invalid value
        
        Notes
        -----
        The method requires:
        - A reference tags file with population information
        - An eigenvectors file from PCA analysis
        - Both files should be previously set in the class instance
        The results are saved in:
        - population_tags: CSV file with population assignments
        - ancestry_fails: List of samples identified as ancestry outliers
        
        The distance-based approach provides more robust outlier detection compared to 
        per-dimension thresholds by considering the overall multivariate distance from 
        population centroids.
        """

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be a float")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be a float")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer")
        if num_pcs <= 0:
            raise ValueError("num_pcs should be a positive integer")
        if not isinstance(fails_dir, Path):
            raise TypeError("fails_dir should be a Path object")
        
        if not fails_dir.exists():
            logger.info("STEP: Identifying ancestry outliers: `fails_dir` does not exist.")
            logger.info(f"STEP: Identifying ancestry outliers: ancestry outliers will be saved in {self.output_path}")
            fails_dir = self.output_path
        
        logger.info("STEP: Identifying ancestry outliers")

        df_tags = pd.read_csv(self.reference_tags, sep="\t", usecols=['#IID', 'SuperPop'])
        df_tags['ID'] = '0'
        df_tags = df_tags[['ID', '#IID', 'SuperPop']]
        df_tags = df_tags.rename(columns={'ID': 'ID1', '#IID': 'ID2', 'SuperPop': 'SuperPop'})

        if self.eigenvectors is None:
            raise ValueError("eigenvectors is not set. Make sure execute_pca() is called before this method and completed successfully.")

        df = pd.read_csv(self.eigenvectors, sep=r"\s+",engine='python', header=None)
        logger.info("STEP: Identifying ancestry outliers: read eigenvec file")

        df = df[[0, 1]]
        df = df.rename(columns = {0: 'ID1', 1:'ID2'})

        df['ID2'] = df['ID2'].astype(str)
        df['ID1'] = df['ID1'].astype(str)

        df = pd.merge(df, df_tags, on=['ID1', 'ID2'], how='left')
        df['SuperPop'] = df['SuperPop'].fillna('StPop', inplace=False)

        df.to_csv((self.output_path / (self.output_name + 'pop_tags.csv')), sep='\t', index=False)

        self.population_tags = self.output_path / (self.output_name + 'pop_tags.csv')

        # filter samples who are ethnicity outliers
        ancestry_fails = self._find_pca_fails(
            output_path  = fails_dir,
            df_tags      = df,
            ref_threshold= ref_threshold,
            stu_threshold= stu_threshold,
            reference_pop= reference_pop,
            num_pcs      = num_pcs,
            distance_metric= distance_metric
        )

        self.ancestry_fails = ancestry_fails

        return
    
    def execute_drop_ancestry_outliers(self, output_dir: Path = Path()) -> None:
        """
        Drop ancestry outliers from the study data by removing samples identified as ancestry outliers
        using PLINK command line tool.
        This method reads a file containing samples identified as ancestry outliers and creates new
        binary PLINK files excluding these samples.

        Parameters
        ----------
        output_dir : Path, optional
            Directory where the cleaned files will be saved. If not provided or doesn't exist,
            files will be saved in self.output_path.
        
        Returns
        -------
        None

        Raises
        ------
        TypeError
            If output_dir is not a Path object.
        
        Notes
        -----
        The method creates new PLINK binary files (.bed, .bim, .fam) with the suffix '-ancestry-cleaned'
        excluding the samples listed in self.ancestry_fails file.
        """

        logger.info("STEP: Dropping ancestry outliers from the study data")

        if not isinstance(output_dir, Path):
            raise TypeError("output_dir should be a Path object")
        
        if not output_dir.exists():
            logger.info("STEP: Dropping ancestry outliers from the study data: `output_dir` does not exist.")
            logger.info(f"STEP: Dropping ancestry outliers from the study data: ancestry outliers will be saved in {self.output_path}")
            output_dir = self.output_path

        if self.ancestry_fails is None:
            raise ValueError("ancestry_fails is not set. Make sure find_ancestry_outliers() is called before this method and completed successfully.")

        with open(self.ancestry_fails, 'r') as f:
            logger.info(f"STEP: Dropping ancestry outliers from the study data: {len(f.readlines())} samples identified as ancestry outliers")

        # Execute PLINK2 command: create cleaned binary files
        run_plink2([
            '--bfile', str(self.input_path / self.input_name),
            '--allow-no-sex',
            '--remove', str(self.ancestry_fails),
            '--make-bed',
            '--out', str(output_dir / (self.output_name))
        ])

        return
    
    def execute_outlier_analyzer_pipeline(
        self, 
        pca: int = 10, 
        maf: float = 0.01, 
        ref_threshold: float = 3.0, 
        stu_threshold: float = 3.0, 
        reference_pop: str = 'SAS', 
        num_pcs: int = 2, 
        fails_dir: Path = Path(), 
        output_dir: Path = Path(),
        distance_metric: Union[str, float] = 'infinity'
    ) -> None:
        """
        Execute the complete genomic outlier analysis pipeline.
        
        This method runs the full workflow for identifying and removing ancestry outliers:
        1. Performs PCA on merged data
        2. Identifies ancestry outliers based on distance thresholds
        3. Removes outliers from the study data
        
        Parameters
        ----------
        pca : int, default=10
            Number of principal components to calculate
        maf : float, default=0.01
            Minor allele frequency threshold for PCA
        ref_threshold : float, default=3.0
            Distance threshold for reference population samples
        stu_threshold : float, default=3.0
            Distance threshold for study population samples
        reference_pop : str, default='EUR'
            Name of the reference population to compare against
        num_pcs : int, default=2
            Number of principal components to use in outlier detection
        fails_dir : Path, optional
            Directory to save failed samples information
        output_dir : Path, optional
            Directory to save cleaned output files
        distance_metric : str or float, default='infinity'
            Distance metric for outlier detection:
            - 'infinity' or 'chebyshev' → Chebyshev distance
            - numeric p >= 1 → Minkowski distance with order p
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If parameters are not of the expected type
        ValueError
            If parameter values are invalid
        
        Notes
        -----
        This is a convenience method that orchestrates the complete pipeline.
        It calls execute_pca(), find_ancestry_outliers(), and execute_drop_ancestry_outliers()
        in sequence. Results are stored in the class attributes.
        """
        
        logger.info("Starting Genomic Outlier Analysis Pipeline")
        
        # Step 1: Perform PCA
        self.execute_pca(pca=pca, maf=maf)
        
        # Step 2: Identify ancestry outliers
        self.find_ancestry_outliers(
            ref_threshold  =ref_threshold,
            stu_threshold  =stu_threshold,
            reference_pop  =reference_pop,
            num_pcs        =num_pcs,
            fails_dir      =fails_dir,
            distance_metric=distance_metric
        )
        
        # Step 3: Remove ancestry outliers from study data
        self.execute_drop_ancestry_outliers(output_dir=output_dir)
        
        logger.info("Genomic Outlier Analysis Pipeline completed successfully")
        
        return

    def _find_pca_fails(
        self, 
        output_path: Path, 
        df_tags: pd.DataFrame, 
        ref_threshold: float, 
        stu_threshold: float, 
        reference_pop: str, 
        num_pcs: int = 2, 
        distance_metric: Union[str, float] = 'infinity'
    ) -> Path:
        """
        Identifies ancestry outliers based on PCA results using distance-based thresholds:
        one for reference population and another for study population.

        Parameters
        ----------
        output_path : Path
            Path where the output file will be saved
        df_tags : pd.DataFrame
            DataFrame containing subject IDs and population tags
        ref_threshold : float
            Distance threshold for identifying outliers relative to reference population mean
        stu_threshold : float
            Distance threshold for identifying outliers relative to study population mean
        reference_pop : str
            Reference population name as it appears in df_tags
        num_pcs : int, optional
            Number of principal components to use in the analysis (default is 2)
        distance_metric : str or float, optional
            Distance metric to use:
            - 'infinity' or 'chebyshev' → Chebyshev distance (L∞ norm)
            - numeric p >= 1 → Minkowski distance with order p (e.g., 2 for Euclidean)
            Default is 'infinity'

        Returns
        -------
        Path
            Path to the output file containing the IDs of subjects identified as ancestry outliers

        Raises
        ------
        TypeError
            If ref_threshold, stu_threshold are not numeric
            If reference_pop is not a string
            If num_pcs is not an integer
            If distance_metric is not a string or numeric value
        ValueError
            If ref_threshold, stu_threshold are not positive
            If num_pcs is less than 1
            If num_pcs is greater than available PCs in eigenvec file
            If distance_metric has invalid value

        Notes
        -----
        The method identifies outliers using a distance-based approach that considers
        the overall multivariate distance from population centroids, rather than 
        per-dimension thresholds. This provides more robust outlier detection.
        
        Subjects are identified as outliers if they exceed BOTH:
        1. The reference population distance threshold (ref_threshold)
        2. The study population distance threshold (stu_threshold)
        """

        if not isinstance(ref_threshold, (float, int)):
            raise TypeError("ref_threshold should be an integer or float value")
        if not isinstance(stu_threshold, (float, int)):
            raise TypeError("stu_threshold should be an integer or float value")
        if stu_threshold<=0:
            raise ValueError("stu_threshold should be a positive value")
        if ref_threshold<=0:
            raise ValueError("ref_threshold should be a positive value")
        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(num_pcs, int):
            raise TypeError("num_pcs should be an integer value")
        if num_pcs<1:
            raise ValueError("num_pcs should be a positive integer")

        # filters reference subjects
        mask1 = (df_tags['SuperPop']==reference_pop)
        # filters subjects from study data
        mask2 = (df_tags['SuperPop']=='StPop')

        # generates two data frames with filtered subjects
        df_ref = df_tags[mask1].reset_index(drop=True)
        df_stu = df_tags[mask2].reset_index(drop=True)

        if self.eigenvectors is None:
            raise ValueError("eigenvectors is not set. Make sure execute_pca() is called before this method and completed successfully.")

        # read .eigenvec file
        df_eigenvec = pd.read_csv(
            self.eigenvectors,
            #header=None,
            sep   =r"\s+",
            engine='python'
        )

        if num_pcs>df_eigenvec.shape[1]-2:
            raise ValueError("num_pcs should be less than or equal to the number of principal components in the .eigenvec file")
        
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:2+num_pcs]].copy()

        # renames columns for consistency
        new_col_names = []
        for k in range(2+num_pcs):
            if k<2:
                new_col_names.append(f"ID{k+1}")
            else:
                new_col_names.append(f"pc_{k-1}")
        df_eigenvec.columns = new_col_names

        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)
        df_eigenvec['ID2'] = df_eigenvec['ID2'].astype(str)

        # merge filtered subjects with its principal components
        df_ref = df_ref.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)
        df_stu = df_stu.merge(df_eigenvec, on=['ID1', 'ID2'])\
            .drop(columns=['SuperPop'], inplace=False)

        # computes mean and standard deviation by columns in reference data
        mean_ref = df_ref[df_ref.columns[2:]].mean()
        std_ref = df_ref[df_ref.columns[2:]].std()

        # computes mean and standard deviation by columns in study data
        mean_stu = df_stu[df_stu.columns[2:]].mean()
        std_stu = df_stu[df_stu.columns[2:]].std()

        # Compute distances from study samples to reference population centroid
        distances_ref = self._compute_distances(
            data_df=df_stu[df_stu.columns[2:]], 
            vec_mean=mean_ref, 
            vec_std=std_ref,
            distance=distance_metric
        )

        # Compute distances from study samples to study population centroid
        distances_stu = self._compute_distances(
            data_df=df_stu[df_stu.columns[2:]], 
            vec_mean=mean_stu, 
            vec_std=std_stu,
            distance=distance_metric
        )

        # Identify outliers based on distance thresholds
        ref_outliers_mask = distances_ref > ref_threshold
        stu_outliers_mask = distances_stu > stu_threshold

        # Subjects must be outliers relative to both reference and study populations
        combined_outliers_mask = ref_outliers_mask & stu_outliers_mask

        # Get the outlier samples
        df_outliers = df_stu[combined_outliers_mask][['ID1', 'ID2']].copy()

        ancestry_fails = output_path / (self.output_name + '_fail-ancestry-qc.txt')

        logger.info(f"STEP: Identifying ancestry outliers using {distance_metric} distance metric")
        logger.info(f"STEP: Reference threshold: {ref_threshold}, Study threshold: {stu_threshold}")
        logger.info(f"STEP: Identifying ancestry outliers: {df_outliers.shape[0]} samples identified as ancestry outliers")

        # save samples considered as ethnicity outliers
        df_outliers.to_csv(
            ancestry_fails,
            header=False,
            index =False,
            sep   ='\t'
        )

        return ancestry_fails
    
    def _compute_distances(
        self,
        data_df: Union[pd.DataFrame, np.ndarray],
        vec_mean: Union[pd.Series, np.ndarray],
        vec_std: Union[pd.Series, np.ndarray],
        distance: Union[str, float] = 'infinity'
    ) -> np.ndarray:
        """
        Compute standardized Minkowski- or Chebyshev-type distances between each
        row in data_df and the reference centroid (vec_mean).

        Parameters
        ----------
        data_df : DataFrame or ndarray
            Samples as rows, features as columns.
        vec_mean : Series or ndarray
            Reference mean vector.
        vec_std : Series or ndarray
            Reference standard deviation vector.
        distance : str or float
            Distance metric:
            - 'infinity' or 'chebyshev' → Chebyshev distance
            - numeric p >= 1 → Minkowski distance with order p

        Returns
        -------
        np.ndarray
            Vector of distances (length = number of samples).
        """

        # Convert inputs to numpy arrays
        X = np.asarray(data_df)
        mu = np.asarray(vec_mean)
        sigma = np.asarray(vec_std)

        # Standardize data (avoid division by zero)
        sigma_safe = np.where(sigma == 0, 1, sigma)
        X_std = (X - mu) / sigma_safe

        # Determine which metric to use
        if isinstance(distance, str):
            distance = distance.lower()
            if distance in ['infinity', 'chebyshev']:
                metric = 'chebyshev'
                kwargs = {}
            else:
                raise ValueError("distance must be 'infinity' or 'chebyshev' when passed as string")
        elif isinstance(distance, (float, int)):
            if distance == float('inf'):
                metric = 'chebyshev'
                kwargs = {}
            elif distance >= 1:
                metric = 'minkowski'
                kwargs = {'p': distance}
            else:
                raise ValueError("numeric distance must be >= 1 or float('inf')")
        else:
            raise TypeError("distance must be a string or numeric value")

        # Compute distances from standardized data to origin
        if metric == 'chebyshev':
            distances = dist.cdist(X_std, [np.zeros_like(mu)], metric='chebyshev').flatten()
        else:
            # For minkowski distance, pass p parameter explicitly
            p_value = kwargs.get('p', 2)
            distances = dist.cdist(X_std, [np.zeros_like(mu)], metric='minkowski', p=p_value).flatten()

        return distances


class AncestryQC:

    def __init__(
        self, 
        input_path: Path, 
        input_name: str, 
        output_path: Path, 
        output_name: str, 
        high_ld_regions_file: Path, 
        reference_files: dict = dict(), 
        recompute_merge: bool = True, 
        build: str = '38', 
        rename_snps: bool = False
    ) -> None:
        """Initialize AncestryQC class.
        
        This class performs ancestry quality control analysis on genetic data by merging it with 1000 Genomes reference data
        and running principal component analysis.

        Parameters:
        -----------
        input_path: Path 
            Path to directory containing input files
        input_name: str 
            Base name of input files (without extension) 
        output_path: Path 
            Path to directory where output files will be saved
        output_name: str 
            Base name for output files
        high_ld_regions_file: Path 
            Path to file containing high LD regions to exclude
        reference_files: dict (optional) 
            Dictionary with paths to reference files. Must contain 'bim', 'bed', 'fam' and 'psam' keys. 
            If not provided, will download 1000 Genomes reference files. Defaults to empty dict.
        recompute_merge: bool (optional): 
            Whether to recompute merge with reference even if merged files exist. Defaults to True.
        build: str (optional) 
            Genome build version, either '37' or '38'. Defaults to '38'.
        rename_snps: bool (optional): 
            Whether to rename SNPs to avoid duplicates during merge. Defaults to False.
        
        Raises
        ------
        TypeError
            If input arguments are not of expected types
        ValueError
            If build is not '37' or '38'
        FileNotFoundError
            If input_path or output_path do not exist
        
        Notes
        -----
        Creates the following directory structure under output_path:
        
        - ancestry_qc_results/
            - merging/
            - plots/ 
            - fail_samples/
            - clean_files/
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(high_ld_regions_file, Path):
            raise TypeError("high_ld_regions_file should be a Path object")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not isinstance(output_name, str):
            raise TypeError("output_name should be a string")
        if not isinstance(recompute_merge, bool):
            raise TypeError("recompute_merge should be a boolean")
        if not isinstance(build, str):
            raise TypeError("build should be a string")
        if build not in ['37', '38']:
            raise ValueError("build should be either '37' or '38'")
        if not isinstance(rename_snps, bool):
            raise TypeError("rename_snps should be a boolean")
        
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not high_ld_regions_file.is_file():
            logger.info(f"High LD file not found at {high_ld_regions_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions(build=build)
            ld_fetcher.get_ld_regions()

            ld_regions = ld_fetcher.ld_regions
            if ld_regions is None:
                raise ValueError("Failed to fetch high LD regions file")
            high_ld_regions_file = ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_regions_file}")
        
        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.output_name= output_name
        # Convert reference_files values to Path objects if they are strings
        if reference_files:
            self.reference_files = {k: Path(v) if isinstance(v, str) else v for k, v in reference_files.items()}
        else:
            self.reference_files = reference_files
        self.high_ld_regions_file = high_ld_regions_file
        self.recompute_merge = recompute_merge
        self.build = build
        self.rename_snps = rename_snps

        if not reference_files:

            logger.info(f"No reference files provided. Fetching 1000 Genomes reference data for build {self.build}")

            fetcher = Fetcher1000Genome(build=self.build)
            fetcher.get_1000genomes()
            fetcher.get_1000genomes_binaries()

            self.reference_files = {
                'bim': fetcher.bim_file,
                'bed': fetcher.bed_file,
                'fam': fetcher.fam_file,
                'psam': fetcher.psam_file
            }
        
        # Validate that all reference files exist (whether user-provided or fetched)
        self._validate_reference_files()

        self.eigenvalues: Optional[Path] = None
        self.einvectors: Optional[Path] = None
        self.ancestry_fails: Optional[Path] = None
        self.population_tags: Optional[Path] = None

        self.results_dir = self.output_path / 'ancestry_qc_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.merging_dir = self.results_dir / 'merging'
        self.merging_dir.mkdir(parents=True, exist_ok=True)

        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

        self.fail_samples_dir = self.results_dir / 'fail_samples'
        self.fail_samples_dir.mkdir(parents=True, exist_ok=True)

        self.clean_dir = self.results_dir / 'clean_files'
        self.clean_dir.mkdir(parents=True, exist_ok=True)

        pass

    def merge_reference_study(self, ind_pair: list = [50, 5, 0.2]) -> None:
        """
        Merge reference and study data by applying quality control filters and merging steps.
        This method performs a series of quality control steps to merge study data with reference data:
        1. Filters problematic SNPs
        2. Performs LD pruning
        3. Fixes chromosome mismatches
        4. Fixes position mismatches  
        5. Fixes allele flips
        6. Removes remaining mismatches
        7. Merges the datasets
        
        Parameters
        ----------
        ind_pair : list, default [50, 5, 0.2]
            Parameters for LD pruning: [window size, step size, r2 threshold]
        
        Returns
        -------
        None
        
        Notes
        -----
        If recompute_merge is False, the method will skip the merging process and expect
        merged data to already exist in the merging directory.
        
        Raises
        ------
        TypeError
            If ind_pair is not a list
        """

        if not isinstance(ind_pair, list):
            raise TypeError("ind_pair should be a list")
        
        if not self.recompute_merge:
            logger.info("STEP: Merging study and reference data: recompute_merge is set to False. Skipping merging step")
            logger.info(f"STEP: Merging study and reference data: merged data is expected to be in {self.merging_dir}")
            return

        rgm = ReferenceGenomicMerger(
            input_path= self.input_path,
            input_name= self.input_name,
            output_path= self.merging_dir, 
            output_name= self.output_name,
            high_ld_regions_file =self.high_ld_regions_file, 
            reference_files = self.reference_files,
            build= self.build
        )

        rgm.execute_merging_pipeline(ind_pair=ind_pair)

        return
    
    def _clean_merging_dir(self) -> None:
        """
        Cleans up the merging directory by removing unnecessary files.
        This method removes all files in the merging directory except:
        - Files containing '-merged' in their name
        - Log files with '.log' extension
        The cleanup helps manage disk space and removes intermediate files that are no longer needed
        after the merging process is complete.
        
        Returns:
        --------

            None
        """
        

        for file in self.merging_dir.iterdir():
            if file.is_file() and '-merged' not in file.name and file.suffix != '.log':
                file.unlink()

        return
    
    def _validate_reference_files(self) -> None:
        """
        Validate that all reference files exist and have correct extensions.
        If validation fails, automatically fetch 1000 Genomes reference files.
        
        Raises
        ------
        ValueError
            If reference_files dictionary is empty or missing required keys after fetching
        FileNotFoundError
            If any of the reference files do not exist after fetching
        ValueError
            If files don't have expected extensions after fetching
        
        Notes
        -----
        Expected reference files: 'bed', 'bim', 'fam' (and optionally 'psam')
        If user-provided files are invalid, automatically downloads and sets up 1000 Genomes reference files
        """
        try:
            # First try to validate user-provided reference files
            self._validate_user_reference_files()
            logger.info("User-provided reference files validation successful")
        except (ValueError, FileNotFoundError, TypeError) as e:
            logger.warning(f"User-provided reference files validation failed: {e}")
            logger.info(f"Automatically fetching 1000 Genomes reference files for build {self.build}")
            
            # Import here to avoid circular imports
            from ..core.get_references import Fetcher1000Genome
            
            # Fetch 1000 Genomes reference files
            fetcher = Fetcher1000Genome(build=self.build)
            fetcher.get_1000genomes()
            fetcher.get_1000genomes_binaries()

            # Update reference_files with fetched files
            self.reference_files = {
                'bim': fetcher.bim_file,
                'bed': fetcher.bed_file,
                'fam': fetcher.fam_file,
                'psam': fetcher.psam_file
            }
            
            logger.info("Successfully fetched and set 1000 Genomes reference files")
            
            # Validate the fetched files to ensure they're properly set up
            try:
                self._validate_user_reference_files()
                logger.info("Fetched reference files validation successful")
            except Exception as fetch_error:
                raise RuntimeError(f"Failed to fetch valid reference files: {fetch_error}")
    
    def _validate_user_reference_files(self) -> None:
        """
        Validate user-provided reference files without automatic fetching.
        
        Raises
        ------
        ValueError
            If reference_files dictionary is empty or missing required keys
        FileNotFoundError
            If any of the reference files do not exist
        TypeError
            If reference file paths are not Path objects
        ValueError
            If files don't have expected extensions
        """
        if not self.reference_files:
            raise ValueError("reference_files dictionary cannot be empty")
        
        # Required reference file types
        required_files = {'bed', 'bim', 'fam'}
        optional_files = {'psam'}
        
        # Check if all required files are present
        missing_files = required_files - set(self.reference_files.keys())
        if missing_files:
            raise ValueError(f"Missing required reference files: {missing_files}")
        
        # Expected file extensions
        expected_extensions = {
            'bed': '.bed',
            'bim': '.bim', 
            'fam': '.fam',
            'psam': '.psam'
        }
        
        # Validate each reference file
        for file_type, file_path in self.reference_files.items():
            if file_type not in (required_files | optional_files):
                logger.warning(f"Unknown reference file type '{file_type}', skipping validation")
                continue
                
            if not isinstance(file_path, Path):
                raise TypeError(f"Reference file '{file_type}' must be a Path object, got {type(file_path)}")
            
            # Check if file exists
            if not file_path.exists():
                raise FileNotFoundError(f"Reference {file_type} file does not exist: {file_path}")
            
            # Check if it's a file (not a directory)
            if not file_path.is_file():
                raise ValueError(f"Reference {file_type} path is not a file: {file_path}")
            
            # Check file extension
            if file_type in expected_extensions:
                expected_ext = expected_extensions[file_type]
                if not file_path.name.endswith(expected_ext):
                    logger.warning(f"Reference {file_type} file doesn't have expected extension '{expected_ext}': {file_path}")

    def execute_pca(
        self, 
        ref_population: str, 
        pca: int = 10, 
        maf: float = 0.01, 
        num_pca: int = 10, 
        ref_threshold: float = 4, 
        stu_threshold: float = 4, 
        distance_metric: Union[str, float] = 'infinity'
    ) -> None:
        """
        Performs Principal Component Analysis (PCA) on genetic data and identifies ancestry outliers.

        This method executes a complete PCA workflow including:
        1. Running the PCA analysis
        2. Identifying ancestry outliers using distance-based detection
        3. Removing identified outliers

        Parameters
        ----------
        ref_population : str
            Reference population identifier for ancestry comparison
        pca : int, optional
            Number of principal components to calculate (default=10)
        maf : float, optional
            Minor allele frequency threshold for filtering (default=0.01)
        num_pca : int, optional
            Number of principal components to use in outlier detection (default=10)
        ref_threshold : float, optional
            Distance threshold for reference population outlier detection (default=4)
        stu_threshold : float, optional
            Distance threshold for study population outlier detection (default=4)
        aspect_ratio : str or float, optional
            Aspect ratio for PCA plots (default='equal')
        distance_metric : str or float, optional
            Distance metric to use for outlier detection:
            - 'infinity' or 'chebyshev' → Chebyshev distance (L∞ norm)
            - numeric p >= 1 → Minkowski distance with order p (e.g., 2 for Euclidean)
            Default is 'infinity' (Chebyshev distance)
        explained_variance_threshold : float, optional
            Threshold for reporting significant principal components based on explained variance (default=0.01)

        Returns
        -------
        None
            Results are saved to specified output directories

        Notes
        -----
        The method uses the GenomicOutlierAnalyzer class to perform the analysis and 
        saves results in the directories specified during class initialization.
        The distance-based outlier detection provides more robust identification of
        ancestry outliers compared to per-dimension thresholds.
        """

        # Make sure the reference tag path is valid before creating the analyzer
        if 'psam' not in self.reference_files or not isinstance(self.reference_files['psam'], Path):
            raise ValueError("Reference files dictionary must contain a valid 'psam' Path")
        
        goa = GenomicOutlierAnalyzer(
            input_path= self.input_path, 
            input_name= self.input_name,
            merged_file= self.merging_dir / (self.output_name + '-merged'),
            reference_tags= self.reference_files['psam'],
            output_path= self.results_dir, 
            output_name= self.output_name
        )

        logger.info("STEP: Running PCA analysis with parameters:")
        logger.info(f"  - Reference population: {ref_population}")
        logger.info(f"  - PCs to calculate: {pca}, PCs for outlier detection: {num_pca}")
        logger.info(f"  - MAF threshold: {maf}")
        logger.info(f"  - Outlier thresholds: ref={ref_threshold}, study={stu_threshold}")
        logger.info(f"  - Distance metric: {distance_metric}")

        goa.execute_pca(pca=pca, maf=maf)
        goa.find_ancestry_outliers(
            ref_threshold=ref_threshold, 
            stu_threshold=stu_threshold, 
            reference_pop=ref_population, 
            num_pcs      =num_pca, 
            fails_dir    =self.fail_samples_dir,
            distance_metric=distance_metric
        )
        goa.execute_drop_ancestry_outliers(output_dir=self.clean_dir)

        self.eigenvectors   = goa.eigenvectors
        self.eigenvalues    = goa.eigenvalues
        self.population_tags= goa.population_tags
        self.ancestry_fails = goa.ancestry_fails

        return

    def execute_ancestry_qc_pipeline(self, ancestry_params: dict) -> None:
        """
        Execute complete ancestry QC pipeline.
        
        This method runs the full ancestry quality control workflow including:
        1. Merging reference and study data
        2. Cleaning intermediate files
        3. Running PCA analysis and outlier detection
        
        Parameters
        ----------
        ancestry_params : dict
            Dictionary containing pipeline parameters. Required keys:
            - ind_pair : list - LD pruning parameters [window, step, r2]
            - reference_pop : str - Reference population name
            - pca : int - Number of PCs to calculate
            - maf : float - MAF threshold
            - num_pcs : int - Number of PCs for outlier detection
            - ref_threshold : float - Reference outlier threshold
            - stu_threshold : float - Study outlier threshold
            - aspect_ratio : str or float - Plot aspect ratio
            - distance_metric : str or float - Distance metric for outliers
            - explained_variance_threshold : float - Variance threshold for reporting
            
        Returns
        -------
        None
        """
        
        ancestry_qc_steps = {
            'merge_study_reference'    : (self.merge_reference_study, {"ind_pair":ancestry_params['ind_pair']}),
            'delete_intermediate_files': (self._clean_merging_dir, {}),
            'pca_analysis'             : (self.execute_pca, 
                {
                    "ref_population": ancestry_params['reference_pop'],
                    "pca":ancestry_params['pca'],
                    "maf":ancestry_params['maf'],
                    "num_pca":ancestry_params['num_pcs'],
                    "ref_threshold":ancestry_params['ref_threshold'],
                    "stu_threshold":ancestry_params['stu_threshold'],
                    "distance_metric":ancestry_params['distance_metric']
                }
            ),
        }

        step_description = {
            'merge_study_reference'    : "Merging reference genome with study genome",
            'delete_intermediate_files': "Deleting intermediate files generated during merging",
            'pca_analysis'             : "Running PCA analysis to perform ancestry QC"
        }

        logger.info("=" * 70)
        logger.info("Starting Ancestry QC Pipeline")
        logger.info("=" * 70)
        
        for step_num, (name, (func, params)) in enumerate(ancestry_qc_steps.items(), 1):
            logger.info(f"\nSTEP {step_num}/{len(ancestry_qc_steps)}: {step_description[name]}")
            logger.info("-" * 70)
            try:
                func(**params)
                logger.info(f"✓ Step {step_num} completed successfully")
            except Exception as e:
                logger.error(f"✗ Step {step_num} failed: {str(e)}")
                raise
        
        logger.info("=" * 70)
        logger.info("Ancestry QC Pipeline completed successfully")
        logger.info("=" * 70)

        return


class AncestryQCReport:

    def __init__(
        self, 
        output_path: Path,  
        einvectors: Path,  
        eigenvalues: Path, 
        ancestry_fails: Path, 
        population_tags: Path
    ) -> None:
        """
        Initialize ReportAncestryCheck class for generating ancestry QC reports and visualizations.
        
        Parameters
        ----------
        output_path : Path
            Path to output directory for reports and plots
        population_tags : Path
            Path to population tags file
        einvectors : Path
            Path to eigenvectors file from PCA
        eigenvalues : Path
            Path to eigenvalues file from PCA
        ancestry_fails : Path
            Path to ancestry fails file
            
        Raises
        ------
        TypeError
            If output_path is not a Path object
            If output_name is not a string
        """
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(population_tags, Path):
            raise TypeError("population_tags should be a Path object")
        if not isinstance(einvectors, Path):
            raise TypeError("einvectors should be a Path object")
        if not isinstance(eigenvalues, Path):
            raise TypeError("eigenvalues should be a Path object")
        if not isinstance(ancestry_fails, Path):
            raise TypeError("ancestry_fails should be a Path object")
        
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not population_tags.is_file():
            raise FileNotFoundError("population_tags file does not exist")
        if not einvectors.is_file():
            raise FileNotFoundError("einvectors file does not exist")
        if not eigenvalues.is_file():
            raise FileNotFoundError("eigenvalues file does not exist")
        if not ancestry_fails.is_file():
            raise FileNotFoundError("ancestry_fails file does not exist")
        
        self.output_path: Path = output_path
        self.population_tags: Path = population_tags
        self.einvectors: Path = einvectors
        self.eigenvalues: Path = eigenvalues
        self.ancestry_fails: Path = ancestry_fails

    def report_ancestry_qc(
            self,
            reference_pop: str,
            aspect_ratio: Union[Literal['auto', 'equal'], float]='equal',
            format: str='svg'
        ) -> None:

        self.draw_pca_plot(
            reference_pop=reference_pop,
            aspect_ratio=aspect_ratio,
            exclude_outliers=False,
            plot_dir=self.output_path,
            plot_name='pca_plot',
            format=format
        )

        self.draw_pca_plot(
            reference_pop=reference_pop,
            aspect_ratio=aspect_ratio,
            exclude_outliers=True,
            plot_dir=self.output_path,
            plot_name='pca_plot_no_outliers',
            format=format
        )

        self.report_pca(

        )

        pass

    def draw_pca_plot(
        self, 
        reference_pop: str, 
        aspect_ratio: Union[Literal['auto', 'equal'], float], 
        exclude_outliers: bool = False, 
        plot_dir: Path = Path(), 
        plot_name: str = 'pca_plot', 
        format: str = 'svg'
    ) -> None:
        """
        Generate 2D and 3D PCA plots from eigenvector data and population tags.
        This method creates two PCA visualization plots:
        - A 2D scatter plot showing PC1 vs PC2 colored by super-population
        - A 3D scatter plot showing PC1 vs PC2 vs PC3 colored by super-population
        
        Parameters
        ----------
        reference_pop : str
            Reference population identifier for zoomed plots
        aspect_ratio : Union[Literal['auto', 'equal'], float]
            Aspect ratio for the plot axes. Can be 'auto', 'equal', or a numeric value
        exclude_outliers : bool, default=False
            Whether to exclude ancestry outliers from the plots
        plot_dir : Path, optional
            Directory path where plots will be saved. Defaults to current directory.
            If directory doesn't exist, plots will be saved in self.output_path
        plot_name : str, optional
            Base name for the plot files. Defaults to 'pca_plot.svg'.
            Final filenames will be prefixed with '2D-' and '3D-'
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If plot_dir is not a Path object
            If plot_name is not a string
            If reference_pop is not a string
        ValueError
            If required attributes (population_tags, einvectors, eigenvalues) are not set
        
        Notes
        -----
        Requires the following class attributes to be set:
        - self.population_tags : Path to population tags file (tab-separated)
        - self.einvectors : Path to eigenvectors file (space-separated)
        - self.eigenvalues : Path to eigenvalues file
        - self.ancestry_fails : Path to ancestry fails file (if exclude_outliers=True)
        The population tags file should contain columns 'ID1', 'ID2', and 'SuperPop'
        The eigenvectors file should contain the principal components data
        """

        logger.info("STEP: Generating PCA plots")

        if not isinstance(reference_pop, str):
            raise TypeError("reference_pop should be a string")
        if not isinstance(plot_dir, Path):
            raise TypeError("plot_dir should be a Path object")
        if not isinstance(plot_name, str):
            raise TypeError("plot_name should be a string")
        if not isinstance(exclude_outliers, bool):
            raise TypeError("exclude_outliers should be a boolean")
        
        if not plot_dir.exists():
            logger.info('STEP: Generating PCA plots: `plot_dir` does not exist.')
            logger.info(f'STEP: Generating PCA plots: pca plots will be saved in {self.output_path}')
            plot_dir = self.output_path

        if self.population_tags is None:
            raise ValueError("population_tags is not set. Please set it using the population_tags attribute or pass it during initialization.")

        # add population tags to pca output
        df_tags = pd.read_csv(self.population_tags, sep='\t')
        df_tags['ID1'] = df_tags['ID1'].astype(str)

        if self.einvectors is None:
            raise ValueError("einvectors is not set. Please set it using the einvectors attribute or pass it during initialization.")
        if self.eigenvalues is None:
            raise ValueError("eigenvalues is not set. Please set it using the eigenvalues attribute or pass it during initialization.")
        
        # load .eigenval file and calculate variance explained by the first two PCs
        df_eigenval = pd.read_csv(
            self.eigenvalues,
            header=None,
            sep   =r"\s+",
            engine='python'
        )

        total_variance = df_eigenval[0].sum()
        pc1_var = df_eigenval[0][0]
        pc2_var = df_eigenval[0][1]
        pc3_var = df_eigenval[0][2]

        pc1_var_perc = round((pc1_var / total_variance) * 100, 2)
        pc2_var_perc = round((pc2_var / total_variance) * 100, 2)
        pc3_var_perc = round((pc3_var / total_variance) * 100, 2)

        # load .eigenvec file and keep the first three principal components
        df_eigenvec = pd.read_csv(
            self.einvectors,
            #header=None,
            sep   =r"\s+",
            engine='python'
        )
        df_eigenvec = df_eigenvec[df_eigenvec.columns[:5]].copy()
        df_eigenvec.columns = ['ID1', 'ID2', 'pc_1', 'pc_2', 'pc_3']
        df_eigenvec['ID1'] = df_eigenvec['ID1'].astype(str)

        if exclude_outliers:
            # load ancestry outliers
            if self.ancestry_fails is None:
                raise ValueError("ancestry_fails is not set. Please set it using the ancestry_fails attribute or pass it during initialization.")
            logger.info("STEP: Generating PCA plots: excluding ancestry outliers")

            if self.ancestry_fails.stat().st_size == 0:
                logger.info("No ancestry outliers found; 'exclude_outliers' plot will show all samples.")
                df_outliers = pd.DataFrame(columns=['ID1', 'ID2'])
            else:
                df_outliers = pd.read_csv(self.ancestry_fails, sep=r'\s+', header=None, engine='python')
                df_outliers.columns = ['ID1', 'ID2']
                df_outliers['ID1'] = df_outliers['ID1'].astype(str)
                df_outliers['ID2'] = df_outliers['ID2'].astype(str)

            df_eigenvec = df_eigenvec.merge(df_outliers, on=['ID1', 'ID2'], how='left', indicator=True)
            df_eigenvec = df_eigenvec[df_eigenvec['_merge'] == 'left_only'].drop(columns=['_merge'])

            plot_name = f'no-outliers-{plot_name}'

        # merge to get data with tagged populations
        df = pd.merge(df_eigenvec, df_tags, on=['ID1', 'ID2'])

        # generates a 2D scatter plot
        fig, ax = plt.subplots(figsize=(10,10))
        sns.scatterplot(data=df, x='pc_1', y='pc_2', hue='SuperPop', ax=ax, marker='.', s=70)
        ax.set_aspect(aspect_ratio, adjustable='datalim')
        plt.xlabel(f'PC_1 ({pc1_var_perc}%)')
        plt.ylabel(f'PC_2 ({pc2_var_perc}%)')
        fig.tight_layout()
        fig.savefig(plot_dir / f'2D-aspect-{aspect_ratio}-{plot_name}.{format}', format=format, dpi=400)

        fig.clf()
        plt.close()

        fig3, ax3 = plt.subplots(figsize=(10,10))
        df_zoom = df[(df['SuperPop'] == 'StPop') | (df['SuperPop'] == reference_pop)].reset_index(drop=True)
        sns.scatterplot(data=df_zoom, x='pc_1', y='pc_2', hue='SuperPop', ax=ax3, marker='.', s=70)
        ax3.set_aspect(aspect_ratio, adjustable='datalim')
        plt.xlabel(f'PC_1 ({pc1_var_perc}%)')
        plt.ylabel(f'PC_2 ({pc2_var_perc}%)')
        fig3.tight_layout()
        fig3.savefig(plot_dir / f'2D-zoom-aspect-{aspect_ratio}-{plot_name}.{format}', format=format, dpi=400)

        # generates a 3D scatter plot
        fig2= plt.figure()
        ax  = fig2.add_subplot(111, projection='3d')

        grouped = df.groupby('SuperPop')
        for s, group in grouped:
            ax.scatter(
                group['pc_1'],
                group['pc_2'],
                group['pc_3'],
                label=s
            )

        ax.legend()
        plt.savefig(plot_dir / f'3D-{plot_name}.{format}', format=format, dpi=400)
        plt.close()

        return
    
    def report_pca(self, threshold: float = 0.01)-> None:
        """
        Generate PCA report including scree plot and variance explained analysis.
        
        Parameters
        ----------
        threshold : float, default=0.01
            Threshold for determining significant principal components (as fraction, e.g., 0.01 = 1%)
        
        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If eigenvalues attribute is not set
        TypeError
            If threshold is not a float
        
        Notes
        -----
        Creates two output files:
        - Scree plot with eigenvalues and cumulative variance
        - TSV file with detailed PCA statistics
        """
        
        if not isinstance(threshold, (float, int)):
            raise TypeError("threshold should be a float")
        if threshold <= 0 or threshold >= 1:
            raise ValueError("threshold should be between 0 and 1")
        
        if self.eigenvalues is None:
            raise ValueError("eigenvalues is not set. Please set it using the eigenvalues attribute or pass it during initialization.")

        eigenvalues = np.loadtxt(self.eigenvalues)

        # -------------------------------
        # 2. Compute variance explained
        # -------------------------------
        var_explained = eigenvalues / np.sum(eigenvalues)
        cum_var = np.cumsum(var_explained)

        # -------------------------------
        # 3. Scree plot and cumulative variance in subplots
        # -------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(14,5))

        # Scree plot
        axes[0].plot(range(1, len(eigenvalues)+1), eigenvalues, 'o-', color='blue')
        axes[0].set_xlabel('Principal Component')
        axes[0].set_ylabel('Eigenvalue')
        axes[0].set_title('Scree Plot (Eigenvalues)')
        axes[0].grid(True)

        # Cumulative variance plot
        axes[1].plot(range(1, len(eigenvalues)+1), cum_var, 'o-', color='green')
        axes[1].set_xlabel('Principal Component')
        axes[1].set_ylabel('Cumulative Variance Explained')
        axes[1].set_title('Cumulative Variance Explained')
        axes[1].grid(True)

        plt.tight_layout()
        plt.savefig(self.output_path / '_scree_plot.png', dpi=600)
        plt.close()

        # -------------------------------
        # 4. Determine significant PCs
        # -------------------------------
        significant_threshold = var_explained > threshold

        # -------------------------------
        # 5. Create table
        # -------------------------------
        df = pd.DataFrame({
            'PC': range(1, len(eigenvalues)+1),
            'Eigenvalue': eigenvalues,
            'Variance_Explained': var_explained,
            'Cumulative_Variance': cum_var,
            'Significant_Threshold': significant_threshold
        })

        df.to_csv(self.output_path / '_pca_report.tsv', sep='\t', index=False)

        # Optional: summary
        print("\nSummary:")
        print(f"PCs > {threshold*100}% variance: {np.sum(significant_threshold)}")

        return
    
    def _set_population_tags(
        self, 
        psam_path: Path, 
        study_fam_path: Path
    ) -> pd.DataFrame:
        """
        Sets population tags for genetic data by combining information from a PSAM file and a study FAM file.

        This method processes population information from reference data (PSAM file) and study data (FAM file), 
        combining them into a single DataFrame with consistent column naming and structure.

        Parameters
        ----------
        psam_path : Path
            Path to the PSAM file containing reference population information.
        study_fam_path : Path
            Path to the FAM file containing study individual IDs.

        Returns
        -------
        pd.DataFrame
            Combined DataFrame containing:
                - ID1: Family or group identifier (0 for reference data)
                - ID2: Individual identifier
                - SuperPop: Population tag ('StPop' for study individuals, actual population for reference data)

        Notes
        -----
        The PSAM file should contain at least '#IID' and 'SuperPop' columns.
        The FAM file should be space-separated with no header.
        """

        # Read population information from the .psam file
        df_psam = pd.read_csv(
            psam_path,
            sep='\t',
            usecols=['#IID', 'SuperPop']
        )

        # Set an ID column and rename columns for consistency
        df_psam['ID'] = 0
        df_psam = df_psam[['ID', '#IID', 'SuperPop']]
        df_psam.columns = ['ID1', 'ID2', 'SuperPop']

        # read individual IDs from the study .fam file
        df_fam = pd.read_csv(
            study_fam_path,
            sep=' ',
            header=None,
            index_col=False
        )

        # select relevant columns, assign a placeholder population tag, and rename columns
        df_fam = df_fam[df_fam.columns[:2]].copy()
        df_fam['SuperPop'] = 'StPop'
        df_fam.columns = ['ID1', 'ID2', 'SuperPop']

        # concatenate the two DataFrames to merge the information
        return pd.concat([df_fam, df_psam], axis=0)