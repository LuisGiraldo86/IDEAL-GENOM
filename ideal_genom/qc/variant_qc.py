"""
Python module to perform variant quality control
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from pathlib import Path
from typing import Optional

from ..core.executor import run_plink, run_plink2
from ..core.utils import get_available_memory, get_optimal_threads

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class VariantQC:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, output_name: str) -> None:
        """
        Initialize the VariantQC class.
        This class handles quality control for genetic variants data stored in PLINK binary format
        (.bed, .bim, .fam files).
        
        Parameters:
        -----------
        input_path: Path 
            Directory path containing input PLINK files
        input_name: str 
            Base name of input PLINK files (without extension)
        output_path: Path 
            Directory path where output files will be saved
        output_name: str 
            Base name for output files
        
        Raises:
        -------
        TypeError: 
            If input_path/output_path are not Path objects or if input_name/output_name are not strings
        FileNotFoundError: 
            If input_path/output_path don't exist or required PLINK files are not found
        
        Attributes:
        -----------
        input_path: Path 
            Path to input directory
        output_path: Path 
            Path to output directory  
        input_name: str 
            Base name of input files
        output_name: str 
            Base name for output files
        hwe_results: 
            Storage for Hardy-Weinberg equilibrium test results
        results_dir: Path 
            Directory for all QC results
        fails_dir: Path 
            Directory for failed samples
        clean_dir: Path 
            Directory for cleaned files
        plots_dir: Path 
            Directory for QC plots
        """

        if not isinstance(input_path, Path) or not isinstance(output_path, Path):
            raise TypeError("input_path and output_path should be of type Path")
        if not isinstance(input_name, str) or not isinstance(output_name, str):
            raise TypeError("input_name and output_name should be of type str")
        
        if not input_path.exists() or not output_path.exists():
            raise FileNotFoundError("input_path or output_path is not a valid path")
        if not (input_path / f"{input_name}.bed").exists():
            raise FileNotFoundError(".bed file not found")
        if not (input_path / f"{input_name}.fam").exists():
            raise FileNotFoundError(".fam file not found")
        if not (input_path / f"{input_name}.bim").exists():
            raise FileNotFoundError(".bim file not found")
        
        self.input_path = input_path
        self.output_path= output_path
        self.input_name = input_name
        self.output_name= output_name

        self.hwe_results = None

        # create results folder
        self.results_dir = self.output_path / 'variant_qc_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # create fails folder
        self.fails_dir = self.results_dir / 'fail_markers'
        self.fails_dir.mkdir(parents=True, exist_ok=True)

        # create clean files folder
        self.clean_dir = self.results_dir / 'clean_files'
        self.clean_dir.mkdir(parents=True, exist_ok=True)
        
        # create figures folder
        self.plots_dir = self.results_dir / 'plots'
        self.plots_dir.mkdir(parents=True, exist_ok=True)

    def execute_missing_data_rate(self, chr_y: int = 24) -> None:
        """
        Executes missing data rate analysis using PLINK for male and female subjects separately.
        This method performs two PLINK operations:
        1. Generates .vmiss and .smiss files for male subjects on chromosome Y
        2. Generates .vmiss and .smiss files for all subjects excluding chromosome Y
        
        Parameters
        ----------
        chr_y : int, default=24
            Chromosome Y number in the dataset. Must be between 0 and 26.
        
        Returns
        -------
        None
        
        Raises
        ------
        TypeError
            If chr_y is not an integer
        ValueError
            If chr_y is not between 0 and 26
        
        Notes
        -----
        The method uses 2/3 of available system memory for PLINK operations.
        Output files are generated in the results directory with the following naming pattern:
        - {output_name}-missing-males-only.lmiss/.imiss : For male subjects
        - {output_name}-missing-not-y.lmiss/.imiss : For non-Y chromosome data
        The results are stored in self.males_missing_data and self.females_missing_data as Path objects.
        """

        # check type for chr_y
        if not isinstance(chr_y, int):
            raise TypeError("chr_y should be of type integer.")
        
        if chr_y < 0 or chr_y > 26:
            raise ValueError("chr_y should be between 1 and 26")

        logger.info("Identifying markers with excessive missing rate...")

        memory = get_available_memory()
        threads = get_optimal_threads()

        # generates  .lmiss and .imiss files for male subjects
        run_plink2([
            '--bfile', str(self.input_path / self.input_name),
            '--missing',
            '--filter-males',
            '--chr', str(chr_y),
            '--out', str(self.results_dir / (self.output_name+'-missing-males-only')),
            '--memory', str(memory),
            '--threads', str(threads)
        ])

        # generates .lmiss and. imiss files for female subjects
        run_plink2([
            '--bfile', str(self.input_path / self.input_name),
            '--missing',
            '--not-chr', str(chr_y),
            '--out', str(self.results_dir / (self.output_name+'-missing-not-y')),
            '--memory', str(memory),
            '--threads', str(threads)
        ])

        self.males_missing_data = self.results_dir / (self.output_name+'-missing-males-only.vmiss')
        self.females_missing_data = self.results_dir / (self.output_name+'-missing-not-y.vmiss')

        return

    def execute_different_genotype_call_rate(self) -> None:
        """
        Execute test for different genotype call rates between cases and controls using PLINK.

        This method performs the following operations:
        1. Calculates available memory for PLINK execution
        2. Runs PLINK's --test-missing command to identify markers with significantly different
        missing rates between cases and controls
        3. Generates a .missing file with the results

        The method uses approximately 2/3 of available system memory for PLINK execution.

        Returns
        -------
            None

        Notes
        -----
        - Creates a .missing file in the results directory
        - Sets self.case_control_missing path attribute
        """

        logger.info("Identifying markers with different genotype call rates between cases and controls...")

        # Get the virtual memory details
        memory = get_available_memory()
        threads = get_optimal_threads()

        # generates .missing file
        run_plink([
            '--bfile', str(self.input_path / self.input_name),
            '--test-missing',
            '--out', str(self.results_dir / (self.output_name+'-case-control-missing')),
            '--memory', str(memory),
            '--threads', str(threads)
        ])

        self.case_control_missing = self.results_dir / (self.output_name+'-case-control-missing.missing')

        return
    
    def get_fail_variants(self, marker_call_rate_thres: float = 0.2, case_controls_thres: float = 1e-5) -> pd.DataFrame:
        """
        Identify and consolidate failing variants based on multiple quality control criteria.
        This method combines the results of three QC checks:
        1. Variants with high missing data rates
        2. Variants with significantly different genotype call rates between cases and controls
        3. Variants failing Hardy-Weinberg equilibrium test
        
        Parameters
        ----------
        marker_call_rate_thres : float, optional
            Threshold for failing variants based on missing data rate (default: 0.2)
        case_controls_thres : float, optional
            P-value threshold for differential missingness between cases and controls (default: 1e-5)
        hwe_threshold : float, optional
            P-value threshold for Hardy-Weinberg equilibrium test (default: 5e-8)
        
        Returns
        -------
        pd.DataFrame
            A summary DataFrame containing:
            - Counts of variants failing each QC criterion
            - Number of variants failing multiple criteria (duplicates)
            - Total number of unique failing variants
        
        Notes
        -----
        - Results are also written to a tab-separated file 'fail_markers.txt'
        - Variants failing multiple criteria are only counted once in the final output file
        """

        if not isinstance(marker_call_rate_thres, float):
            raise TypeError("marker_call_rate_thres should be of type float.")
        if not isinstance(case_controls_thres, float):
            raise TypeError("case_controls_thres should be of type float.")

        # ==========================================================================================================
        #                                             MARKERS WITH MISSING DATA 
        # ==========================================================================================================

        fail_missing_data = self.report_missing_data(
            filename_male  =self.males_missing_data, 
            filename_female=self.females_missing_data,
            threshold      =marker_call_rate_thres
        )

        # ==========================================================================================================
        #                                             MARKERS WITH DIFFERENT GENOTYPE CALL RATE
        # ==========================================================================================================

        fail_genotype = self.report_different_genotype_call_rate(
            filename =self.case_control_missing, 
            threshold=case_controls_thres, 
        )

        fails = pd.concat([fail_missing_data, fail_genotype], axis=0, ignore_index=True)

        summary = fails['Failure'].value_counts().reset_index()
        num_dup = fails.duplicated(subset=['SNP']).sum()

        totals = summary.select_dtypes(include="number").sum() - num_dup
        dups_row = pd.DataFrame({'Failure':['Duplicated SNPs'], 'count':[-num_dup]})
        total_row = pd.DataFrame({col: [totals[col] if col in totals.index else "Total"] for col in summary.columns})

        fails = fails.drop_duplicates(subset='SNP', keep='first', inplace=False)

        fails = fails.drop(columns=['Failure'], inplace=False)

        fails.to_csv(self.fails_dir / 'fail_markers.txt', header=False, index=False)

        df_report = pd.concat([summary, dups_row, total_row], ignore_index=True)
        df_report.to_csv(self.results_dir / 'variant_qc_summary.tsv', sep='\t', index=False)

        return df_report

    def execute_drop_variants(self, maf: float = 5e-8, geno: float = 0.1, hwe: float = 5e-8) -> None:
        """
        Execute variant filtering based on quality control parameters using PLINK.

        This method removes variants that fail quality control criteria including minor allele frequency (MAF),
        genotype missingness rate, and Hardy-Weinberg equilibrium (HWE) test.

        Parameters
        ----------
        maf : float, optional
            Minor allele frequency threshold. Variants with MAF below this value are removed.
            Default is 5e-8.
        geno : float, optional
            Maximum per-variant missing genotype rate. Variants with missing rate above this 
            value are removed. Default is 0.1 (10%).
        hwe : float, optional 
            Hardy-Weinberg equilibrium test p-value threshold. Variants with HWE p-value below
            this are removed. Default is 5e-8.

        Returns
        -------
        None
            Creates quality controlled PLINK binary files (.bed, .bim, .fam) in the clean directory
            with suffix '-variantQCed'.
        """

        logger.info("Removing markers failing quality control...")
        if not isinstance(maf, float):
            raise TypeError("maf should be of type float.")
        if not isinstance(geno, float):
            raise TypeError("geno should be of type float.")
        if not isinstance(hwe, float):
            raise TypeError("hwe should be of type float.") 

        # create cleaned binary files
        run_plink([
            '--bfile', str(self.input_path / self.input_name),
            '--exclude', str(self.fails_dir / 'fail_markers.txt'),
            '--autosome',
            '--maf', str(maf),
            '--hwe', str(hwe), 'midp',
            '--geno', str(geno),
            '--make-bed',
            '--out', str(self.clean_dir / (self.output_name))
        ])

        return

    def report_missing_data(self, filename_male: Path, filename_female: Path, threshold: float) -> pd.DataFrame:
        """
        Analyze and report missing data rates for male and female subjects.
        This method processes missing data information from separate files for male and female subjects,
        creates visualizations of missing data distributions, and identifies SNPs that fail the missing
        data threshold for each sex group.
        
        Parameters
        ----------
        directory : str
            Path to the directory containing the input files
        filename_male : str
            Name of the file containing missing data information for male subjects (.lmiss format)
        filename_female : str
            Name of the file containing missing data information for female subjects (.lmiss format)
        threshold : float
            Maximum allowed missing data rate (between 0 and 1)
        y_axis_cap : int, optional
            Upper limit for y-axis in histogram plots (default is 10)
        
        Returns
        -------
        pd.DataFrame
            A DataFrame containing SNPs that fail the missing data threshold for either sex,
            with columns ['SNP', 'Failure'] where 'Failure' indicates the failing category
        
        Notes
        -----
        The method generates two histogram plots saved as 'missing_data_male' and 
        'missing_data_female' showing the distribution of missing data rates for each sex.
        """
   
        if not isinstance(threshold, float):
            raise TypeError("threshold should be of type float.")
        if not isinstance(filename_male, Path) or not isinstance(filename_female, Path):
            raise TypeError("filename_male and filename_female should be of type Path.")

        # load .vmiss file for male subjects
        df_males = pd.read_csv(
            filename_male,
            sep=r"\s+"
        )
        
        ## filter male subjects
        fail_males = df_males[df_males['F_MISS']>=threshold].reset_index(drop=True)
        fail_males = fail_males[['ID']].copy()
        fail_males['Failure'] = 'Missing data rate on males'

        # load .vmiss file for female subjects
        df_females = pd.read_csv(
            filename_female,
            sep=r"\s+"
        )
        
        ## filter female subjects
        fail_females = df_females[df_females['F_MISS']>=threshold].reset_index(drop=True)
        fail_females = fail_females[['ID']].copy()
        fail_females['Failure'] = 'Missing data rate on females'

        # concatenate female and male subjects who failed QC
        fails = pd.concat([fail_females, fail_males], axis=0)

        return fails

    def report_different_genotype_call_rate(self, filename: Path, threshold: float) -> pd.DataFrame:
        """
        Reports markers with different genotype call rates based on a given threshold.
        This function reads a .missing file, filters markers with a different genotype call rate
        below the specified threshold, and returns a DataFrame containing these markers.

        Parameters:
        -----------
            directory (str): The directory where the .missing file is located.
            filename (str): The name of the .missing file.
            threshold (float): The threshold for filtering markers based on the P-value.
        
        Returns:
        --------
            pd.DataFrame: A DataFrame containing markers with different genotype call rates
                          below the specified threshold. The DataFrame has two columns:
                          'SNP' and 'Failure', where 'Failure' is set to 'Different genotype call rate'.
        """

        if not isinstance(threshold, float):
            raise TypeError("threshold should be of type float.")
        if not isinstance(filename, Path):
            raise TypeError("filename should be of type Path.")

        # load .missing file
        df_diffmiss = pd.read_csv(
            filename,
            sep=r"\s+"
        )

        # filter markers with different genotype call rate
        fail_diffmiss = df_diffmiss[df_diffmiss['P']<threshold].reset_index(drop=True)
        fail_diffmiss = fail_diffmiss[['SNP']].copy()
        fail_diffmiss['Failure'] = 'Different genotype call rate'

        return fail_diffmiss
    
    def execute_variant_qc_pipeline(self, variant_params: dict) -> None:
        """Execute a comprehensive variant quality control pipeline.
        
        This method runs a series of quality control steps on genetic variants,
        including missing data analysis, genotype calling assessment, Hardy-Weinberg
        equilibrium testing, and variant filtering based on specified thresholds.
        
        Parameters
        ----------
        variant_params : dict
            Dictionary containing quality control parameters with the following keys:
            
            - 'chr-y' : bool or str
                Flag for chromosome Y analysis
            - 'miss_data_rate' : float
                Threshold for missing data rate filtering
            - 'diff_genotype_rate' : float
                Threshold for differential genotype call rate between cases/controls
            - 'hwe' : float
                Hardy-Weinberg equilibrium p-value threshold
            - 'maf' : float
                Minor allele frequency threshold for variant filtering
            - 'geno' : float
                Genotype call rate threshold for variant filtering
        
        Returns
        -------
        None
            This method performs quality control operations in-place and does not
            return any values.
        
        Notes
        -----
        The pipeline executes the following steps in order:
        
        1. Missing data rate computation (sex-stratified analysis)
        2. Case/control nonrandom missingness test
        3. Hardy-Weinberg equilibrium test
        4. Identification of variants failing QC thresholds
        5. Removal of variants that failed quality control
        
        Each step prints a colored status message indicating the current operation
        being performed.
        """
        
        variant_qc_steps = {
            'Missing data rate'         : (self.execute_missing_data_rate, {'chr_y': variant_params['chr_y']}),
            'Different genotype'        : (self.execute_different_genotype_call_rate, {}),
            'Get fail variants'         : (self.get_fail_variants, {'marker_call_rate_thres': variant_params['miss_data_rate'], 'case_controls_thres': variant_params['diff_genotype_rate']}),
            'Drop fail variants'        : (self.execute_drop_variants, {'maf': variant_params['maf'], 'geno': variant_params['geno'], 'hwe': variant_params['hwe']}),
        }

        variant_step_description = {
            'Missing data rate'         : 'Compute missing data rate for males and females',
            'Different genotype'        : 'Case/control nonrandom missingness test',
            'Hardy-Weinberg equilibrium': 'Hardy-Weinberg equilibrium test',
            'Get fail variants'         : 'Get variants that failed quality control',
            'Drop fail variants'        : 'Drop variants that failed quality control'
        }

        for name, (func, params) in variant_qc_steps.items():
            print(f"\033[34m{variant_step_description[name]}.\033[0m")
            func(**params)

        return
    
class VariantQCReport:
    """
    Handles visualization and reporting for variant quality control results.
    
    This class is responsible for generating plots and reports from variant QC analyses,
    following the Single Responsibility Principle by separating visualization concerns
    from the main VariantQC logic.
    
    Parameters:
    -----------
    output_path: Path 
        Directory path where plots and reports will be saved
        
    Raises:
    -------
    TypeError: 
        If output_path is not a Path object
    FileNotFoundError: 
        If output_path doesn't exist
        
    Attributes:
    -----------
    output_path: Path 
        Path to output directory where visualizations will be saved
    """

    def __init__(self, output_path: Path) -> None:
        """
        Initialize the VariantQCReport class.
        
        Parameters:
        -----------
        output_path: Path 
            Directory path where plots and reports will be saved
            
        Raises:
        -------
        TypeError: 
            If output_path is not a Path object
        FileNotFoundError: 
            If output_path doesn't exist
        """
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be of type Path")
        
        if not output_path.exists():
            raise FileNotFoundError("output_path is not a valid path")
            
        self.output_path = output_path

    def report_variant_qc(
        self,
        missing_data_rate_male: Path,
        missing_data_rate_female: Path,
        y_axis_cap: Optional[float]=100,
        missing_data_threshold: float=0.1
    ) -> None:
        """
        Generate comprehensive visualization reports for variant quality control results.
        
        This method creates histogram plots for missing data rates (separated by sex) and 
        Hardy-Weinberg equilibrium test results, providing visual assessment of QC metrics.
        
        Parameters
        ----------
        missing_data_rate_male : Path
            Path to the .vmiss file containing missing data statistics for male subjects
        missing_data_rate_female : Path
            Path to the .vmiss file containing missing data statistics for female subjects
        hwe_file : Path
            Path to the .hwe file containing Hardy-Weinberg equilibrium test results
        y_axis_cap : float, optional
            Maximum value for y-axis in histograms. If None, automatically determined. 
            Default is 100.
        missing_data_threshold : float, optional
            Threshold for missing data rate visualization (vertical line). Default is 0.2.
        hwe_threshold : float, optional
            P-value threshold for Hardy-Weinberg equilibrium test (vertical line). 
            Default is 5e-8.
            
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If any of the input files don't exist
        pd.errors.EmptyDataError
            If input files are empty or malformed
            
        Notes
        -----
        - Creates 'missing_data_male.svg' and 'missing_data_female.svg' plots
        - Creates 'hwe-histogram.svg' plot
        - All plots are saved in SVG format with 600 DPI resolution
        - Threshold lines are displayed in red with dashed style
        """
        
        # ==========================================================================================================
        #                                             MISSING DATA RATE
        # ==========================================================================================================

        # load .vmiss file for male subjects
        df_males = pd.read_csv(
            missing_data_rate_male,
            sep=r"\s+"
        )

        self._make_histogram(
            values  = df_males['F_MISS'], 
            output_name='missing_data_male',
            output_path=self.output_path,
            threshold=missing_data_threshold, 
            x_label='Ratio of missing data', 
            title='Missing data for males', 
            y_lim_cap=y_axis_cap
        )

        # load .vmiss file for female subjects
        df_females = pd.read_csv(
            missing_data_rate_female,
            sep=r"\s+"
        )

        self._make_histogram(
            values  = df_females['F_MISS'], 
            output_name='missing_data_female',
            output_path=self.output_path,
            threshold=missing_data_threshold, 
            x_label='Ratio of missing data', 
            title='Missing data for females', 
            y_lim_cap=y_axis_cap
        )

        return

    def _make_histogram(self, values: pd.Series, output_name: str, output_path: Path, threshold: float, x_label: str, title: str, y_lim_cap: Optional[float] = None, format: str = 'svg') -> None:
        """
        Create a histogram plot with a vertical threshold line and save to file.
        
        This private method generates standardized histogram visualizations for QC metrics
        with consistent styling, threshold indicators, and configurable output formats.
        
        Parameters
        ----------
        values : pd.Series
            Series containing the numeric values to plot in the histogram
        output_name : str
            Name of the output file (without extension)
        output_path : Path
            Directory path where the plot will be saved
        threshold : float
            Value where to draw the vertical threshold line
        x_label : str
            Label for the x-axis
        title : str
            Title of the plot
        y_lim_cap : float, optional
            Upper limit for y-axis. If None, automatically determined. Default is None.
        format : str, optional
            Output file format ('svg', 'pdf', 'png', etc.). Default is 'svg'.
        
        Returns
        -------
        None
            This method saves the plot to file but does not return any value.
        
        Notes
        -----
        - Uses 50 bins for histogram with color '#1B9E77'
        - Threshold line is displayed in red with dashed style
        - Plot is saved with 600 DPI resolution
        - Memory is properly managed with plt.close() after saving
        """

        plt.clf()

        fig_path = output_path / f"{output_name}"

        plt.hist(values, bins=50, color='#1B9E77')
        plt.xlabel(x_label)
        plt.ylabel('Number of SNPs')
        plt.ylim(0, y_lim_cap if y_lim_cap else None)
        plt.title(title)

        # Draw the vertical line indicating the cut off threshold
        plt.axvline(x=threshold, linestyle='--', color='red')

        plt.savefig(fig_path, format=format, dpi=600)
        plt.show(block=False)
        plt.close()

        return

class VariantQCCleanUp:

    def __init__(self, output_path: Path) -> None:
        self.output_path = output_path

    def clean_all(self) -> None:

        """Remove intermediate files from output directory.

        This method deletes temporary files created during sample QC steps:
        - Files ending with '.bed', '.bim', '.fam', '.vmiss', '.smiss', '.nosex', '.missing', '.hwe'

        Returns
        -------
        None

        Notes
        -----
        Only removes files if they exist. No error is raised if files are not found.
        """

        logger.info("Cleaning up intermediate files from output directory")

        # Remove intermediate files
        extensions = ['.bed', '.bim', '.fam', '.vmiss', '.smiss', '.nosex', '.missing', '.hwe']
        for ext in extensions:
            for file in self.output_path.glob(f'*{ext}'):
                if file.exists():
                    file.unlink()
                    logger.info(f"Deleted: {file}")

        return
