

import logging

import pandas as pd

from pathlib import Path

from ..core.executor import run_plink
from ..core.get_references import FetcherLDRegions, Fetcher1000Genome
from ..core.utils import get_optimal_threads, get_available_memory
from ..qc.ancestry_qc import ReferenceGenomicMerger

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class FstSummary:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, high_ld_file: Path=Path(), build: str = '38', recompute_merge: bool = True, reference_files: dict = dict()) -> None:
        """
        Initialize FstSummary object for Fst analysis.
        
        Parameters
        ----------
        input_path : Path
            Path to the directory containing input files
        input_name : str
            Name of the input file
        output_path : Path
            Path to the directory where results will be saved
        
        Raises
        ------
        TypeError
            If input types are incorrect for any parameter
        FileNotFoundError
            If input_path or output_path do not exist
        """

        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(input_name, str): 
            raise TypeError("input_name should be a string")
        if not input_path.exists():
            raise FileNotFoundError("input_path does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")
        if not isinstance(build, str):
            raise TypeError("build should be a string")
        if build not in ['37', '38']:
            raise ValueError("build should be either '37' or '38'")
        if not isinstance(reference_files, dict):
            raise TypeError("reference_files should be a dictionary")
        if not isinstance(high_ld_file, Path):
            raise TypeError("high_ld_file should be a Path object")
        if not high_ld_file.is_file():
            logger.info(f"High LD file not found at {high_ld_file}")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions()
            ld_fetcher.get_ld_regions()

            if ld_fetcher.ld_regions is None:
                raise FileNotFoundError("Could not fetch LD regions file.")
                
            high_ld_file = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_file}")

        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.recompute_merge = recompute_merge
        self.high_ld_regions = high_ld_file
        self.build = build

        # Convert reference_files values to Path objects if they are strings
        if reference_files:
            self.reference_files = {k: Path(v) if isinstance(v, str) else v for k, v in reference_files.items()}
        else:
            self.reference_files = reference_files

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

        self.results_dir = self.output_path / 'fst_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.merging_dir = self.results_dir / 'merging'
        self.merging_dir.mkdir(parents=True, exist_ok=True)

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
            output_name= 'cleaned-with-ref',
            high_ld_regions_file=self.high_ld_regions,
            reference_files = self.reference_files,
        )

        rgm.execute_rename_snpid()
        rgm.execute_filter_prob_snps()
        rgm.execute_ld_pruning(ind_pair=ind_pair)
        rgm.execute_fix_chromosome_mismatch()
        rgm.execute_fix_position_mismatch()
        rgm.execute_fix_allele_flip()
        rgm.execute_remove_mismatches()
        rgm.execute_merge_data()

        for file in self.merging_dir.iterdir():
            if file.is_file() and '-merged' not in file.name and file.suffix != '.log':
                file.unlink()

        return
    
    def add_population_tags(self) -> None:
        """
        Add population tags to the merged dataset.
        This method adds population super-population tags from the reference dataset to
        the merged dataset. It reads population information from the reference PSAM file,
        merges it with the study dataset, and assigns 'StPop' (study population) to samples
        not present in the reference dataset.

        Requirements
        ------------
            - Merged dataset files (.bed, .bim, .fam) must exist in the merging directory
            - Reference files dictionary must contain a valid 'psam' Path
        
        Raises
        ------
            FileNotFoundError: If any of the required merged files are not found
            ValueError: If the reference files dictionary doesn't contain a valid 'psam' Path
        
        Side Effects
        ------------
            - Creates a new tab-separated file with population tags at 
              {merging_dir}/cleaned-with-ref-merged-pop-tags.csv
            - Sets self.population_tags to the path of the created file
        
        Returns
        -------
            None
        """
       
        merged_bed = self.merging_dir / 'cleaned-with-ref-merged.bed'
        merged_bim = self.merging_dir / 'cleaned-with-ref-merged.bim'
        merged_fam = self.merging_dir / 'cleaned-with-ref-merged.fam'
        
        if not merged_bed.is_file():
            raise FileNotFoundError(f"Merged data file not found at {merged_bed}")
        if not merged_bim.is_file():
            raise FileNotFoundError(f"Merged BIM file not found at {merged_bim}")
        if not merged_fam.is_file():
            raise FileNotFoundError(f"Merged FAM file not found at {merged_fam}")

        if 'psam' not in self.reference_files or not isinstance(self.reference_files['psam'], Path):
            raise ValueError("Reference files dictionary must contain a valid 'psam' Path")

        reference_tags = self.reference_files['psam']
        
        df_tags = pd.read_csv(reference_tags, sep=r"\s+", engine='python')
        df_tags['ID'] = '0'
        df_tags = df_tags[['ID', '#IID', 'SuperPop']]
        df_tags = df_tags.rename(columns={'ID': 'ID1', '#IID': 'ID2', 'SuperPop': 'SuperPop'})

        logger.info(f"Population tags loaded from {reference_tags}")
        logger.info(f'Population tags columns: {df_tags.columns.tolist()}')

        df_merged_fam = pd.read_csv(merged_fam, sep=r"\s+", header=None, engine='python')
        df_merged_fam = df_merged_fam.rename(columns={0: 'ID1', 1: 'ID2'})
        df_merged_fam['ID1'] = df_merged_fam['ID1'].astype(str)
        df_merged_fam['ID2'] = df_merged_fam['ID2'].astype(str)

        logger.info(f"Merged BIM file loaded from {merged_fam}")
        logger.info(f'Merged BIM file columns: {df_merged_fam.columns.tolist()}')

        df = pd.merge(
            df_merged_fam[['ID1', 'ID2']],
            df_tags,
            on=['ID1', 'ID2'],
            how='left'
        )
        df['SuperPop'] = df['SuperPop'].fillna('StPop')
        logger.info(f"Added population tags to the merged dataset")

        self.population_tags = self.merging_dir / 'cleaned-with-ref-merged-pop-tags.csv'
        df.to_csv(
            self.population_tags,
            index=False,
            sep='\t'
        )

        return
    
    def compute_fst(self) -> None:
        """
        Compute FST (fixation index) statistics between populations.

        This method calculates FST statistics between each super-population in the dataset
        and a study population ('StPop'). The process involves:
        1. Reading population tags from the specified file
        2. For each unique super-population (except 'StPop'):
        - Creating population filter files (keep and within files)
        - Running PLINK commands to filter the dataset and compute FST statistics

        The method requires the following instance variables to be set:
            - population_tags: Path to a file containing population information
            - results_dir: Directory where results will be stored
            - merging_dir: Directory containing the merged genotype data

        Returns:
        --------
        None
        """

        df_tags = pd.read_csv(self.population_tags, sep=r"\s+", engine='python')

        files = dict()

        for pop in df_tags['SuperPop'].unique():
            if pop != 'StPop':
                df_temp = df_tags[(df_tags['SuperPop'] == pop) | (df_tags['SuperPop'] == 'StPop')].reset_index(drop=True)
                df_temp[['ID1', 'ID2']].to_csv(self.results_dir / f'keep-{pop}_StPop.txt', sep='\t', index=False, header=False)
                df_temp.to_csv(self.results_dir / f'within-{pop}_StPop.txt', index=False, header=False, sep='\t',)

                files[pop] = (self.results_dir / f'keep-{pop}_StPop.txt', self.results_dir / f'within-{pop}_StPop.txt')

                logger.info(f"Created keep and within files for population {pop}")

        input_file = self.merging_dir / 'cleaned-with-ref-merged'

        max_threads = get_optimal_threads()
        memory = get_available_memory()

        for key in files.keys():

            keep_file, within_file = files[key]
            output_file = self.results_dir / f'keep-{key}-StPop'

            run_plink([
                '--bfile', str(input_file),
                '--keep', str(keep_file),
                '--threads', str(max_threads),
                '--memory', str(int(memory)),
                '--make-bed',
                '--out', str(output_file)
            ])
            run_plink([
                '--bfile', str(output_file),
                '--fst',
                '--within', str(within_file),
                '--out', str(self.results_dir / f'fst-{key}-StPop')
            ])
        logger.info("Fst computation completed for all populations.")

        return
    
    def report_fst(self) -> pd.DataFrame:
        """
        Generate a report of Fst results.
        This method reads the Fst results from the results directory and generates a summary report.
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing the Fst results summary
        
        Raises
        ------
        FileNotFoundError
            If no Fst result files are found in the results directory.
        """

        df_summary = pd.DataFrame(columns=['SuperPop', 'Fst', 'WeightedFst'])

        # Get a list of all log files in the results directory
        log_files = [f for f in self.results_dir.iterdir() if f.is_file() and f.suffix == '.log']

        if not log_files:
            raise FileNotFoundError(f"No log files found in {self.results_dir}")

        logger.info(f"Found {len(log_files)} log files in {self.results_dir}")

        # Extract the population names from log file names
        # Assuming log files follow the pattern 'fst-{population}-StPop.log'
        files = {}
        for log_file in log_files:
            if log_file.stem.startswith('fst-') and log_file.stem.endswith('-StPop'):
                pop = log_file.stem.split('-')[1]
                files[pop] = log_file
                logger.info(f"Found Fst result file for population {pop} at {log_file}")

        if not files:
            raise FileNotFoundError(f"No Fst result files found in {self.results_dir}")

        for key in files.keys():
        
            log_file = files[key]
            fst = None
            weighted_fst = None
            with open(log_file, 'r') as f:

                lines = f.readlines()
                for line in lines:
                    if line.startswith('Mean Fst'):
                        fst = line.split(':')[1].strip()
                    if line.startswith('Weighted Fst'):
                        weighted_fst = line.split(':')[1].strip()

            if fst is None or weighted_fst is None:
                logger.warning(f"Could not find Fst values in {log_file}, skipping population {key}")
                continue

            df_summary = pd.concat([df_summary, pd.DataFrame({'SuperPop': [key], 'Fst': [fst], 'WeightedFst': [weighted_fst]})], ignore_index=True)

        df_summary.to_csv(
            self.results_dir / 'fst_summary.csv',
            index=False,
            sep='\t'
        )
        logger.info(f"Fst summary report generated at {self.results_dir / 'fst_summary.csv'}")

        for file in self.results_dir.iterdir():
            if file.is_file() and (file.suffix == '.bed' or file.suffix == '.bim' or file.suffix == '.fam'):
                file.unlink()

        return df_summary

    def execute_fst_pipeline(self, fst_params: dict) -> None:
        """
        Execute the complete Fst analysis pipeline.

        This method runs the full Fst workflow:
        1. Merging reference and study data
        2. Tagging samples with their (super-)population
        3. Computing Fst between the study population and each reference super-population
        4. Summarizing the Fst results into a report

        Parameters
        ----------
        fst_params : dict
            Dictionary containing pipeline parameters. Required keys:
            - ind_pair : list - LD pruning parameters [window, step, r2]

        Returns
        -------
        None
        """

        fst_steps = {
            'merge_study_reference': (self.merge_reference_study, {"ind_pair": fst_params['ind_pair']}),
            'add_population_tags':   (self.add_population_tags, {}),
            'compute_fst':            (self.compute_fst, {}),
            'report_fst':             (self.report_fst, {}),
        }

        step_description = {
            'merge_study_reference': "Merging reference genome with study genome",
            'add_population_tags':   "Tagging samples with their (super-)population",
            'compute_fst':            "Computing Fst between study and reference populations",
            'report_fst':             "Summarizing Fst results",
        }

        logger.info("=" * 70)
        logger.info("Starting Fst Pipeline")
        logger.info("=" * 70)

        result = None
        for step_num, (name, (func, params)) in enumerate(fst_steps.items(), 1):
            logger.info(f"\nSTEP {step_num}/{len(fst_steps)}: {step_description[name]}")
            logger.info("-" * 70)
            try:
                result = func(**params)
                logger.info(f"✓ Step {step_num} completed successfully")
            except Exception as e:
                logger.error(f"✗ Step {step_num} failed: {str(e)}")
                raise

        self.fst_summary = result

        logger.info("=" * 70)
        logger.info("Fst Pipeline completed successfully")
        logger.info("=" * 70)

        return