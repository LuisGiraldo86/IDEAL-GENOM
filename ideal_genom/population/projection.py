"""
Module to draw plots based on UMAP dimension reduction
"""

import umap
import warnings
import logging
from typing import Optional, Literal

import pandas as pd
# Set matplotlib to use non-interactive backend before importing pyplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE

from pathlib import Path

from ..core.executor import run_plink2
from ..core.utils import get_optimal_threads, get_available_memory
from ..core.get_references import FetcherLDRegions

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class PCAReduction:

    def __init__(self, input_path: Path, input_name: str, output_path: Path, build: str = '38', high_ld_regions: Optional[Path] = None, generate_plot: bool = True) -> None:

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
        if not isinstance(generate_plot, bool):
            raise TypeError("generate_plot should be a boolean")
        
        if high_ld_regions is not None:
            if not isinstance(high_ld_regions, Path):
                raise TypeError("high_ld_regions should be a Path object")
#            if not high_ld_regions.is_file():
#                raise FileNotFoundError("high_ld_regions file does not exist")
        
        if high_ld_regions is None or (not high_ld_regions.exists()):
            logger.info(f"High LD file not provided.")
            logger.info('High LD file will be fetched from the package')
            
            ld_fetcher = FetcherLDRegions(build=build)
            ld_fetcher.get_ld_regions()

            if ld_fetcher.ld_regions is None:
                raise FileNotFoundError("Could not fetch LD regions file.")
                
            high_ld_regions = ld_fetcher.ld_regions
            logger.info(f"High LD file fetched from the package and saved at {high_ld_regions}")
        
        self.input_path = input_path
        self.input_name = input_name
        self.output_path= output_path
        self.build = build
        self.high_ld_regions = high_ld_regions
        self.generate_plot = generate_plot

        self.pruned_data: Optional[Path] = None
        
        return
    
    def execute_ld_pruning(self, maf: float = 0.001, geno: float = 0.1, mind: float = 0.2, hwe: float = 5e-8, ind_pair:list = [50, 5, 0.2]) -> None:
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
        - Excludes high LD regions specified in `self.high_ld_regions`.
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
            '--bfile', str(self.input_path / self.input_name),
            '--maf', str(maf),
            '--geno', str(geno),
            '--mind', str(mind),
            '--hwe', str(hwe),
            '--exclude', 'range', str(self.high_ld_regions),
            '--indep-pairwise', str(ind_pair[0]), str(ind_pair[1]), str(ind_pair[2]),
            '--threads', str(max_threads),
            '--out', str(self.output_path / self.input_name)
        ])

        # Execute PLINK2 command: prune study data and creates a filtered binary file
        run_plink2([
            '--bfile', str(self.input_path / self.input_name),
            '--extract', str((self.output_path / self.input_name).with_suffix('.prune.in')),
            '--threads', str(max_threads),
            '--maf', str(maf),
            '--geno', str(geno),
            '--mind', str(mind),
            '--hwe', str(hwe),
            '--make-bed',
            '--out', str(self.output_path / (self.input_name + '-pruned'))
        ])

        self.pruned_data = self.output_path / (self.input_name+'-pruned')

        return

    def execute_pca(self, pca: int = 20, maf: float = 0.01) -> None:
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
            '--bfile', str(self.pruned_data),
            '--maf', str(maf),
            '--out', str(self.output_path / (self.input_name + '-pca')),
            '--pca', str(pca),
            '--memory', str(int(memory)),
            '--threads', str(max_threads)
        ])

        self.eigenvectors = self.output_path / (self.input_name+'-pca.eigenvec')
        self.eigenvalues = self.output_path / (self.input_name+'-pca.eigenval')

        return
    
    def execute_pcareduction_pipeline(self, maf: float = 0.001, geno: float = 0.1, mind: float = 0.2, hwe: float = 5e-8, ind_pair:list = [50, 5, 0.2], pca: int = 20, case_control_markers: bool = True) -> None:
        """
        Execute the full preparation pipeline: LD pruning followed by PCA.

        This method sequentially performs LD pruning on the genetic data and then
        computes principal components using PCA. It combines the functionalities of
        `execute_ld_pruning` and `execute_pca` methods.

        Parameters
        ----------
        ind_pair : list
            A list containing three elements for LD pruning:
            - ind_pair[0] (int): Window size in SNPs
            - ind_pair[1] (int): Number of SNPs to shift the window at each step
            - ind_pair[2] (float): r² threshold for pruning
        pca : int, default=20
            Number of principal components to calculate.
            Must be a positive integer.
        maf : float, default=0.01
            Minor allele frequency threshold for filtering variants.
            Must be between 0 and 0.5
        Returns
        -------
        None
        """

        self.execute_ld_pruning(maf=maf, geno=geno, mind=mind, hwe=hwe, ind_pair=ind_pair)
        self.execute_pca(pca=pca, maf=maf)

        plot_pca_2d = Plot2D(
            output_dir=self.output_path
        )
        if self.generate_plot:

            df_eigenvec = pd.read_csv(
                self.eigenvectors,
                sep=r'\s+',
                engine='python')
            
            # Map eigenvector columns to Plot2D expected format
            # Rename #FID to ID1 and IID to ID2 for compatibility with Plot2D
            df_eigenvec = df_eigenvec.rename(columns={'#FID': 'ID1', 'IID': 'ID2'})
            
            plot_pca_2d.prepare_metadata(
                case_control_markers=case_control_markers,
                fam_file= self.input_path / (self.input_name + '.fam')
            )

            plot_pca_2d.generate_plot(
                data=df_eigenvec,
                x_col='PC1',
                y_col='PC2',
                plot_name='pca_2d_plot',
                xlabel='Principal Component 1',
                ylabel='Principal Component 2',
                title='PCA 2D Projection',
            )

        logger.info("Preparation pipeline completed successfully.")

        return

class UMAPReduction:
    """
    Class for performing UMAP dimensionality reduction on PCA eigenvectors.
    
    This class handles UMAP transformation of high-dimensional PCA data into
    2D space for visualization. Use Plot2D class for generating plots.
    """

    def __init__(self, eigenvector: Path, output_path: Path) -> None:
        """
        Initialize UMAPReduction object.
        
        Parameters
        ----------
        eigenvector : Path
            Path to the eigenvector file (.eigenvec) from PCA analysis
        output_path : Path
            Path to the directory where results will be saved

        Raises
        ------
        TypeError
            If input types are incorrect
        FileNotFoundError
            If eigenvector file or output_path do not exist
        
        Notes
        -----
        Creates 'umap_results' directory in the output path.
        """
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(eigenvector, Path):
            raise TypeError("eigenvector should be a Path object")
        if not eigenvector.exists():
            raise FileNotFoundError("eigenvector file does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")

        self.eigenvector = eigenvector
        self.output_path = output_path

        self.results_dir = self.output_path / 'umap_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache for eigenvector data to avoid repeated file loading
        self._eigenvec_cache: Optional[pd.DataFrame] = None
    
    def _load_eigenvectors(self) -> pd.DataFrame:
        """
        Load eigenvector data from file with caching.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns ['ID1', 'ID2', 'pca_1', 'pca_2', ...]
        
        Raises
        ------
        FileNotFoundError
            If eigenvector file does not exist
        """
        if self._eigenvec_cache is not None:
            return self._eigenvec_cache.copy()
        
        eigenvec_file = self.eigenvector
        
        if not eigenvec_file.exists():
            raise FileNotFoundError(f"Eigenvector file not found: {eigenvec_file}")
        
        df_eigenvec = pd.read_csv(
            eigenvec_file, 
            sep=r'\s+',
            engine='python'
        )
        logger.info(f"Eigenvector file loaded from {eigenvec_file}")
        logger.info(f"Eigenvector file has {df_eigenvec.shape[0]} rows and {df_eigenvec.shape[1]} columns")
        
        # Validate eigenvector file format
        if df_eigenvec.shape[1] < 3:
            raise ValueError(
                f"Invalid eigenvector file format. Expected at least 3 columns "
                f"(ID1, ID2, PC1), but found {df_eigenvec.shape[1]} columns. "
                f"This suggests the PCA step may have failed or produced invalid output. "
                f"Please check the PCA preparation step and ensure it completed successfully."
            )
        
        # Rename columns
        num_pc = df_eigenvec.shape[1] - 2
        new_cols = [f"pca_{k}" for k in range(1, num_pc + 1)]
        df_eigenvec.columns = ['ID1', 'ID2'] + new_cols
        
        self._eigenvec_cache = df_eigenvec
        return df_eigenvec.copy()
    
    def fit_transform(self,
                     n_neighbors: int = 15,
                     min_dist: float = 0.1,
                     metric: str = 'euclidean',
                     random_state: Optional[int] = None,
                     n_components: int = 2,
                     umap_kwargs: Optional[dict] = None) -> pd.DataFrame:
        """
        Perform UMAP dimensionality reduction on PCA eigenvectors.
        
        Parameters
        ----------
        n_neighbors : int, default=15
            Number of neighbors for UMAP manifold approximation. Must be positive.
        min_dist : float, default=0.1
            Minimum distance between points in low-dimensional space. Must be non-negative.
        metric : str, default='euclidean'
            Distance metric for UMAP (e.g., 'euclidean', 'cosine', 'manhattan')
        random_state : int, optional
            Random seed for reproducibility. Must be non-negative.
        n_components : int, default=2
            Number of dimensions in the output
        umap_kwargs : dict, optional
            Additional keyword arguments to pass to UMAP constructor.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns ['ID1', 'ID2', 'umap_1', 'umap_2', ...]
        
        Raises
        ------
        TypeError
            If parameters are not of correct type
        ValueError
            If parameter values are invalid
        """
        # Validate parameters
        if not isinstance(n_neighbors, int):
            raise TypeError("n_neighbors should be an integer")
        if n_neighbors <= 0:
            raise ValueError("n_neighbors should be a positive integer")
        if not isinstance(min_dist, float):
            raise TypeError("min_dist should be a float")
        if min_dist < 0:
            raise ValueError("min_dist should be a non-negative float")
        if not isinstance(metric, str):
            raise TypeError("metric should be a string")
        if random_state is not None:
            if not isinstance(random_state, int):
                raise TypeError("random_state should be an integer")
            if random_state < 0:
                raise ValueError("random_state should be a non-negative integer")
        if not isinstance(n_components, int):
            raise TypeError("n_components should be an integer")
        if n_components <= 0:
            raise ValueError("n_components should be a positive integer")
        if umap_kwargs is None:
            umap_kwargs = {}
        if not isinstance(umap_kwargs, dict):
            raise TypeError("umap_kwargs should be a dictionary")
        
        logger.info(f"Performing UMAP reduction")
        logger.info(f"Parameters: n_neighbors={n_neighbors}, min_dist={min_dist}, metric={metric}")
        
        # Load eigenvector data
        df_eigenvec = self._load_eigenvectors()
        
        df_ids = df_eigenvec[['ID1', 'ID2']].copy()
        pca_cols = [col for col in df_eigenvec.columns if col.startswith('pca_')]
        df_vals = df_eigenvec[pca_cols].to_numpy()
        
        # Instantiate UMAP
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state,
            **umap_kwargs
        )
        
        # Perform UMAP transformation
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            umap_proj = reducer.fit_transform(df_vals)
            
            if w:
                for warning in w:
                    logger.warning(f"UMAP warning: {warning.message}")
        
        # Create result DataFrame
        umap_cols = [f"umap_{k}" for k in range(1, n_components + 1)]
        df_umap = pd.DataFrame(data=umap_proj, columns=umap_cols)
        df_result = pd.concat([df_ids.reset_index(drop=True), df_umap], axis=1)
        
        logger.info(f"UMAP reduction completed: {df_result.shape[0]} samples, {n_components} dimensions")
        
        return df_result


class TSNEReduction:
    """
    Class for performing t-SNE dimensionality reduction on PCA eigenvectors.
    
    This class handles t-SNE transformation of high-dimensional PCA data into
    2D or 3D space for visualization. Use Plot2D class for generating plots.
    """

    def __init__(self, eigenvector: Path, output_path: Path) -> None:
        """
        Initialize TSNEReduction object.
        
        Parameters
        ----------
        eigenvector : Path
            Path to the eigenvector file (.eigenvec) from PCA analysis
        output_path : Path
            Path to the directory where results will be saved

        Raises
        ------
        TypeError
            If input types are incorrect
        FileNotFoundError
            If eigenvector file or output_path do not exist
        
        Notes
        -----
        Creates 'tsne_results' directory in the output path.
        """
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(eigenvector, Path):
            raise TypeError("eigenvector should be a Path object")
        if not eigenvector.exists():
            raise FileNotFoundError("eigenvector file does not exist")
        if not output_path.exists():
            raise FileNotFoundError("output_path does not exist")

        self.eigenvector = eigenvector
        self.output_path = output_path

        self.results_dir = self.output_path / 'tsne_results' 
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache for eigenvector data to avoid repeated file loading
        self._eigenvec_cache: Optional[pd.DataFrame] = None
    
    def _load_eigenvectors(self) -> pd.DataFrame:
        """
        Load eigenvector data from file with caching.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns ['ID1', 'ID2', 'pca_1', 'pca_2', ...]
        
        Raises
        ------
        FileNotFoundError
            If eigenvector file does not exist
        """
        if self._eigenvec_cache is not None:
            return self._eigenvec_cache.copy()
        
        eigenvec_file = self.eigenvector
        
        if not eigenvec_file.exists():
            raise FileNotFoundError(f"Eigenvector file not found: {eigenvec_file}")
        
        df_eigenvec = pd.read_csv(
            eigenvec_file, 
            sep=r'\s+',
            engine='python'
        )
        logger.info(f"Eigenvector file loaded from {eigenvec_file}")
        logger.info(f"Eigenvector file has {df_eigenvec.shape[0]} rows and {df_eigenvec.shape[1]} columns")
        
        # Validate eigenvector file format
        if df_eigenvec.shape[1] < 3:
            raise ValueError(
                f"Invalid eigenvector file format. Expected at least 3 columns "
                f"(ID1, ID2, PC1), but found {df_eigenvec.shape[1]} columns. "
                f"This suggests the PCA step may have failed or produced invalid output. "
                f"Please check the PCA preparation step and ensure it completed successfully."
            )
        
        # Rename columns
        num_pc = df_eigenvec.shape[1] - 2
        new_cols = [f"pca_{k}" for k in range(1, num_pc + 1)]
        df_eigenvec.columns = ['ID1', 'ID2'] + new_cols
        
        self._eigenvec_cache = df_eigenvec
        return df_eigenvec.copy()
    
    def fit_transform(self,
                     n_components: int = 2,
                     perplexity: float = 30.0,
                     learning_rate: float = 200.0,
                     n_iter: int = 1000,
                     metric: str = 'euclidean',
                     random_state: Optional[int] = None,
                     early_exaggeration: float = 12.0,
                     init: Literal['pca', 'random'] = 'pca',
                     tsne_kwargs: Optional[dict] = None) -> pd.DataFrame:
        """
        Perform t-SNE dimensionality reduction on PCA eigenvectors.
        
        Parameters
        ----------
        n_components : int, default=2
            Number of dimensions in the output (typically 2 or 3)
        perplexity : float, default=30.0
            Related to number of nearest neighbors. Should be between 5 and 50.
            Larger datasets require larger perplexity.
        learning_rate : float, default=200.0
            Learning rate for t-SNE optimization. Usually between 10.0 and 1000.0.
        n_iter : int, default=1000
            Maximum number of iterations for optimization
        metric : str, default='euclidean'
            Distance metric to use ('euclidean', 'manhattan', 'cosine', etc.)
        random_state : int, optional
            Random seed for reproducibility. Must be non-negative.
        early_exaggeration : float, default=12.0
            Controls how tight natural clusters are in the original space
        init : str, default='pca'
            Initialization method ('pca' or 'random')
        tsne_kwargs : dict, optional
            Additional keyword arguments to pass to TSNE constructor.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns ['ID1', 'ID2', 'tsne_1', 'tsne_2', ...]
        
        Raises
        ------
        TypeError
            If parameters are not of correct type
        ValueError
            If parameter values are invalid
        
        Notes
        -----
        t-SNE is computationally expensive. For large datasets (>10,000 samples),
        consider using perplexity between 30-50 and reducing n_iter if needed.
        """
        # Validate parameters
        if not isinstance(n_components, int):
            raise TypeError("n_components should be an integer")
        if n_components <= 0:
            raise ValueError("n_components should be a positive integer")
        if not isinstance(perplexity, (int, float)):
            raise TypeError("perplexity should be a number")
        if perplexity <= 0:
            raise ValueError("perplexity should be positive")
        if not isinstance(learning_rate, (int, float)):
            raise TypeError("learning_rate should be a number")
        if learning_rate <= 0:
            raise ValueError("learning_rate should be positive")
        if not isinstance(n_iter, int):
            raise TypeError("n_iter should be an integer")
        if n_iter <= 0:
            raise ValueError("n_iter should be a positive integer")
        if not isinstance(metric, str):
            raise TypeError("metric should be a string")
        if random_state is not None:
            if not isinstance(random_state, int):
                raise TypeError("random_state should be an integer")
            if random_state < 0:
                raise ValueError("random_state should be a non-negative integer")
        if not isinstance(early_exaggeration, (int, float)):
            raise TypeError("early_exaggeration should be a number")
        if early_exaggeration <= 0:
            raise ValueError("early_exaggeration should be positive")
        if init not in ['pca', 'random']:
            raise ValueError("init should be either 'pca' or 'random'")
        if tsne_kwargs is None:
            tsne_kwargs = {}
        if not isinstance(tsne_kwargs, dict):
            raise TypeError("tsne_kwargs should be a dictionary")
        
        logger.info(f"Performing t-SNE reduction")
        logger.info(f"Parameters: n_components={n_components}, perplexity={perplexity}, learning_rate={learning_rate}")
        
        # Load eigenvector data
        df_eigenvec = self._load_eigenvectors()
        
        df_ids = df_eigenvec[['ID1', 'ID2']].copy()
        pca_cols = [col for col in df_eigenvec.columns if col.startswith('pca_')]
        df_vals = df_eigenvec[pca_cols].to_numpy()
        
        # Instantiate t-SNE
        reducer = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=n_iter,
            metric=metric,
            random_state=random_state,
            early_exaggeration=early_exaggeration,
            init=init,
            verbose=1,  # Show progress
            **tsne_kwargs
        )
        
        # Perform t-SNE transformation
        logger.info(f"Starting t-SNE optimization (this may take several minutes)...")
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            tsne_proj = reducer.fit_transform(df_vals)
            
            if w:
                for warning in w:
                    logger.warning(f"t-SNE warning: {warning.message}")
        
        # Create result DataFrame
        tsne_cols = [f"tsne_{k}" for k in range(1, n_components + 1)]
        df_tsne = pd.DataFrame(data=tsne_proj, columns=tsne_cols)
        df_result = pd.concat([df_ids.reset_index(drop=True), df_tsne], axis=1)
        
        logger.info(f"t-SNE reduction completed: {df_result.shape[0]} samples, {n_components} dimensions")
        
        return df_result


class Plot2D:
    """
    Class for generating 2D scatter plots with metadata integration.
    
    This class handles the preparation of metadata (color hue files, case-control markers)
    and generates publication-quality 2D scatter plots for dimensionality reduction results.
    """
    
    def __init__(self, output_dir: Path) -> None:
        """
        Initialize Plot2D object.
        
        Parameters
        ----------
        output_dir : Path
            Directory where plots will be saved
        
        Raises
        ------
        TypeError
            If output_dir is not a Path object
        FileNotFoundError
            If output_dir does not exist
        """
        if not isinstance(output_dir, Path):
            raise TypeError("output_dir should be a Path object")
        if not output_dir.exists():
            raise FileNotFoundError(f"output_dir does not exist: {output_dir}")
        
        self.output_dir = output_dir
        self.metadata: Optional[pd.DataFrame] = None
    
    def prepare_metadata(self, 
                        color_hue_file: Optional[Path] = None,
                        case_control_markers: bool = False,
                        fam_file: Optional[Path] = None) -> Optional[pd.DataFrame]:
        """
        Prepare metadata DataFrame from color hue file and/or case-control markers.
        
        Parameters
        ----------
        color_hue_file : Path, optional
            Path to tab-separated file with metadata for coloring. Must have at least
            3 columns: ID1, ID2, and a metadata column.
        case_control_markers : bool, default=False
            Whether to load case-control labels from .fam file
        fam_file : Path, optional
            Path to .fam file containing case-control information. Required if
            case_control_markers is True.
        
        Returns
        -------
        pd.DataFrame
            Metadata DataFrame with columns ['ID1', 'ID2', ...] or None if no metadata
        
        Raises
        ------
        FileNotFoundError
            If specified files don't exist
        TypeError
            If parameters are of incorrect type
        """
        df_color_hue = None
        df_fam = None
        
        # Validate inputs
        if color_hue_file is not None:
            if not isinstance(color_hue_file, Path):
                raise TypeError("color_hue_file should be a Path object")
            if not color_hue_file.is_file():
                raise FileNotFoundError(f"color_hue_file not found at {color_hue_file}")
        
        if not isinstance(case_control_markers, bool):
            raise TypeError("case_control_markers should be a boolean")
        
        if case_control_markers and fam_file is None:
            raise ValueError("fam_file must be provided when case_control_markers is True")
        
        # Load color hue file if provided
        if color_hue_file is not None:
            df_color_hue = pd.read_csv(color_hue_file, sep='\t')
            logger.info(f"Color hue file loaded from {color_hue_file}")
            logger.info(f"Column '{df_color_hue.columns[2]}' will be used for color hue")
            df_color_hue.columns = ["ID1", "ID2", df_color_hue.columns[2]]
            logger.info(f"Color hue file has {df_color_hue.shape[0]} rows and {df_color_hue.shape[1]} columns")
        
        # Load case-control markers if requested
        if case_control_markers and fam_file is not None:
            if not isinstance(fam_file, Path):
                raise TypeError("fam_file should be a Path object")
            if not fam_file.exists():
                raise FileNotFoundError(f".fam file not found at {fam_file}")
            
            df_fam = pd.read_csv(fam_file, sep=r'\s+', engine='python')
            logger.info(f"Case-control labels loaded from {fam_file}")
            
            df_fam.columns = ["ID1", "ID2", "F_ID", "M_ID", "Sex", "Phenotype"]
            recode = {1: 'Control', 2: 'Patient'}
            df_fam["Phenotype"] = df_fam["Phenotype"].map(recode)
            df_fam = df_fam[['ID1', 'ID2', 'Phenotype']].copy()
            logger.info(f"Case-control markers has {df_fam.shape[0]} rows and {df_fam.shape[1]} columns")
        
        # Merge or return appropriate DataFrame
        if df_color_hue is not None and df_fam is not None:
            self.metadata = df_color_hue.merge(df_fam, on=['ID1', 'ID2'], how='inner')
            logger.info(f"Color hue file merged with case-control markers")
            logger.info(f"Merged metadata has {self.metadata.shape[0]} rows and {self.metadata.shape[1]} columns")
        elif df_color_hue is not None:
            self.metadata = df_color_hue
            logger.info(f"Using color hue file as metadata")
        elif df_fam is not None:
            self.metadata = df_fam
            logger.info(f"Using case-control markers as metadata")
        else:
            self.metadata = None
            logger.info(f"No metadata provided")
        
        return self.metadata
    
    def generate_plot(self,
                     data: pd.DataFrame,
                     x_col: str,
                     y_col: str,
                     plot_name: str,
                     hue_col: Optional[str] = None,
                     style_col: Optional[str] = None,
                     title: Optional[str] = None,
                     xlabel: Optional[str] = None,
                     ylabel: Optional[str] = None,
                     figsize: tuple = (5, 5),
                     dpi: int = 500,
                     format: str = 'pdf',
                     marker: str = '.',
                     marker_size: int = 10,
                     alpha: float = 0.5,
                     equal_aspect: bool = True,
                     legend_params: Optional[dict] = None) -> Path:
        """
        Generate a 2D scatter plot with optional metadata coloring and styling.
        
        Parameters
        ----------
        data : pd.DataFrame
            DataFrame containing the 2D coordinates and IDs (must have 'ID1', 'ID2' columns)
        x_col : str
            Column name for x-axis values
        y_col : str
            Column name for y-axis values
        plot_name : str
            Name of the output plot file
        hue_col : str, optional
            Column name for point coloring. If None and metadata exists, uses third column
            or 'Phenotype' if available.
        style_col : str, optional
            Column name for point styling (different markers)
        title : str, optional
            Plot title
        xlabel : str, optional
            X-axis label. If None, uses x_col.
        ylabel : str, optional
            Y-axis label. If None, uses y_col.
        figsize : tuple, default=(5, 5)
            Figure size in inches (width, height)
        dpi : int, default=500
            Resolution for saving the plot
        format : str, default='pdf'
            Output format ('pdf', 'png', 'jpeg', 'svg')
        marker : str, default='.'
            Marker style for scatter plot
        marker_size : int, default=10
            Size of markers
        alpha : float, default=0.5
            Transparency of markers (0-1)
        equal_aspect : bool, default=True
            Whether to set equal aspect ratio
        legend_params : dict, optional
            Parameters for legend customization (bbox_to_anchor, ncols, fontsize, etc.)
        
        Returns
        -------
        Path
            Path to the saved plot file
        
        Raises
        ------
        ValueError
            If required columns are missing or hue_col not found
        TypeError
            If parameters are of incorrect type
        """
        # Validate inputs
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data should be a pandas DataFrame")
        if 'ID1' not in data.columns or 'ID2' not in data.columns:
            raise ValueError("data must contain 'ID1' and 'ID2' columns")
        if x_col not in data.columns:
            raise ValueError(f"Column '{x_col}' not found in data")
        if y_col not in data.columns:
            raise ValueError(f"Column '{y_col}' not found in data")
        if format not in ['pdf', 'png', 'jpeg', 'svg']:
            raise ValueError("format should be one of 'pdf', 'png', 'jpeg', or 'svg'")
        
        # Merge with metadata if available
        if self.metadata is not None:
            df_plot = data.merge(self.metadata, on=['ID1', 'ID2'], how='inner')
            logger.info(f"Data merged with metadata: {df_plot.shape[0]} samples")
            
            # Auto-determine hue_col if not specified
            if hue_col is None:
                if 'Phenotype' in df_plot.columns:
                    hue_col = 'Phenotype'
                elif self.metadata.shape[1] > 2:
                    hue_col = self.metadata.columns[2]
        else:
            df_plot = data.copy()
            logger.info(f"No metadata to merge. Plotting {df_plot.shape[0]} samples")
        
        # Validate hue_col if specified
        if hue_col is not None and hue_col not in df_plot.columns:
            raise ValueError(f"hue_col '{hue_col}' not found in data")
        
        # Validate style_col if specified
        if style_col is not None and style_col not in df_plot.columns:
            raise ValueError(f"style_col '{style_col}' not found in data")
        
        logger.info(f"Generating 2D plot: {plot_name}")
        logger.info(f"Plot dimensions: {x_col} vs {y_col}")
        
        # Create plot
        sns.set_context(font_scale=0.9)
        fig, ax = plt.subplots(figsize=figsize)
        
        # Prepare scatter plot kwargs
        scatter_kwargs = {
            'data': df_plot,
            'x': x_col,
            'y': y_col,
            'hue': hue_col,
            'marker': marker,
            's': marker_size,
            'alpha': alpha,
            'ax': ax
        }
        
        if style_col is not None:
            scatter_kwargs['style'] = style_col
        
        # Generate scatter plot
        sns.scatterplot(**scatter_kwargs)
        
        # Configure legend
        if self.metadata is not None and ax.get_legend() is not None:
            if legend_params is None:
                legend_params = {
                    'bbox_to_anchor': (0., 1.02, 1., .102),
                    'loc': 'lower left',
                    'ncols': 3,
                    'mode': 'expand',
                    'borderaxespad': 0.,
                    'fontsize': 7,
                    'markerscale': 2
                }
            plt.legend(**legend_params)
        
        # Set labels and title
        ax.set_xlabel(xlabel if xlabel else x_col, fontsize=7)
        ax.set_ylabel(ylabel if ylabel else y_col, fontsize=7)
        if title:
            ax.set_title(title, fontsize=9)
        
        # Set tick label size
        ax.tick_params(axis='both', labelsize=7)
        
        # Set aspect ratio
        if equal_aspect:
            ax.set_aspect('equal', adjustable='datalim')
        
        plt.tight_layout()
        
        # Save plot
        output_path = self.output_dir / plot_name
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
        
        logger.info(f"Plot saved to {output_path}")
        
        return output_path


class Plot3D:
    pass


class DimensionalityReductionPipeline:
    """
    Pipeline for running PCA preparation and dimensionality reduction workflows.
    
    This class orchestrates the complete workflow from raw genetic data to
    dimensionality reduction visualizations, including:
    1. PCA preparation (LD pruning + PCA)
    2. Optional UMAP reduction
    3. Optional t-SNE reduction
    4. Automated plotting with metadata
    """
    
    def __init__(self, 
                 input_path: Path,
                 input_name: str,
                 output_path: Path,
                 build: str = '38',
                 high_ld_regions_file: Optional[Path] = None,
                 generate_plot: bool = True) -> None:
        """
        Initialize the dimensionality reduction pipeline.
        
        Parameters
        ----------
        input_path : Path
            Path to directory containing input genetic data files (.bed/.bim/.fam)
        input_name : str
            Base name of input files (without extension)
        output_path : Path
            Path to directory where all results will be saved
        build : str, default='38'
            Genome build version ('37' or '38')
        high_ld_regions_file : Path, optional
            Path to file containing high LD regions. If None, will be fetched automatically.
        generate_plot : bool, default=True
            Whether to generate plots automatically
        
        Raises
        ------
        TypeError
            If input types are incorrect
        FileNotFoundError
            If input_path or output_path don't exist
        ValueError
            If build is not '37' or '38'
        """
        if not isinstance(input_path, Path):
            raise TypeError("input_path should be a Path object")
        if not isinstance(output_path, Path):
            raise TypeError("output_path should be a Path object")
        if not isinstance(input_name, str):
            raise TypeError("input_name should be a string")
        if not input_path.exists():
            raise FileNotFoundError(f"input_path does not exist: {input_path}")
        if not output_path.exists():
            raise FileNotFoundError(f"output_path does not exist: {output_path}")
        if build not in ['37', '38']:
            raise ValueError("build should be either '37' or '38'")
        
        self.input_path = input_path
        self.input_name = input_name
        self.output_path = output_path
        self.build = build
        self.high_ld_regions = high_ld_regions_file
        self.generate_plot = generate_plot
        
        # Create main results directory
        self.results_dir = self.output_path / 'dimensionality_reduction_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.pca_preparator: Optional[PCAReduction] = None
        self.umap_reducer: Optional[UMAPReduction] = None
        self.tsne_reducer: Optional[TSNEReduction] = None
        self.plotter: Optional[Plot2D] = None
        
        # Store results
        self.eigenvector_file: Optional[Path] = None
        self.eigenvalue_file: Optional[Path] = None
        self.umap_results: Optional[pd.DataFrame] = None
        self.tsne_results: Optional[pd.DataFrame] = None
        
        logger.info(f"Pipeline initialized for {input_name}")
        logger.info(f"Results will be saved to {self.results_dir}")
    
    def _is_parameter_grid(self, params: Optional[dict]) -> bool:
        """
        Check if any parameter values are lists, indicating a parameter grid search.
        
        Parameters
        ----------
        params : dict, optional
            Parameter dictionary to check
            
        Returns
        -------
        bool
            True if any parameter values are lists (grid search), False otherwise
        """
        if params is None:
            return False
        return any(isinstance(v, list) for v in params.values())
    
    def execute_pca_preparation(self,
                           maf: float = 0.001,
                           geno: float = 0.1,
                           mind: float = 0.2,
                           hwe: float = 5e-8,
                           ind_pair: list = [50, 5, 0.2],
                           pca: int = 20,
                           case_control_markers: bool = False) -> Path:
        """
        Run PCA preparation: LD pruning and principal component analysis.
        
        Parameters
        ----------
        maf : float, default=0.001
            Minor allele frequency threshold
        geno : float, default=0.1
            Genotype missingness threshold
        mind : float, default=0.2
            Sample missingness threshold
        hwe : float, default=5e-8
            Hardy-Weinberg equilibrium p-value threshold
        ind_pair : list, default=[50, 5, 0.2]
            LD pruning parameters: [window_size, step_size, r2_threshold]
        pca : int, default=20
            Number of principal components to calculate
        
        Returns
        -------
        Path
            Path to the generated eigenvector file
        
        Notes
        -----
        This step is required before running UMAP or t-SNE reductions.
        """
        logger.info("=" * 80)
        logger.info("STEP 1: PCA Preparation")
        logger.info("=" * 80)
        
        # Initialize PCA preparator
        self.pca_preparator = PCAReduction(
            input_path=self.input_path,
            input_name=self.input_name,
            output_path=self.results_dir,
            build=self.build,
            high_ld_regions=self.high_ld_regions
        )
        
        # Run the preparation pipeline
        self.pca_preparator.execute_pcareduction_pipeline(
            maf=maf,
            geno=geno,
            mind=mind,
            hwe=hwe,
            ind_pair=ind_pair,
            pca=pca,
            case_control_markers=case_control_markers
        )
        
        # Store eigenvector and eigenvalue file paths
        self.eigenvector_file = self.pca_preparator.eigenvectors
        self.eigenvalue_file = self.pca_preparator.eigenvalues
        
        logger.info(f"PCA preparation completed")
        logger.info(f"Eigenvector file: {self.eigenvector_file}")
        logger.info(f"Eigenvalue file: {self.eigenvalue_file}")
        
        return self.eigenvector_file
    
    def execute_umap(self,
                 n_neighbors: int = 15,
                 min_dist: float = 0.1,
                 metric: str = 'euclidean',
                 random_state: Optional[int] = None,
                 n_components: int = 2,
                 umap_kwargs: Optional[dict] = None) -> pd.DataFrame:
        """
        Run UMAP dimensionality reduction.
        
        Parameters
        ----------
        n_neighbors : int, default=15
            Number of neighbors for UMAP
        min_dist : float, default=0.1
            Minimum distance between points
        metric : str, default='euclidean'
            Distance metric
        random_state : int, optional
            Random seed for reproducibility
        n_components : int, default=2
            Number of output dimensions
        umap_kwargs : dict, optional
            Additional UMAP parameters
        
        Returns
        -------
        pd.DataFrame
            UMAP results with columns ['ID1', 'ID2', 'umap_1', 'umap_2', ...]
        
        Raises
        ------
        RuntimeError
            If PCA preparation hasn't been run yet
        """
        if self.eigenvector_file is None:
            raise RuntimeError("Must run PCA preparation before UMAP. Call run_pca_preparation() first.")
        
        logger.info("=" * 80)
        logger.info("STEP 2: UMAP Reduction")
        logger.info("=" * 80)
        
        # Initialize UMAP reducer
        self.umap_reducer = UMAPReduction(
            eigenvector=self.eigenvector_file,
            output_path=self.results_dir
        )
        
        # Run UMAP
        self.umap_results = self.umap_reducer.fit_transform(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state,
            n_components=n_components,
            umap_kwargs=umap_kwargs
        )
        
        # Save results
        output_file = self.results_dir / 'umap_results' / 'umap_coordinates.tsv'
        self.umap_results.to_csv(output_file, sep='\t', index=False)
        logger.info(f"UMAP results saved to {output_file}")
        
        return self.umap_results
    
    def execute_tsne(self,
                 n_components: int = 2,
                 perplexity: float = 30.0,
                 learning_rate: float = 200.0,
                 n_iter: int = 1000,
                 metric: str = 'euclidean',
                 random_state: Optional[int] = None,
                 early_exaggeration: float = 12.0,
                 init: Literal['pca', 'random'] = 'pca',
                 tsne_kwargs: Optional[dict] = None) -> pd.DataFrame:
        """
        Run t-SNE dimensionality reduction.
        
        Parameters
        ----------
        n_components : int, default=2
            Number of output dimensions
        perplexity : float, default=30.0
            t-SNE perplexity parameter
        learning_rate : float, default=200.0
            Learning rate for optimization
        n_iter : int, default=1000
            Number of optimization iterations
        metric : str, default='euclidean'
            Distance metric
        random_state : int, optional
            Random seed for reproducibility
        early_exaggeration : float, default=12.0
            Early exaggeration parameter
        init : {'pca', 'random'}, default='pca'
            Initialization method
        tsne_kwargs : dict, optional
            Additional t-SNE parameters
        
        Returns
        -------
        pd.DataFrame
            t-SNE results with columns ['ID1', 'ID2', 'tsne_1', 'tsne_2', ...]
        
        Raises
        ------
        RuntimeError
            If PCA preparation hasn't been run yet
        """
        if self.eigenvector_file is None:
            raise RuntimeError("Must run PCA preparation before t-SNE. Call run_pca_preparation() first.")
        
        logger.info("=" * 80)
        logger.info("STEP 3: t-SNE Reduction")
        logger.info("=" * 80)
        
        # Initialize t-SNE reducer
        self.tsne_reducer = TSNEReduction(
            eigenvector=self.eigenvector_file,
            output_path=self.results_dir
        )
        
        # Run t-SNE
        self.tsne_results = self.tsne_reducer.fit_transform(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=n_iter,
            metric=metric,
            random_state=random_state,
            early_exaggeration=early_exaggeration,
            init=init,
            tsne_kwargs=tsne_kwargs
        )
        
        # Save results
        output_file = self.results_dir / 'tsne_results' / 'tsne_coordinates.tsv'
        self.tsne_results.to_csv(output_file, sep='\t', index=False)
        logger.info(f"t-SNE results saved to {output_file}")
        
        return self.tsne_results
    
    def generate_plots(self,
                      color_hue_file: Optional[Path] = None,
                      case_control_markers: bool = False,
                      fam_file: Optional[Path] = None,
                      plot_format: str = 'pdf',
                      dpi: int = 500,
                      figsize: tuple = (5, 5),) -> dict:
        """
        Generate plots for all completed reductions.
        
        Parameters
        ----------
        color_hue_file : Path, optional
            Path to metadata file for coloring
        case_control_markers : bool, default=False
            Whether to use case-control markers
        fam_file : Path, optional
            Path to .fam file (required if case_control_markers=True)
        plot_format : str, default='pdf'
            Output format ('pdf', 'png', 'jpeg', 'svg')
        dpi : int, default=500
            Resolution for plots
        figsize : tuple, default=(5, 5)
            Figure size in inches
        include_pca : bool, default=True
            Whether to generate PCA plots
        
        Returns
        -------
        dict
            Dictionary mapping method names to plot file paths
        
        Raises
        ------
        RuntimeError
            If no reductions have been run
        """
        has_pca = hasattr(self, 'eigenvector_file') and self.eigenvector_file is not None and self.eigenvector_file.exists()
        if self.umap_results is None and self.tsne_results is None and not has_pca:
            raise RuntimeError("No reductions have been run. Run PCA preparation and/or UMAP/t-SNE first.")
        
        logger.info("=" * 80)
        logger.info("STEP 4: Generating Plots")
        logger.info("=" * 80)
        
        # Create plots directory
        plots_dir = self.results_dir / 'plots'
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize plotter
        self.plotter = Plot2D(output_dir=plots_dir)
        
        # Prepare metadata once
        self.plotter.prepare_metadata(
            color_hue_file=color_hue_file,
            case_control_markers=case_control_markers,
            fam_file=fam_file
        )
        
        plot_paths = {}
        
        # Generate UMAP plot
        if self.umap_results is not None:
            logger.info("Generating UMAP plot...")
            umap_path = self.plotter.generate_plot(
                data=self.umap_results,
                x_col='umap_1',
                y_col='umap_2',
                plot_name=f'umap_plot.{plot_format}',
                xlabel='UMAP 1',
                ylabel='UMAP 2',
                figsize=figsize,
                dpi=dpi,
                format=plot_format
            )
            plot_paths['umap'] = umap_path
        
        # Generate t-SNE plot
        if self.tsne_results is not None:
            logger.info("Generating t-SNE plot...")
            tsne_path = self.plotter.generate_plot(
                data=self.tsne_results,
                x_col='tsne_1',
                y_col='tsne_2',
                plot_name=f'tsne_plot.{plot_format}',
                xlabel='t-SNE 1',
                ylabel='t-SNE 2',
                figsize=figsize,
                dpi=dpi,
                format=plot_format
            )
            plot_paths['tsne'] = tsne_path
        
        logger.info("All plots generated successfully")
        return plot_paths
    
    def execute_dimensionality_reduction_pipeline(self,
                         # PCA parameters
                         pca_params: Optional[dict] = None,
                         force_pca_recompute: bool = False,
                         # UMAP parameters
                         run_umap: bool = True,
                         umap_params: Optional[dict] = None,
                         # t-SNE parameters
                         run_tsne: bool = True,
                         tsne_params: Optional[dict] = None,
                         # Plotting parameters
                         color_hue_file: Optional[Path] = None,
                         case_control_markers: bool = False,
                         fam_file: Optional[Path] = None,
                         plot_format: str = 'pdf',
                         dpi: int = 500,
                         include_pca: bool = True,
                         # Grid-specific parameters
                         save_all_coordinates: bool = True,
                         generate_all_plots: bool = True,
                         grid_summary: bool = True) -> dict:
        """
        Run the complete dimensionality reduction pipeline with automatic parameter grid detection.
        
        This method automatically detects whether parameters contain single values or lists.
        If lists are detected, it runs a parameter grid search exploring all combinations.
        Otherwise, it runs a single analysis with the provided parameters.
        
        Parameters
        ----------
        pca_params : dict, optional
            Parameters for PCA preparation (maf, geno, mind, hwe, ind_pair, pca)
        force_pca_recompute : bool, default=False
            If True, recompute PCA even if files already exist. If False, skip
            PCA computation if eigenvector and eigenvalue files are found.
        run_umap : bool, default=True
            Whether to run UMAP reduction
        umap_params : dict, optional
            Parameters for UMAP. Can contain single values or lists for grid search.
            Example single: {'n_neighbors': 15, 'min_dist': 0.1}
            Example grid: {'n_neighbors': [10, 15, 30], 'min_dist': [0.1, 0.5]}
        run_tsne : bool, default=True
            Whether to run t-SNE reduction
        tsne_params : dict, optional
            Parameters for t-SNE. Can contain single values or lists for grid search.
            Example single: {'perplexity': 30, 'learning_rate': 200}
            Example grid: {'perplexity': [20, 30, 50], 'learning_rate': [100, 200]}
        color_hue_file : Path, optional
            Metadata file for plot coloring
        case_control_markers : bool, default=False
            Whether to use case-control markers in plots
        fam_file : Path, optional
            Path to .fam file. If not provided, will automatically look for
            {input_name}.fam in the input_path directory
        plot_format : str, default='pdf'
            Output format for plots
        dpi : int, default=500
            Resolution for plots
        include_pca : bool, default=True
            Whether to generate PCA plots
        save_all_coordinates : bool, default=True
            For grid search: whether to save coordinate files for all parameter combinations
        generate_all_plots : bool, default=True
            For grid search: whether to generate plot files for all parameter combinations
        grid_summary : bool, default=True
            For grid search: whether to generate summary table of all parameter combinations
        
        Returns
        -------
        dict
            Results summary with file paths and metadata. For grid searches, includes
            information about all parameter combinations explored.
        
        Examples
        --------
        Single analysis:
        >>> pipeline = DimensionalityReductionPipeline(...)
        >>> results = pipeline.run_full_pipeline(
        ...     umap_params={'n_neighbors': 15, 'min_dist': 0.1},
        ...     tsne_params={'perplexity': 30}
        ... )
        
        Parameter grid search (automatically detected):
        >>> results = pipeline.run_full_pipeline(
        ...     umap_params={
        ...         'n_neighbors': [10, 15, 30], 
        ...         'min_dist': [0.1, 0.5],
        ...         'random_state': [42]
        ...     },
        ...     tsne_params={
        ...         'perplexity': [20, 30, 50],
        ...         'random_state': [42]
        ...     }
        ... )
        """
        logger.info("=" * 80)
        logger.info("STARTING FULL DIMENSIONALITY REDUCTION PIPELINE")
        logger.info("=" * 80)
        
        # Set default parameters
        pca_params = pca_params or {}
        umap_params = umap_params or {}
        tsne_params = tsne_params or {}
        
        # Check if this is a parameter grid search
        is_umap_grid = self._is_parameter_grid(umap_params) if run_umap else False
        is_tsne_grid = self._is_parameter_grid(tsne_params) if run_tsne else False
        
        if is_umap_grid or is_tsne_grid:
            logger.info("Parameter grid search detected - running grid exploration")
            return self._run_parameter_grid_pipeline(
                pca_params=pca_params,
                force_pca_recompute=force_pca_recompute,
                run_umap=run_umap,
                umap_params=umap_params if is_umap_grid else None,
                run_tsne=run_tsne, 
                tsne_params=tsne_params if is_tsne_grid else None,
                color_hue_file=color_hue_file,
                case_control_markers=case_control_markers,
                fam_file=fam_file,
                plot_format=plot_format,
                dpi=dpi,
                save_all_coordinates=save_all_coordinates,
                generate_all_plots=generate_all_plots,
                grid_summary=grid_summary
            )
        else:
            logger.info("Single parameter set detected - running standard pipeline")
            return self._run_single_parameter_pipeline(
                pca_params=pca_params,
                force_pca_recompute=force_pca_recompute,
                run_umap=run_umap,
                umap_params=umap_params,
                run_tsne=run_tsne,
                tsne_params=tsne_params,
                color_hue_file=color_hue_file,
                case_control_markers=case_control_markers,
                fam_file=fam_file,
                plot_format=plot_format,
                dpi=dpi,
                include_pca=include_pca
            )
    
    def _run_single_parameter_pipeline(self,
                                     pca_params: dict,
                                     force_pca_recompute: bool,
                                     run_umap: bool,
                                     umap_params: dict,
                                     run_tsne: bool,
                                     tsne_params: dict,
                                     color_hue_file: Optional[Path],
                                     case_control_markers: bool,
                                     fam_file: Optional[Path],
                                     plot_format: str,
                                     dpi: int,
                                     include_pca: bool) -> dict:
        """Run pipeline with single parameter sets (original behavior)."""
        
        results = {
            'input': self.input_name,
            'output_dir': str(self.results_dir),
            'steps_completed': [],
            'files': {}
        }
        
        # Step 1: PCA Preparation (always required)
        expected_eigenvec = self.results_dir / f"{self.input_name}-pca.eigenvec"
        expected_eigenval = self.results_dir / f"{self.input_name}-pca.eigenval"
        
        if not force_pca_recompute and expected_eigenvec.exists() and expected_eigenval.exists():
            logger.info("PCA files already exist, skipping PCA preparation")
            logger.info(f"Using existing eigenvector file: {expected_eigenvec}")
            logger.info(f"Using existing eigenvalue file: {expected_eigenval}")
            self.eigenvector_file = expected_eigenvec
            self.eigenvalue_file = expected_eigenval
            eigenvec_path = expected_eigenvec
        else:
            if force_pca_recompute:
                logger.info("Force recompute enabled, running PCA preparation")
            else:
                logger.info("PCA files not found, running PCA preparation")
            # Extract individual PCA parameters from pca_params dictionary
            pca_kwargs = {
                'maf': pca_params.get('maf', 0.001),
                'geno': pca_params.get('geno', 0.1),
                'mind': pca_params.get('mind', 0.2),
                'hwe': pca_params.get('hwe', 5e-8),
                'ind_pair': pca_params.get('ind_pair', [50, 5, 0.2]),
                'pca': pca_params.get('pca', 20),
                'case_control_markers': case_control_markers
            }
            eigenvec_path = self.execute_pca_preparation(**pca_kwargs)
            
        results['steps_completed'].append('pca_preparation')
        results['files']['eigenvector'] = str(eigenvec_path)
        results['files']['eigenvalue'] = str(self.eigenvalue_file)
        
        # Step 2: UMAP (optional)
        if run_umap:
            self.execute_umap(**umap_params)
            results['steps_completed'].append('umap')
            results['files']['umap_coordinates'] = str(
                self.results_dir / 'umap_results' / 'umap_coordinates.tsv'
            )

        # Step 3: t-SNE (optional)
        if run_tsne:
            self.execute_tsne(**tsne_params)
            results['steps_completed'].append('tsne')
            results['files']['tsne_coordinates'] = str(
                self.results_dir / 'tsne_results' / 'tsne_coordinates.tsv'
            )
        
        # Step 4: Generate plots
        has_pca = hasattr(self, 'eigenvector_file') and self.eigenvector_file is not None
        if run_umap or run_tsne or (include_pca and has_pca):
            # Auto-detect fam file from input path if not provided
            if fam_file is None:
                potential_fam_file = self.input_path / f"{self.input_name}.fam"
                if potential_fam_file.exists():
                    fam_file = potential_fam_file
                    logger.info(f"Auto-detected fam file: {fam_file}")
                elif case_control_markers:
                    logger.warning("case_control_markers=True but no fam file found. Plots will be generated without case-control markers.")
            
            plot_paths = self.generate_plots(
                color_hue_file=color_hue_file,
                case_control_markers=case_control_markers,
                fam_file=fam_file,
                plot_format=plot_format,
                dpi=dpi
            )
            results['steps_completed'].append('plotting')
            results['files']['plots'] = {k: str(v) for k, v in plot_paths.items()}
        
        logger.info("=" * 80)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY")
        logger.info("=" * 80)
        logger.info(f"Steps completed: {', '.join(results['steps_completed'])}")
        logger.info(f"Results saved to: {self.results_dir}")
        
        return results
    
    def _run_parameter_grid_pipeline(self,
                                   pca_params: dict,
                                   force_pca_recompute: bool,
                                   run_umap: bool,
                                   umap_params: Optional[dict],
                                   run_tsne: bool,
                                   tsne_params: Optional[dict],
                                   color_hue_file: Optional[Path],
                                   case_control_markers: bool,
                                   fam_file: Optional[Path],
                                   plot_format: str,
                                   dpi: int,
                                   save_all_coordinates: bool,
                                   generate_all_plots: bool,
                                   grid_summary: bool) -> dict:
        """Run pipeline with parameter grid search."""
        
        # First ensure PCA is completed
        expected_eigenvec = self.results_dir / f"{self.input_name}-pca.eigenvec"
        expected_eigenval = self.results_dir / f"{self.input_name}-pca.eigenval"
        
        if not force_pca_recompute and expected_eigenvec.exists() and expected_eigenval.exists():
            logger.info("PCA files already exist, skipping PCA preparation")
            logger.info(f"Using existing eigenvector file: {expected_eigenvec}")
            logger.info(f"Using existing eigenvalue file: {expected_eigenval}")
            self.eigenvector_file = expected_eigenvec
            self.eigenvalue_file = expected_eigenval
        else:
            if force_pca_recompute:
                logger.info("Force recompute enabled, running PCA preparation")
            else:
                logger.info("PCA files not found, running PCA preparation")
            # Extract individual PCA parameters from pca_params dictionary
            pca_kwargs = {
                'maf': pca_params.get('maf', 0.001),
                'geno': pca_params.get('geno', 0.1),
                'mind': pca_params.get('mind', 0.2),
                'hwe': pca_params.get('hwe', 5e-8),
                'ind_pair': pca_params.get('ind_pair', [50, 5, 0.2]),
                'pca': pca_params.get('pca', 20),
                'case_control_markers': case_control_markers
            }
            self.execute_pca_preparation(**pca_kwargs)
        
        # Now run the parameter grid
        return self.execute_parameter_grid(
            umap_grid=umap_params if (run_umap and umap_params) else None,
            tsne_grid=tsne_params if (run_tsne and tsne_params) else None,
            save_coordinates=save_all_coordinates,
            generate_plots=generate_all_plots,
            color_hue_file=color_hue_file,
            case_control_markers=case_control_markers,
            fam_file=fam_file,
            plot_format=plot_format
        )
    
    def execute_parameter_grid(self,
                          umap_grid: Optional[dict] = None,
                          tsne_grid: Optional[dict] = None,
                          plot_params: Optional[dict] = None,
                          save_coordinates: bool = True,
                          generate_plots: bool = True,
                          color_hue_file: Optional[Path] = None,
                          case_control_markers: bool = False,
                          fam_file: Optional[Path] = None,
                          plot_format: str = 'pdf') -> dict:
        """
        Run systematic parameter grid exploration for UMAP and/or t-SNE.
        
        This method explores all combinations of specified parameters, saving
        coordinates and generating plots for each combination. Results are
        organized with clear naming conventions for easy comparison.
        
        Parameters
        ----------
        umap_grid : dict, optional
            Dictionary with parameter names as keys and lists of values as values.
            Example: {'n_neighbors': [15, 30], 'min_dist': [0.1, 0.5]}
        tsne_grid : dict, optional
            Dictionary with parameter names as keys and lists of values as values.
            Example: {'perplexity': [20, 50], 'learning_rate': [100, 200]}
        plot_params : dict, optional
            Additional parameters for plot generation (figsize, dpi, etc.)
        save_coordinates : bool, default=True
            Whether to save coordinate files for each combination
        generate_plots : bool, default=True
            Whether to generate plot files for each combination
        color_hue_file : Path, optional
            Metadata file for plot coloring
        case_control_markers : bool, default=False
            Whether to use case-control markers in plots
        fam_file : Path, optional
            Path to .fam file for case-control markers
        plot_format : str, default='pdf'
            Output format for plots
        
        Returns
        -------
        dict
            Summary of all parameter combinations and results
            
        Raises
        ------
        RuntimeError
            If PCA preparation hasn't been run yet
        ValueError
            If neither umap_grid nor tsne_grid is provided
        
        Examples
        --------
        >>> pipeline = DimensionalityReductionPipeline(...)
        >>> pipeline.run_pca_preparation()
        >>> results = pipeline.run_parameter_grid(
        ...     umap_grid={
        ...         'n_neighbors': [15, 30], 
        ...         'min_dist': [0.1, 0.5],
        ...         'random_state': [42]
        ...     },
        ...     tsne_grid={
        ...         'perplexity': [20, 50],
        ...         'random_state': [42]
        ...     }
        ... )
        """
        import itertools
        
        if self.eigenvector_file is None:
            raise RuntimeError("Must run PCA preparation before parameter grid. Call run_pca_preparation() first.")
        
        if umap_grid is None and tsne_grid is None:
            raise ValueError("Must provide at least one of umap_grid or tsne_grid")
        
        logger.info("=" * 80)
        logger.info("STARTING PARAMETER GRID EXPLORATION")
        logger.info("=" * 80)
        
        # Set default plot parameters
        if plot_params is None:
            plot_params = {'figsize': (6, 6), 'dpi': 500}
        
        # Auto-detect fam file if needed
        if fam_file is None and case_control_markers:
            potential_fam_file = self.input_path / f"{self.input_name}.fam"
            if potential_fam_file.exists():
                fam_file = potential_fam_file
                logger.info(f"Auto-detected fam file: {fam_file}")
        
        results = {
            'umap_combinations': [],
            'tsne_combinations': [],
            'total_combinations': 0,
            'files_created': []
        }
        
        # UMAP parameter grid
        if umap_grid is not None:
            logger.info("Running UMAP parameter grid...")
            
            # Generate all UMAP parameter combinations
            param_names = list(umap_grid.keys())
            param_values = list(umap_grid.values())
            umap_combinations = list(itertools.product(*param_values))
            
            logger.info(f"UMAP: {len(umap_combinations)} parameter combinations to explore")
            
            # Create UMAP grid results directory
            umap_grid_dir = self.results_dir / 'umap_grid_results'
            umap_grid_dir.mkdir(parents=True, exist_ok=True)
            
            for i, combination in enumerate(umap_combinations):
                # Create parameter dictionary for this combination
                params = dict(zip(param_names, combination))
                
                # Create readable parameter string for naming
                param_str = "_".join([f"{k}{v}" for k, v in params.items() if k != 'random_state'])
                config_name = f"umap_{param_str}"
                
                logger.info(f"UMAP {i+1}/{len(umap_combinations)}: {config_name}")
                
                # Run UMAP with these parameters
                umap_result = self.execute_umap(**params)
                
                # Save coordinates if requested
                if save_coordinates:
                    coords_file = umap_grid_dir / f'{config_name}_coordinates.tsv'
                    umap_result.to_csv(coords_file, sep='\t', index=False)
                    results['files_created'].append(str(coords_file))
                
                # Generate plot if requested
                if generate_plots:
                    # Initialize plotter for this combination
                    plotter = Plot2D(output_dir=umap_grid_dir)
                    plotter.prepare_metadata(
                        color_hue_file=color_hue_file,
                        case_control_markers=case_control_markers,
                        fam_file=fam_file
                    )
                    
                    #plot_title = f"UMAP: {', '.join([f'{k}={v}' for k, v in params.items() if k != 'random_state'])}"
                    plot_name = f'{config_name}_plot.{plot_format}'
                    
                    plot_path = plotter.generate_plot(
                        data=umap_result,
                        x_col='umap_1',
                        y_col='umap_2',
                        plot_name=plot_name,
                        #title=plot_title,
                        xlabel='UMAP 1',
                        ylabel='UMAP 2',
                        format=plot_format,
                        **plot_params
                    )
                    results['files_created'].append(str(plot_path))
                
                # Store combination info
                combination_info = params.copy()
                combination_info.update({
                    'config_name': config_name,
                    'n_samples': len(umap_result)
                })
                results['umap_combinations'].append(combination_info)
        
        # t-SNE parameter grid
        if tsne_grid is not None:
            logger.info("Running t-SNE parameter grid...")
            
            # Generate all t-SNE parameter combinations
            param_names = list(tsne_grid.keys())
            param_values = list(tsne_grid.values())
            tsne_combinations = list(itertools.product(*param_values))
            
            logger.info(f"t-SNE: {len(tsne_combinations)} parameter combinations to explore")
            
            # Create t-SNE grid results directory
            tsne_grid_dir = self.results_dir / 'tsne_grid_results'
            tsne_grid_dir.mkdir(parents=True, exist_ok=True)
            
            for i, combination in enumerate(tsne_combinations):
                # Create parameter dictionary for this combination
                params = dict(zip(param_names, combination))
                
                # Create readable parameter string for naming
                param_str = "_".join([f"{k}{v}" for k, v in params.items() if k != 'random_state'])
                config_name = f"tsne_{param_str}"
                
                logger.info(f"t-SNE {i+1}/{len(tsne_combinations)}: {config_name}")
                
                # Run t-SNE with these parameters
                tsne_result = self.execute_tsne(**params)
                
                # Save coordinates if requested
                if save_coordinates:
                    coords_file = tsne_grid_dir / f'{config_name}_coordinates.tsv'
                    tsne_result.to_csv(coords_file, sep='\t', index=False)
                    results['files_created'].append(str(coords_file))
                
                # Generate plot if requested
                if generate_plots:
                    # Initialize plotter for this combination
                    plotter = Plot2D(output_dir=tsne_grid_dir)
                    plotter.prepare_metadata(
                        color_hue_file=color_hue_file,
                        case_control_markers=case_control_markers,
                        fam_file=fam_file
                    )
                    
                    plot_title = f"t-SNE: {', '.join([f'{k}={v}' for k, v in params.items() if k != 'random_state'])}"
                    plot_name = f'{config_name}_plot.{plot_format}'
                    
                    plot_path = plotter.generate_plot(
                        data=tsne_result,
                        x_col='tsne_1',
                        y_col='tsne_2',
                        plot_name=plot_name,
                        title=plot_title,
                        xlabel='t-SNE 1',
                        ylabel='t-SNE 2',
                        format=plot_format,
                        **plot_params
                    )
                    results['files_created'].append(str(plot_path))
                
                # Store combination info
                combination_info = params.copy()
                combination_info.update({
                    'config_name': config_name,
                    'n_samples': len(tsne_result)
                })
                results['tsne_combinations'].append(combination_info)
        
        # Calculate total combinations
        results['total_combinations'] = len(results['umap_combinations']) + len(results['tsne_combinations'])
        
        # Save method-specific summaries to avoid cross-method NaN columns
        if results['umap_combinations']:
            umap_df = pd.DataFrame(results['umap_combinations'])
            umap_df.insert(0, 'method', 'UMAP')
            umap_summary_file = self.results_dir / 'umap_grid_summary.tsv'
            umap_df.to_csv(umap_summary_file, sep='\t', index=False)
            results['files_created'].append(str(umap_summary_file))
            logger.info(f"UMAP grid summary saved to: {umap_summary_file}")

        if results['tsne_combinations']:
            tsne_df = pd.DataFrame(results['tsne_combinations'])
            tsne_df.insert(0, 'method', 't-SNE')
            tsne_summary_file = self.results_dir / 'tsne_grid_summary.tsv'
            tsne_df.to_csv(tsne_summary_file, sep='\t', index=False)
            results['files_created'].append(str(tsne_summary_file))
            logger.info(f"t-SNE grid summary saved to: {tsne_summary_file}")
        
        logger.info("=" * 80)
        logger.info("PARAMETER GRID EXPLORATION COMPLETED")
        logger.info("=" * 80)
        logger.info(f"Total combinations explored: {results['total_combinations']}")
        logger.info(f"Files created: {len(results['files_created'])}")
        logger.info(f"Results saved to: {self.results_dir}")
        
        return results