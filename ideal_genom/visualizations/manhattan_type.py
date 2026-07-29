"""Generate Manhattan, Miami and Brisbane plots for genome-wide association studies (GWAS).

This module provides functions to create Manhattan, Miami and Brisbane plots to visualize
the significance and density of genetic variants across the genome.

Features:
- Data processing and visualization of GWAS summary statistics.
- Annotation of SNPs with gene information from various sources.
- Highlighting and labeling of specific SNPs of interest.

"""

import logging
import os
import warnings
from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
import textalloc as ta
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from ..utilities.annotations import annotate_snp

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

_NON_AUTOSOME_ORDER = {'X': 23, 'Y': 24, 'XY': 25, 'MT': 26, 'M': 26}


def _chromosome_sort_key(chrom) -> tuple:
    """Map a chromosome label to a sort key reflecting standard genomic order (1-22, X, Y, XY, MT/M).

    Accepts ints, numeric strings, and labels with an optional leading 'chr' (case-insensitive).
    Unrecognized labels sort after all known chromosomes, ordered alphabetically among themselves.
    """

    label = str(chrom).strip()
    if label.lower().startswith('chr'):
        label = label[3:]
    label_upper = label.upper()

    if label_upper in _NON_AUTOSOME_ORDER:
        return (0, _NON_AUTOSOME_ORDER[label_upper], '')

    try:
        return (0, int(label), '')
    except ValueError:
        return (1, 0, label_upper)


def _compute_chr_offsets(data: pd.DataFrame, chr_col: str, pos_col: str) -> pd.DataFrame:
    """Compute, per chromosome, the cumulative genomic offset in standard genomic order.

    Returns a DataFrame with columns [chr_col, 'cumulativechrlength'], where
    'cumulativechrlength' is the sum of the lengths of all chromosomes ordered before it
    (1-22, X, Y, XY, MT/M, then anything unrecognized alphabetically).
    """

    chr_grouped = data.groupby(chr_col).agg(chrlength=(pos_col, 'max')).reset_index()
    chr_grouped = chr_grouped.sort_values(
        by=chr_col, key=lambda s: s.map(_chromosome_sort_key)
    ).reset_index(drop=True)

    chr_grouped['cumulativechrlength'] = chr_grouped['chrlength'].cumsum() - chr_grouped['chrlength']

    return chr_grouped[[chr_col, 'cumulativechrlength']]


def compute_relative_pos(
    data: pd.DataFrame, 
    chr_col: str = 'CHR', 
    pos_col: str = 'POS', 
    p_col: str = 'p'
) -> pd.DataFrame:
    
    """Compute the relative position of probes/SNPs across chromosomes and add a -log10(p-value) column.

    Parameters
    ----------
    data : pandas.DataFrame 
        Input DataFrame containing genomic data.
    chr_col : str
        Column name for chromosome identifiers. Default is 'CHR'.
    pos_col : str 
        Column name for base pair positions. Default is 'POS'.
    p_col : str 
        Column name for p-values. Default is 'p'.
    
    Returns
    -------
    pandas.DataFrame 
        DataFrame with additional columns for relative positions and -log10(p-values).

    Raises
    ------
    TypeError 
        If input data is not a pandas DataFrame.
    ValueError 
        If `chr_col`, `pos_col` or `p_col` columns are not found in the DataFrame.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

    # Compute per-chromosome cumulative offsets in genomic order (not lexicographic,
    # which would e.g. place chr10 before chr2), and merge back onto the original data
    chr_offsets = _compute_chr_offsets(data, chr_col=chr_col, pos_col=pos_col)
    data = pd.merge(data, chr_offsets, on=chr_col)

    # Sort by chromosome (genomic order) and position
    data = data.sort_values(
        by=[chr_col, pos_col], key=lambda s: s.map(_chromosome_sort_key) if s.name == chr_col else s
    )

    # Add the relative position of the probe/snp
    data['rel_pos'] = data[pos_col] + data['cumulativechrlength']

    # Drop cumulative chromosome length column
    data = data.drop(columns=['cumulativechrlength'])

    data['log10p']= -np.log10(data[p_col])

    return data

def find_chromosomes_center(
    data: pd.DataFrame, 
    chr_col: str = 'CHR', 
    chr_pos_col: str = 'rel_pos'
) -> pd.DataFrame:
    
    """Calculate the center positions of chromosomes in a given DataFrame. 
    
    This function takes a DataFrame containing chromosome data and calculates the center position 
    for each chromosome based on the specified chromosome column and chromosome position column.
    
    Parameters
    ----------
    data : pandas.DataFrame
        The input DataFrame containing chromosome data.
    chr_col : str, optional
        The name of the column representing chromosome identifiers (default is 'CHR').
    chr_pos_col : str, optional
        The name of the column representing relative positions within chromosomes (default is 'rel_pos').
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns 'CHR' and 'center', where 'CHR' contains chromosome identifiers 
        and 'center' contains the calculated center positions for each chromosome.

    Raises
    ------
    TypeError
        If the input data is not a pandas DataFrame.
    ValueError
        If the specified chromosome or chromosome position columns are not found in the DataFrame.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if chr_pos_col not in data.columns:
        raise ValueError(f"Column '{chr_pos_col}' not found in the input DataFrame.")

    axis_center = (
        data.groupby(chr_col)[chr_pos_col]
        .agg(lambda s: np.round((s.max() + s.min()) / 2, 0))
        .rename('center')
        .reset_index()
    )

    return axis_center

def manhattan_process_data(
    data_df: pd.DataFrame, 
    chr_col: str = 'CHR', 
    pos_col: str = 'POS', 
    p_col: str = 'p'
) -> dict:
    
    """Processes the input DataFrame to prepare data for a Manhattan plot.

    The function computes the relative positions of SNPs across chromosomes,
    calculates the -log10(p-values), and finds the center positions of each chromosome.

    Parameters
    ----------
    data_df : pandas.DataFrame
        The input DataFrame containing genomic data.
    chr_col : str (optional)
        The column name for chromosome data. Defaults to 'CHR'.
    pos_col : str (optional)
        The column name for position data. Defaults to 'POS'.
    p_col : str (optional)
        The column name for p-value data. Defaults to 'p'.

    Raises
    ------
    TypeError
        If `data_df` is not a pandas DataFrame.
        If `chr_col`, `pos_col`, or `p_col` are not strings.
    ValueError
        If the specified columns (`chr_col`, `pos_col`, `p_col`) are not found in the DataFrame.

    Returns
    -------
    dict
        A dictionary containing processed data for the Manhattan plot with the following keys:

        - data : pandas.DataFrame
            The processed DataFrame with relative positions and log-transformed p-values.
        - axis : dict
            The center positions of each chromosome for plotting.
        - maxp : float
            The maximum log-transformed p-value.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if not isinstance(chr_col, str):
        raise TypeError("chr_col must be a string.")
    if not isinstance(pos_col, str):
        raise TypeError("pos_col must be a string.")
    if not isinstance(p_col, str):
        raise TypeError("p_col must be a string.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

    data = compute_relative_pos(
        data_df, 
        chr_col=chr_col, 
        pos_col=pos_col, 
        p_col  =p_col
    )

    axis_center = find_chromosomes_center(data, chr_col=chr_col)

    maxp = np.ceil(data['log10p'].max(skipna=True))

    manhattan_data = {
        'data': data.reset_index(drop=True),
        'axis': axis_center,
        'maxp': maxp
    }

    return manhattan_data

def manhattan_draw(
    data_df: pd.DataFrame, 
    snp_col: str, 
    chr_col: str, 
    pos_col: str, 
    p_col: str, 
    plot_dir: str, 
    to_highlight: pd.DataFrame = pd.DataFrame(), 
    highlight_hue: str = 'hue', 
    to_annotate: pd.DataFrame = pd.DataFrame(), 
    gen_col: Optional[str] = None, 
    build: str = '38', 
    anno_source = 'ensembl', 
    gtf_path: Optional[str] = None, 
    save_name: str = 'manhattan_plot.png', 
    upper_cap: Optional[float] = None, 
    genome_line: float = 5e-8, 
    suggestive_line: float = 1e-5, 
    yaxis_margin: float = 10, 
    dpi: int = 400
) -> bool:
    
    """Generate a Manhattan plot for genomic data.
    
    This function creates a Manhattan plot for visualizing the statistical significance of genetic
    variants across the genome, typically used in Genome-Wide Association Studies (GWAS).
    The plot shows -log10(p-values) against genomic positions, organized by chromosome.
    
    Parameters
    ----------
    data_df : pandas.DataFrame
        DataFrame containing the SNP data to be plotted.
    snp_col : str
        Column name in data_df that contains the SNP identifiers.
    chr_col : str
        Column name in data_df that contains the chromosome information.
    pos_col : str
        Column name in data_df that contains the position information.
    p_col : str
        Column name in data_df that contains the p-values.
    plot_dir : str
        Directory path where the plot will be saved.
    to_highlight : pandas.DataFrame, optional
        DataFrame containing SNPs to be highlighted in the plot.
    highlight_hue : str, default='hue'
        Column name in to_highlight to be used for color-coding highlighted SNPs.
    to_annotate : pandas.DataFrame, optional
        DataFrame containing SNPs to be annotated in the plot.
    gen_col : str, optional
        Column name for gene information to be used in annotation.
    build : str, default='38'
        Genome build version ('38' for GRCh38, etc.).
    anno_source : str, default='ensembl'
        Source for gene annotation information.
    gtf_path : str, optional
        Path to a GTF file for gene annotation if not using online sources.
    save_name : str, default='manhattan_plot.png'
        Filename for the saved plot.
    upper_cap : float, optional
        Upper limit for log10(p-values) to be capped at.
    genome_line : float, default=5e-8
        P-value threshold for genome-wide significance line.
    suggestive_line : float, default=1e-5
        P-value threshold for suggestive significance line.
    yaxis_margin : float, default=10
        Additional margin to add to the y-axis maximum.
    dpi : int, default=400
        Resolution of the saved image in dots per inch.
    
    Returns
    -------
    bool
        True if the plot was successfully created and saved.
    
    Raises
    ------
    TypeError
        If input data is not a pandas DataFrame or to_highlight/to_annotate are not of correct type.
    ValueError
        If required columns are not found in the input DataFrame.
    FileNotFoundError
        If the specified plot directory does not exist.

    Notes
    -----
    Saves a Manhattan plot image to the specified directory.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if snp_col not in data_df.columns:
        raise ValueError(f"Column '{snp_col}' not found in the input DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_df.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

    if to_highlight is not None:
        if not isinstance(to_highlight, pd.DataFrame):
            raise TypeError("to_highlight must be a list of SNP identifiers.")
    if to_annotate is not None:
        if not isinstance(to_annotate, pd.DataFrame):
            raise TypeError("to_annotate must be a data frame of SNP identifiers.")
    if not isinstance(upper_cap, (type(None), float, int)):
        raise TypeError("upper_cap must be a float or None.")
    
    if not os.path.exists(plot_dir):
        raise FileNotFoundError(f"Directory '{plot_dir}' not found.")
    
    chr_colors           = ['#66c2a5', '#fc8d62']
    ylab                 = "-log10(p)"
    xlab                 = "Chromosome"
    genome_line_color    = "#636363"
    suggestive_line_color= "#377eb8"

    # format data to draw manhattan plot
    plot_data = manhattan_process_data(
        data_df=data_df,
        chr_col=chr_col,
        pos_col=pos_col,
        p_col  =p_col
    )

    if upper_cap:
        # Cap the log10(p) values at the specified upper limit
        plot_data['data'] = plot_data['data'][plot_data['data']['log10p'] <= upper_cap].reset_index(drop=True)
        plot_data['maxp'] = np.ceil(upper_cap)

    max_x_axis = plot_data['data']['rel_pos'].max()

    # Create the figure
    fig= plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111)

    # Suppress warnings about the number of chromosomes and just two colors
    warnings.filterwarnings("ignore", category=UserWarning)

    ax = sns.scatterplot(
        x      =plot_data['data']['rel_pos'], 
        y      =plot_data['data']['log10p'],
        hue    =plot_data['data'][chr_col], 
        palette=chr_colors, 
        ax     =ax, 
        s      =3, 
        legend =False,
        edgecolor='none'
    )

    # set axis labels and font size
    ax.set_ylabel(ylab, fontsize=7)
    ax.set_xlabel(xlab, fontsize=7)

    # set axis limits
    ax.set_xlim(0, max_x_axis+1000)
    ax.set_ylim(0, plot_data['maxp'] + yaxis_margin)

    # set x-axis ticks and labels
    x_ticks=plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis'][chr_col].astype(str).tolist()

    ax.set_xticks(ticks=x_ticks)
    ax.set_xticklabels(x_labels)

    # set ticks font size
    ax.tick_params(axis='both', labelsize=7)

    # add suggestive line
    if suggestive_line is not None:
        ax.axhline(
            -np.log10(suggestive_line), 
            color    =suggestive_line_color, 
            linestyle='dashed', 
            lw       =0.5
        )
    
    # add genome-wide line
    if genome_line is not None:
        ax.axhline(
            -np.log10(genome_line), 
            color    =genome_line_color, 
            linestyle='solid', 
            lw       =0.5
        )

    # highlight SNPs
    if to_highlight is not None:

        if to_highlight.empty is not True: 

            df_subset = plot_data['data'].merge(to_highlight, on=snp_col, how='inner')

            if df_subset.empty is not True:
            
                if highlight_hue in df_subset.columns:
                    ax = sns.scatterplot(
                        data     =df_subset, 
                        x        ='rel_pos',
                        y        ='log10p',
                        s        =20,
                        ax       =ax,
                        edgecolor="black",
                        hue      =highlight_hue,
                        palette='Set1',
                    )
                    ax.legend(loc='best', fontsize=9, title=None)
                else:
                    ax = sns.scatterplot(
                            data     =df_subset, 
                            x        ='rel_pos',
                            y        ='log10p',
                            s        =20,
                            ax       =ax,
                            edgecolor="black"
                        )
                
    plt.tight_layout()

    # annotate SNPs   
    if to_annotate is not None and to_annotate.empty is not True:

        if gen_col is not None:

            variants_toanno = plot_data['data'][plot_data['data'][snp_col].isin(to_annotate[snp_col])]\
                .reset_index(drop=True)
            variants_toanno = variants_toanno.merge(to_annotate, on=snp_col, how='left')
            variants_toanno = variants_toanno.rename(columns={gen_col:"GENENAME"})

        else:

            variants_toanno = plot_data['data'][plot_data['data'][snp_col].isin(to_annotate[snp_col])]\
                .reset_index(drop=True)
            
            if (variants_toanno.empty is not True):

                logger.info(f"Type of gtf_path is {type(gtf_path)}")
            
                variants_toanno = annotate_snp(
                    insumstats=variants_toanno,
                    chrom     =chr_col,
                    pos       =pos_col,
                    build     =build,
                    source    =anno_source,
                    gtf_path  =gtf_path  # type: ignore
                ).rename(columns={"GENE":"GENENAME"})

        ax = manhattan_type_annotate(
            axes           =ax, 
            data           =plot_data['data'], 
            variants_toanno=variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )

    # save the plot

    plt.savefig(
        os.path.join(plot_dir, save_name), dpi=dpi
    )
    plt.show()

    return True

def miami_process_data(
    data_top: pd.DataFrame, 
    data_bottom: pd.DataFrame, 
    chr_col: str, 
    pos_col: str, 
    p_col: str
) -> dict:
    
    """Processes Miami plot data by preparing, computing relative positions, and splitting the data.

    For each part of the Miami plot (top and bottom), this function computes the relative positions of SNPs,
    calculates the -log10(p-values), and finds the center positions of chromosomes.

    Parameters
    ----------
    data_top : pandas.DataFrame
        The top part of the data to be processed.
    data_bottom : pandas.DataFrame 
        The bottom part of the data to be processed.
    
    Returns
    -------
    dict 
        A dictionary containing the processed data with the following keys:
            - 'upper' : pandas.DataFrame 
                DataFrame containing the top part of the processed data.
            - 'lower' : pandas.DataFrame 
                DataFrame containing the bottom part of the processed data.
            - 'axis' : pandas.DataFrame
                DataFrame with the center positions of the chromosomes.
            - 'maxp' : float
                The maximum -log10(p-value) in the data.

    Raises
    ------
    TypeError
        If `data_top` or `data_bottom` is not a pandas DataFrame.
    ValueError
        If the specified columns (`chr_col`, `pos_col`, `p_col`) are not found in the DataFrames.
    """

    if not isinstance(data_top, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if not isinstance(data_bottom, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_top.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_top.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in data_top.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")

    data_top['split_by']   = 'top'
    data_bottom['split_by']= 'bottom'

    data = pd.concat([data_top, data_bottom], axis=0, ignore_index=True)

    data = compute_relative_pos(data, chr_col=chr_col, pos_col=pos_col, p_col=p_col)

    axis_center = find_chromosomes_center(data, chr_col=chr_col, chr_pos_col='rel_pos')

    maxp = np.ceil(data['log10p'].max(skipna=True))

    miami_data = {
        'upper': data[data['split_by'] == 'top'].reset_index(drop=True),
        'lower': data[data['split_by'] == 'bottom'].reset_index(drop=True),
        'axis': axis_center,
        'maxp': maxp
    }

    return miami_data

def manhattan_type_annotate(
    axes: Axes, data: pd.DataFrame, 
    variants_toanno: pd.DataFrame, 
    max_x_axis: float, 
    suggestive_line: float, 
    genome_line: float
) -> Axes:
    
    """Annotates a Manhattan plot with gene names.

    This function uses the textalloc library to place gene names on a Manhattan plot.

    Parameters
    ----------
    axes : Axes (matplotlib.axes.Axes)
        The matplotlib axes object where the Manhattan plot is drawn.
    data : pandas.DataFrame
        DataFrame containing the scatter plot data with columns 'rel_pos' and 'log10p'.
    variants_toanno : pandas.DataFrame
        DataFrame containing the variants to annotate with columns 'rel_pos', 'log10p', and 'GENENAME'.
    max_x_axis : float
        The maximum value for the x-axis.
    suggestive_line : float
        The y-value for the suggestive significance line.
    genome_line : float
        The y-value for the genome-wide significance line.

    Returns
    -------
    matplotlib.axes.Axes.
        The matplotlib axes object with annotations.

    Raises
    ------
    TypeError
        If the input parameters are not of the expected types.
    """

    if not isinstance(axes, Axes):
        raise TypeError("axes must be a matplotlib Axes object.")
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if not isinstance(variants_toanno, pd.DataFrame):
        raise TypeError("variants_toanno must be a pandas DataFrame.")
    if not isinstance(max_x_axis, (int, float, np.integer, np.floating)):
        raise TypeError("max_x_axis must be an integer or float.")
    if not isinstance(suggestive_line, float):
        raise TypeError("suggestive_line must be a float.")
    if not isinstance(genome_line, float):
        raise TypeError("genome_line must be a float.")
    
    x_lines_coor = np.linspace(0, max_x_axis, 1000).tolist() # list with a grid of x-coordinates for the lines
    
    x = variants_toanno['rel_pos'].tolist()      # x-coordinates for adjustment
    y = variants_toanno['log10p'].tolist()        # y-coordinates for adjustment
    texts = variants_toanno['GENENAME'].tolist()  # text annotations for adjustment

    allocate = ta.allocate(
            axes,            
            x        =x,     # x-coordinates of the data point to annotate
            y        =y,     # y-coordinates of the data point to annotate
            text_list=texts, # list of text to annotate
            x_scatter=data['rel_pos'].to_numpy(), # all scatter points x-coordinates
            y_scatter=data['log10p'].to_numpy(),  # all scatter points y-coordinates
            linecolor='black',                      # color of the line connecting the text to the data point
            textsize =10,                            # size of the text (Default to Nature standard)
            bbox     =dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='#f0f0f0', alpha=0.5),
            x_lines  = [x_lines_coor, x_lines_coor],
            y_lines  = [[-np.log10(suggestive_line)]*len(x_lines_coor), [-np.log10(genome_line)]*len(x_lines_coor)],
            avoid_label_lines_overlap =True,
            avoid_crossing_label_lines=True,
            min_distance=0.01,
            max_distance=0.4,
            margin      =0.01,
            rotation    =90,
            nbr_candidates=300,
            priority_strategy=42,
            plot_kwargs=dict(linestyle=':')
        )
        
    return axes

def miami_draw_anno_lines(renderer: RendererBase, axes: Axes, texts: list, variants_toanno: pd.DataFrame):
    
    """Draws annotation lines from text labels to their corresponding data points on a plot.

    Parameters
    ----------
    renderer : RendererBase 
        The renderer used to draw the plot.
    axes : Axes (matplotlib.axes.Axes) 
        The axes on which the plot is drawn.
    texts : list
        A list of text objects to annotate.
    variants_toanno : pandas.DataFrame
        A DataFrame containing the data points to annotate, with columns 'GENENAME', 'rel_pos', and 'log10p'.

    Returns
    -------
    matplotlib.axes.Axes 
        The axes with the annotation lines drawn.
    """

    from math import dist

    for k in range(len(texts)):

        text_obj = texts[k]
        bbox = text_obj.get_window_extent(renderer=renderer)
        
        # Transform to data coordinates the corners of the text bbox
        data_coords = axes.transData.inverted().transform(bbox.corners()) # list with coordinates of the text box

        anno = text_obj.get_text()

        data_x = variants_toanno.loc[variants_toanno['GENENAME'] == anno, 'rel_pos'].values[0]  # type: ignore
        data_y = variants_toanno.loc[variants_toanno['GENENAME'] == anno, 'log10p'].values[0] # type: ignore

        closest_point = data_coords[0]
        for point in data_coords:
            if dist([data_x, data_y], point) < dist([data_x, data_y], closest_point):
                closest_point = point 

        axes.plot(
            [data_x, closest_point[0]], 
            [data_y, closest_point[1]], 
            color='black', 
            linestyle=':', 
            linewidth=1
        )

    return axes

def miami_draw(
    df_top: pd.DataFrame, 
    df_bottom: pd.DataFrame, 
    snp_col: str, 
    chr_col: str, 
    pos_col: str, 
    p_col: str, 
    plots_dir: str, 
    top_highlights: list = [], 
    top_annotations: pd.DataFrame=pd.DataFrame(), 
    bottom_highlights: list = [], 
    bottom_annotations: pd.DataFrame = pd.DataFrame(), 
    top_gen_col: Optional[str] = None, 
    bottom_gen_col: Optional[str] = None, 
    gtf_path: Optional[str] = None, 
    source: str = "ensemble", 
    build: str = '38', 
    save_name: str = 'miami_plot.png', 
    legend_top: str = 'top GWAS', 
    legend_bottom: str = 'bottom GWAS', 
    dpi: int = 400
) -> bool:
    
    """Draws a Miami plot (a combination of two Manhattan plots) for visualizing GWAS results.

    Parameters
    ----------
    df_top : pandas.DataFrame
        DataFrame containing the top plot data.
    df_bottom : pandas.DataFrame
        DataFrame containing the bottom plot data.
    snp_col : str
        Column name for SNP identifiers.
    chr_col : str
        Column name for chromosome identifiers.
    pos_col : str
        Column name for base pair positions.
    p_col : str
        Column name for p-values.
    plots_dir : str
        Directory where the plot will be saved.
    top_highlights : list, optional
        List of SNPs to highlight in the top plot.
    top_annotations : list, optional
        List of SNPs to annotate in the top plot.
    bottom_highlights : list, optional
        List of SNPs to highlight in the bottom plot.
    bottom_annotations : list, optional
        List of SNPs to annotate in the bottom plot.
    gtf_path : str, optional
        Path to the GTF file for gene annotation. If None, the file will be downloaded.
    save_name : str, optional
        Name of the file to save the plot as. Default is 'miami_plot.jpeg'.

    Returns
    -------
    bool
        True if the plot is successfully created and saved, False otherwise.

    Raises
    ------
    TypeError
        If `df_top`, `df_bottom`, `top_annotations`, or `bottom_annotations` are not pandas DataFrames.
    TypeError
        If `top_highlights` or `bottom_highlights` are not lists of SNP identifiers.
    TypeError
        If `save_name`, `legend_top`, or `legend_bottom` are not strings.
    ValueError
        If required columns (`chr_col`, `pos_col`, `p_col`) are not found in the DataFrames.

    Notes
    -----
    Saves a Miami plot image to the specified directory.
    """
    
    if not isinstance(df_top, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if not isinstance(df_bottom, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in df_top.columns or chr_col not in df_bottom.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in df_top.columns or pos_col not in df_bottom.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if p_col not in df_top.columns or p_col not in df_bottom.columns:
        raise ValueError(f"Column '{p_col}' not found in the input DataFrame.")
    if not os.path.exists(plots_dir):
        raise FileNotFoundError(f"Directory '{plots_dir}' not found.")
    
    annot_high = [top_highlights,  bottom_highlights]
    for lst in annot_high:
        if not isinstance(lst, list):
            raise TypeError("Annotation lists must be lists of SNP identifiers.")
        for val in lst:
            if not isinstance(val, str):
                raise TypeError("Annotation lists must be lists of SNP identifiers.")
    
    if not isinstance(top_annotations, pd.DataFrame):
        raise TypeError("Annotation data must be a pandas DataFrame.")
    if not isinstance(bottom_annotations, pd.DataFrame):
        raise TypeError("Annotation data must be a pandas DataFrame.")
    
    if not isinstance(save_name, str):
        raise TypeError("save_name must be a string.")
    if not isinstance(legend_top, str):
        raise TypeError("legend_top must be a string.")
    if not isinstance(legend_bottom, str):
        raise TypeError("legend_bottom must be a string.")

    chr_colors           = ['#66c2a5', '#fc8d62']
    upper_ylab           = "-log10(p)" 
    lower_ylab           = "-log10(p)" 
    genome_line          = 5e-8
    genome_line_color    = "red"
    suggestive_line      = 1e-5 
    suggestive_line_color= "blue"

    # Suppress warnings about the number of chromosomes and just two colors
    warnings.filterwarnings("ignore", category=UserWarning)

    # format data to draw miami plot
    plot_data = miami_process_data(df_top, df_bottom, chr_col=chr_col, pos_col=pos_col, p_col=p_col)

    # Set axis labels for upper and lower plot
    def format_ylabel(label):
        return f"{label}\n-log10(p)" if label != "-log10(p)" else r"-log10(p)"
    
    upper_ylab = format_ylabel(upper_ylab)
    lower_ylab = format_ylabel(lower_ylab)

    max_x_axis = max(plot_data['upper']['rel_pos'].max(), plot_data['lower']['rel_pos'].max())+10

    # Create the figure
    fig = plt.figure(figsize=(20, 13))

    # Create the upper plot

    ax_upper = plt.subplot(211)

    sns.scatterplot(
        x      =plot_data['upper']['rel_pos'], 
        y      =plot_data['upper']['log10p'],
        hue    =plot_data['upper'][chr_col], 
        palette=chr_colors, 
        ax     =ax_upper, 
        s      =3, 
        legend =False,
        edgecolor='none'
    )
    ax_upper.set_ylabel(upper_ylab)
    ax_upper.set_xlim(0, max_x_axis)

    x_ticks =plot_data['axis']['center'].tolist()
    x_labels=plot_data['axis'][chr_col].astype(str).tolist()

    ax_upper.set_xticks(ticks=x_ticks)  # Set x-ticks
    ax_upper.set_xticklabels(x_labels)
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_upper.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.7)
    
    if genome_line is not None:
        ax_upper.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.7)

    # Create the lower plot
    ax_lower = plt.subplot(212)
    sns.scatterplot(
        x      =plot_data['lower']['rel_pos'], 
        y      =plot_data['lower']['log10p'],
        hue    =plot_data['lower'][chr_col], 
        palette=chr_colors, 
        ax     =ax_lower, 
        s      =3, 
        legend =False,
        edgecolor='none'
    )
    ax_lower.set_ylabel(lower_ylab)
    ax_lower.set_ylim(plot_data['maxp'], 0)  # Reverse y-axis
    ax_lower.set_xlim(0, max_x_axis)

    ax_lower.set_xticks(ticks=x_ticks) # Set x-ticks
    ax_lower.set_xticklabels([])
    ax_lower.xaxis.set_ticks_position('top')
    
    # Add genome-wide and suggestive lines
    if suggestive_line is not None:
        ax_lower.axhline(-np.log10(suggestive_line), color=suggestive_line_color, linestyle='solid', lw=0.7)
    
    if genome_line is not None:
        ax_lower.axhline(-np.log10(genome_line), color=genome_line_color, linestyle='dashed', lw=0.7)
    
    top = set(top_highlights)
    bottom = set(bottom_highlights)

    both = []
    top_only = []
    bottom_only = []

    if len(top_highlights)>0 or len(bottom_highlights)>0:

        both = list(top.intersection(bottom))
        top_only = list(top.difference(bottom))
        bottom_only = list(bottom.difference(top))

        plot_data['upper']['type'] = None
        plot_data['lower']['type'] = None

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(both), 'type'] = 'on_both'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(both), 'type'] = 'on_both'

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(top_only), 'type'] = 'top_only'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(top_only), 'type'] = 'top_only'

        plot_data['upper'].loc[plot_data['upper'][snp_col].isin(bottom_only), 'type'] = 'bottom_only'
        plot_data['lower'].loc[plot_data['lower'][snp_col].isin(bottom_only), 'type'] = 'bottom_only'

        custom_hue_colors = {
        "on_both"    : "#7B68EE",  
        "top_only"   : "#008080",
        "bottom_only": "#FF1493",
        }

        ax_upper = sns.scatterplot(
            data    =plot_data['upper'][~plot_data['upper']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_upper,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
            edgecolor='none'
        )

        ax_lower = sns.scatterplot(
            data    =plot_data['lower'][~plot_data['lower']['type'].isnull()].reset_index(drop=True),
            x       ='rel_pos', 
            y       ='log10p', 
            ax      =ax_lower,
            hue     ='type',
            palette =custom_hue_colors,
            size    =10,
            legend  =False,
            edgecolor='none'
        )

    # adjust layout
    plt.tight_layout()

    # render and draw the figure without showing it
    canvas = FigureCanvas(fig)
    fig.set_canvas(canvas)
    
    fig.canvas.draw()
    r = canvas.get_renderer()

    if top_annotations.shape[0]>0:
        
        if top_gen_col is not None:
            
            top_variants_toanno = plot_data['upper'][plot_data['upper'][snp_col].isin(top_annotations[snp_col])]\
                .reset_index(drop=True)
            top_variants_toanno = top_variants_toanno.merge(top_annotations, on=snp_col, how='left')
            top_variants_toanno = top_variants_toanno.rename(columns={top_gen_col:"GENENAME"})

        else:

            top_variants_toanno = plot_data['upper'][plot_data['upper'][snp_col].isin(top_annotations[snp_col])]\
                .reset_index(drop=True)
                
            # get gene names for upper plot
            if (top_variants_toanno.empty is not True):
                top_variants_toanno = annotate_snp(
                    top_variants_toanno,
                    chrom  =chr_col,
                    pos    =pos_col,
                    build  =build,
                    source =source,
                    gtf_path=gtf_path # type: ignore
                ).rename(columns={"GENE":"GENENAME"})

        # annotate upper plot
        ax_upper = manhattan_type_annotate(
            axes           =ax_upper, 
            data           =plot_data['upper'], 
            variants_toanno=top_variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )

    if bottom_annotations.shape[0]>0:

        if bottom_gen_col is not None:

            bottom_variants_toanno = plot_data['lower'][plot_data['lower'][snp_col].isin(bottom_annotations[snp_col])]\
                .reset_index(drop=True)
            bottom_variants_toanno = bottom_variants_toanno.merge(bottom_annotations, on=snp_col, how='left')
            bottom_variants_toanno = bottom_variants_toanno.rename(columns={bottom_gen_col:"GENENAME"})
        
        else:

            bottom_variants_toanno = plot_data['lower'][plot_data['lower'][snp_col].isin(bottom_annotations[snp_col])]\
            .reset_index(drop=True)

            # get gane names for bottom plot
            if (bottom_variants_toanno.empty is not True):
                bottom_variants_toanno = annotate_snp(
                    bottom_variants_toanno,
                    chrom  =chr_col,
                    pos    =pos_col,
                    build  =build,
                    source =source,
                    gtf_path=gtf_path # type: ignore
                ).rename(columns={"GENE":"GENENAME"})

        
        # annotate lower plot
        ax_lower= manhattan_type_annotate(
            axes           =ax_lower, 
            data           =plot_data['lower'], 
            variants_toanno=bottom_variants_toanno, 
            max_x_axis     =max_x_axis, 
            suggestive_line=suggestive_line, 
            genome_line    =genome_line
        )

    from matplotlib.lines import Line2D

    if len(top_highlights)>0 or len(bottom_highlights)>0:

        if len(both)>0:

            custom_dots = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#7B68EE", markersize=5, label='Top Hits both'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#008080", markersize=5, label=f'Top Hits {legend_top}'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF1493", markersize=5, label=f'Top Hits {legend_bottom}'),
            ]

        else:
            custom_dots = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#008080", markersize=5, label=f'Top Hits {legend_top}'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF1493", markersize=5, label=f'Top Hits {legend_bottom}'),
            ]

        # Add custom legend
        ax_upper.legend(handles=custom_dots, title='Legend', loc="best", fontsize=10)
    
    ax_lower.set_xlabel("Chromosome")
    ax_upper.set_xlabel("")
    
    # save ad show the plot
    plt.savefig(os.path.join(plots_dir, save_name), dpi=dpi)
    plt.show()

    return True

def brisbane_process_data(
    data_df: pd.DataFrame,
    chr_col: str = 'CHR',
    pos_col: str = 'POS',
    window_kb: float = 100
) -> dict:

    """Bin variants into fixed-size genomic windows and compute per-window SNP density.

    Parameters
    ----------
    data_df : pandas.DataFrame
        Input DataFrame containing genomic data.
    chr_col : str, optional
        Column name for chromosome identifiers. Defaults to 'CHR'.
    pos_col : str, optional
        Column name for base pair positions. Defaults to 'POS'.
    window_kb : float, optional
        Window size in kilobases used to bin variants along each chromosome. Defaults to 100.

    Returns
    -------
    dict
        A dictionary with the following keys:

        - data : pandas.DataFrame
            One row per non-empty window, with columns `chr_col`, 'window_start', 'count'
            (number of variants in the window) and 'rel_pos' (genome-wide relative position
            of the window's midpoint).
        - axis : pandas.DataFrame
            The center positions of each chromosome for x-axis tick placement.
        - max_count : int
            The maximum per-window SNP count.

    Raises
    ------
    TypeError
        If `data_df` is not a pandas DataFrame.
    ValueError
        If `chr_col` or `pos_col` are not found in the DataFrame, or if `window_kb` is not
        a positive number.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if not isinstance(window_kb, (int, float)) or window_kb <= 0:
        raise ValueError(f"window_kb must be a positive number, got {window_kb}.")

    window_bp = window_kb * 1000

    # Per-chromosome genomic offsets, in standard genomic order (reused from Manhattan plotting)
    chr_offsets = _compute_chr_offsets(data_df, chr_col=chr_col, pos_col=pos_col)

    # Vectorized windowed SNP counting: one groupby over the whole dataset, no per-chromosome loop
    window_start = ((data_df[pos_col] // window_bp) * window_bp).rename('window_start')
    density = pd.concat([data_df[chr_col], window_start], axis=1)
    density = density.groupby([chr_col, 'window_start']).size().rename('count').reset_index()

    density = density.merge(chr_offsets, on=chr_col, how='left')
    density['rel_pos'] = density['window_start'] + window_bp / 2 + density['cumulativechrlength']
    density = density.drop(columns=['cumulativechrlength'])

    density = density.sort_values(
        by=[chr_col, 'window_start'], key=lambda s: s.map(_chromosome_sort_key) if s.name == chr_col else s
    ).reset_index(drop=True)

    axis_center = find_chromosomes_center(density, chr_col=chr_col, chr_pos_col='rel_pos')

    return {
        'data': density,
        'axis': axis_center,
        'max_count': int(density['count'].max())
    }

def brisbane_draw(
    data_df: pd.DataFrame,
    chr_col: str,
    pos_col: str,
    plot_dir: str,
    window_kb: float = 100,
    save_name: str = 'brisbane_plot.png',
    chr_colors: list = ['#000000', '#E69F00'],
    point_size: float = 8,
    show_mean_line: bool = True,
    show_median_line: bool = True,
    mean_line_color: str = 'red',
    median_line_color: str = 'red',
    yaxis_margin: float = 2,
    ytick_step: Optional[float] = 1,
    figsize: tuple = (15, 6),
    dpi: int = 400
) -> bool:

    """Generate a Brisbane plot (genome-wide density of significant signals) for genomic data.

    This function bins variants into fixed-size windows along each chromosome and plots
    the number of variants per window against genomic position, following the style used
    in Yengo et al. (2022, Nature): windows are drawn as squares alternating in color by
    chromosome, with horizontal lines marking the genome-wide mean and/or median density.
    Unlike a Manhattan plot, the y-axis represents SNP count per window rather than
    -log10(p-value).

    Parameters
    ----------
    data_df : pandas.DataFrame
        DataFrame containing the SNP data to be plotted.
    chr_col : str
        Column name in data_df that contains the chromosome information.
    pos_col : str
        Column name in data_df that contains the position information.
    plot_dir : str
        Directory path where the plot will be saved.
    window_kb : float, optional
        Window size in kilobases used to bin variants along each chromosome. Defaults to 100.
    save_name : str, default='brisbane_plot.png'
        Filename for the saved plot.
    chr_colors : list, default=['#000000', '#E69F00']
        Two (or more) colors cycled across chromosomes to distinguish adjacent chromosomes.
    point_size : float, default=8
        Marker size for each window's point.
    show_mean_line : bool, default=True
        Whether to draw a horizontal line at the genome-wide mean window density.
    show_median_line : bool, default=True
        Whether to draw a horizontal line at the genome-wide median window density.
    mean_line_color : str, default='red'
        Color of the mean-density line.
    median_line_color : str, default='red'
        Color of the median-density line.
    yaxis_margin : float, default=2
        Additional margin added to the y-axis maximum.
    ytick_step : float, optional
        Spacing between y-axis ticks. Defaults to 1 (a tick at every integer density value).
        Pass None to fall back to matplotlib's automatic tick placement, or a larger value
        to space ticks further apart.
    figsize : tuple, default=(15, 6)
        Figure size (width, height) in inches. Reduce the height to bring y-axis ticks
        physically closer together without changing their values (a more compact plot).
    dpi : int, default=400
        Resolution of the saved image in dots per inch.

    Returns
    -------
    bool
        True if the plot was successfully created and saved.

    Raises
    ------
    TypeError
        If input data is not a pandas DataFrame.
    ValueError
        If required columns are not found in the input DataFrame.
    FileNotFoundError
        If the specified plot directory does not exist.

    Notes
    -----
    Saves a Brisbane plot image to the specified directory.
    """

    if not isinstance(data_df, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")
    if chr_col not in data_df.columns:
        raise ValueError(f"Column '{chr_col}' not found in the input DataFrame.")
    if pos_col not in data_df.columns:
        raise ValueError(f"Column '{pos_col}' not found in the input DataFrame.")
    if not os.path.exists(plot_dir):
        raise FileNotFoundError(f"Directory '{plot_dir}' not found.")

    # Suppress warnings about the number of chromosomes exceeding the number of palette colors
    warnings.filterwarnings("ignore", category=UserWarning)

    plot_data = brisbane_process_data(data_df, chr_col=chr_col, pos_col=pos_col, window_kb=window_kb)
    density = plot_data['data']

    max_x_axis = density['rel_pos'].max()

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax = sns.scatterplot(
        x        =density['rel_pos'],
        y        =density['count'],
        hue      =density[chr_col],
        palette  =chr_colors,
        marker   ='s',
        s        =point_size,
        ax       =ax,
        legend   =False,
        edgecolor='none'
    )

    ax.set_ylabel(f"Density of GWS SNPs within {window_kb:g} kb", fontsize=7)
    ax.set_xlabel("Chromosome", fontsize=7)

    ax.set_xlim(0, max_x_axis + 1000)
    ax.set_ylim(0, plot_data['max_count'] + yaxis_margin)

    if ytick_step is not None:
        ax.yaxis.set_major_locator(mticker.MultipleLocator(ytick_step))

    x_ticks = plot_data['axis']['center'].tolist()
    x_labels = plot_data['axis'][chr_col].astype(str).tolist()
    ax.set_xticks(ticks=x_ticks)
    ax.set_xticklabels(x_labels)
    ax.tick_params(axis='both', labelsize=7)

    if show_mean_line:
        mean_density = density['count'].mean()
        ax.axhline(
            mean_density,
            color    =mean_line_color,
            linestyle='solid',
            lw       =1,
            label    =f"Mean genome-wide density: {mean_density:.2g}"
        )

    if show_median_line:
        median_density = density['count'].median()
        ax.axhline(
            median_density,
            color    =median_line_color,
            linestyle='dashed',
            lw       =1,
            label    =f"Median genome-wide density: {median_density:.2g}"
        )

    if show_mean_line or show_median_line:
        ax.legend(loc='upper left', fontsize=8, frameon=False)

    plt.tight_layout()

    # save the plot
    plt.savefig(
        os.path.join(plot_dir, save_name), dpi=dpi
    )
    plt.show()

    return True