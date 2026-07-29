import logging
import requests

import pandas as pd

from pathlib import Path

from ..core.utils import download_file, unzip_file_flat, extract_gz_file

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def get_trumpet_quantitative_example() -> Path:
    """
    Downloads and extracts a quantitative GWAS summary statistics example file for generating a trumpet plot.
    This function checks if the BBJ Height GWAS summary statistics file already exists in the local data directory.
    If not, it downloads the file from the Human Database (DBCLS), extracts it from the zip archive,
    decompresses the gzipped file, and saves it to the data/sumstats directory.
    
    Returns
    -------
        Path: Path to the extracted and decompressed GWAS summary statistics file
               (2019_BBJ_Height_autosomes_BOLT.txt)

    Side effects
    ------------
        - Creates directories if they don't exist
        - Downloads files from the internet if necessary
        - Extracts and decompresses files
    """

    library_path = Path(__file__).resolve().parent.parent

    URL = "https://humandbs.dbcls.jp/files/hum0014/hum0014.v15.ht.v1.zip"
    FILENAME = "hum0014.v15.ht.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    example_path = LOCAL_PATH / "2019_BBJ_Height_autosomes_BOLT.txt"

    if example_path.exists():
        logger.info(f"File already exists: {example_path}")
        return example_path
    
    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz = unzip_file_flat(local_filename, "hum0014.v15.ht.v1/2019_BBJ_Height_autosomes_BOLT.txt.gz", LOCAL_PATH, remove_zip=True)

    logger.info(f"Decompressing file: {extracted_gz}")
    uncompressed_file = extract_gz_file(extracted_gz, LOCAL_PATH, remove_gz=True)

    return uncompressed_file


def get_top_loci_trumpet_quantitative() -> Path:
    """
    Downloads and processes top genetic loci data from a GWAS study on height.

    This function retrieves a dataset from the supplementary materials of a 
    Biobank Japan (BBJ) height GWAS study published in Nature Communications 
    (2019). It downloads an Excel file containing significant loci, extracts 
    relevant columns, filters out X chromosome variants, and saves the 
    processed data as a CSV file.

    Returns
    -------
        Path: Path to the processed CSV file containing filtered top genetic loci.
        
    Notes
    -----
        - The function checks if the processed CSV already exists before downloading
        - The original Excel file is deleted after processing
        - Only autosomal variants are kept (X chromosome variants are filtered out)
        - The processed file contains 'Variants' and 'New_Locus' columns
    """

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12276-5/MediaObjects/41467_2019_12276_MOESM5_ESM.xlsx"
    output_filename = "2019_BBJ_Height_autosomes_BOLT_loci.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_top = pd.read_excel(output_path, engine='openpyxl')

    df_top = df_top[['Unnamed: 7', 'Unnamed: 9']].copy()
    df_top = df_top.rename(columns={'Unnamed: 7': 'Variants', 'Unnamed: 9': 'New_Locus'}).dropna().iloc[1:, :].copy()

    mask_autosome = df_top['Variants'].str.contains('X')

    df_top = df_top[~mask_autosome].reset_index(drop=True)

    df_top.to_csv(output_path.with_suffix('.csv'), index=False, sep='\t')
    logger.info(f"Saved top hits to: {output_path.with_suffix('.csv')}")

    output_path.unlink()

    return output_csv


def get_top_cond_trumpet_quantitative() -> Path:
    """
    Downloads and processes conditional analysis data for height-associated variants from BBJ study.
    
    This function retrieves conditional analysis data for height-associated variants from a BBJ 
    (BioBank Japan) 2019 study. It downloads an Excel file from the provided URL, extracts 
    relevant columns (rsID, candidate genes, chromosome, and position), filters for autosomal 
    variants, and saves the processed data as a tab-separated CSV file.
    
    The function checks if the processed CSV file already exists before downloading to avoid 
    redundant operations. After processing, the original Excel file is deleted to save space.
    
    Returns
    -------
        Path: Path to the processed CSV file containing the top conditional analysis hits
              for height-associated variants.
    
    Side effects
    ------------
        - Downloads an Excel file from an external URL if not already processed
        - Creates a CSV file in the data/sumstats directory
        - Removes the downloaded Excel file after processing
    """

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12276-5/MediaObjects/41467_2019_12276_MOESM6_ESM.xlsx"
    output_filename = "2019_BBJ_Height_autosomes_BOLT_cond.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_top = pd.read_excel(output_path, engine='openpyxl', header=1)

    df_top = df_top[['rsID', 'Candiate gene(s)', 'CHRa', 'POSa']].copy()
    df_top = df_top.rename(columns={'Candiate gene(s)': 'Gene', 'CHRa': 'CHR', 'POSa': 'POS'}).dropna()

    mask_autosome = (df_top['CHR']!='X')

    df_top = df_top[mask_autosome].reset_index(drop=True)

    df_top.to_csv(output_path.with_suffix('.csv'), index=False, sep='\t')
    logger.info(f"Saved top hits to: {output_path.with_suffix('.csv')}")

    output_path.unlink()

    return output_csv


def get_yengo_height_independent_signals() -> Path:
    """
    Downloads and processes the independent COJO signals from the GIANT height GWAS.

    This function retrieves Supplementary Table 5 from Yengo et al. (2022), which lists the
    12,111 conditionally independent (COJO) genome-wide-significant SNPs identified in the
    GIANT/UK Biobank meta-analysis of human height, along with a precomputed signal-density
    metric (number of other COJO SNPs within 100 kb) used for the paper's Brisbane plot
    (Fig. 2).

    Yengo, L., Vedantam, S., Marouli, E. et al. A saturated map of common genetic variants
    associated with human height. *Nature* **610**, 704-712 (2022).
    https://doi.org/10.1038/s41586-022-05275-y

    Returns
    -------
        Path: Path to the processed CSV file containing the independent signals
              (columns: CHR, POS, SIGNAL_DENSITY, AFR, EAS, HIS). Coordinates are on
              genome build hg19/GRCh37.

    Notes
    -----
        - The function checks if the processed CSV already exists before downloading
        - The original Excel file is deleted after processing
        - Only autosomal variants are present in the source table (no X/Y/MT signals)
    """

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05275-y/MediaObjects/41586_2022_5275_MOESM6_ESM.xlsx"
    output_filename = "2022_Yengo_Height_COJO_independent_signals.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    output_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_signals = pd.read_excel(output_path, engine='openpyxl')

    df_signals = df_signals.rename(columns={
        'Chromosome': 'CHR',
        'Position (hg37)': 'POS',
        'Signal Density (number of other COJO SNPs within 100 kb)': 'SIGNAL_DENSITY',
        'Is this signal contributed by African ancesries (1=Yes; 0 = No)': 'AFR',
        'Is this signal contributed by East-Asian ancesries (1=Yes; 0 = No)': 'EAS',
        'Is this signal contributed by Hispanic Ethnicities (1=Yes; 0 = No)': 'HIS',
    })

    df_signals.to_csv(output_csv, index=False, sep='\t')
    logger.info(f"Saved independent signals to: {output_csv}")

    output_path.unlink()

    return output_csv


def get_trumpet_binary_example() -> Path:
    """
    Downloads and prepares a GWAS summary statistics file for analysis.
    This function retrieves a specific GWAS summary statistics file (GCST90278092)
    from the EBI FTP server. It checks if the file already exists locally before
    downloading. After download, it decompresses the gzipped file.

    Returns
    -------
        Path: Path to the uncompressed GWAS summary statistics file.
    
    Notes
    -----
        - The file is stored in the 'data/sumstats' directory relative to the package root
        - If the file already exists locally, no download is performed
        - The original gzipped file is removed after decompression
    """

    library_path = Path(__file__).resolve().parent.parent

    URL = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90278001-GCST90279000/GCST90278092/harmonised/GCST90278092.h.tsv.gz"
    FILENAME = "GCST90278092.h.tsv.gz"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    example_path = LOCAL_PATH / "GCST90278092.h.tsv"

    if example_path.exists():
        logger.info(f"File already exists: {example_path}")
        return example_path
    
    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

#
    logger.info(f"Decompressing file: {local_filename}")
    uncompressed_file = extract_gz_file(local_filename, LOCAL_PATH, remove_gz=True)

    return uncompressed_file


def get_height_japanese_gwas() -> tuple:
    """
    Downloads and prepares Japanese GWAS (Genome-Wide Association Studies) data for height.

    This function retrieves data from the Japanese Biobank (BBJ) for both female and male subjects.
    If the data files already exist locally, it returns their paths without downloading.
    Otherwise, it downloads the zip archive, extracts the relevant files, and decompresses them.

    Returns
    -------
        tuple: A tuple containing two Path objects:
            - Path to the female height GWAS data file
            - Path to the male height GWAS data file

    Notes
    -----
        - The data is downloaded from the Japanese Database of Clinical and Genomic Data, specifically the hum0014 dataset version 6 with 158k samples.
        - Akiyama, Masato, et al. Characterizing rare and low-frequency height-asssociated variants in the Japanese population. *Nature Communications* **10**, (2019): 4393. https://doi.org/10.1038/s41467-019-12276-5
    
    """

    library_path = Path(__file__).resolve().parent.parent

    URL = r"https://humandbs.dbcls.jp/files/hum0014/hum0014.v6.158k.v1.zip"
    FILENAME = "hum0014.v6.158k.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    female_path = LOCAL_PATH / "Female_2017_BMI_BBJ_autosome.txt"
    male_path = LOCAL_PATH / "Male_2017_BMI_BBJ_autosome.txt"

    if female_path.exists() and male_path.exists():
        logger.info(f"Files {female_path} and {male_path} already exist")
        return female_path, male_path

    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_f = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Female_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=False)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_m = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Male_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=True)

    logger.info(f"Decompressing file: {extracted_gz_f}")
    uncompressed_file_f = extract_gz_file(extracted_gz_f, LOCAL_PATH, remove_gz=True)

    logger.info(f"Decompressing file: {extracted_gz_m}")
    uncompressed_file_m = extract_gz_file(extracted_gz_m, LOCAL_PATH, remove_gz=True)

    return uncompressed_file_f, uncompressed_file_m


def get_bmi_japanese_gwas_sex_dependent() -> Path:
    """
    Download and process sex-dependent BMI GWAS data from the Japanese Biobank (BBJ).
    
    This function:
    1. Downloads a supplementary Excel file from a 2017 Nature Genetics paper
    2. Extracts relevant SNP information including sex-specific associations
    3. Filters for SNPs with rs IDs
    4. Converts the data to CSV format with columns: SNP, Sex, CHR, POS
    5. Removes the original Excel file after processing
    
    The study is the BioBank Japan GWAS which analyzed sex-dependent genetic 
    effects on BMI.
    
    Returns
    -------
        Path: Path to the processed CSV file containing the sex-dependent BMI 
              GWAS results
    
    Note
    ----
        If the CSV file already exists, it will not be re-downloaded or processed.
    """

    library_path = Path(__file__).resolve().parent.parent

    url = r"https://static-content.springer.com/esm/art%3A10.1038%2Fng.3951/MediaObjects/41588_2017_BFng3951_MOESM6_ESM.xlsx"
    output_filename = "2017_BBJ_bmi_supplementary.xlsx"

    output_path = library_path / "data" / "sumstats" / output_filename

    output_csv = output_path.with_suffix('.csv')
    if output_csv.exists():
        logger.info(f"File already exists: {output_csv}")
        return output_csv

    response = requests.get(url, stream=True)  # Stream to handle large files
    if response.status_code == 200:
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
                file.write(chunk)
        logger.info(f"Downloaded file: {output_path}")
    else:
        logger.info(f"Failed to download file. Status code: {response.status_code}")

    df_top = pd.read_excel(output_path, engine='openpyxl', header=[0,1,2,3], sheet_name='S.Table4')

    df_top = df_top[df_top.columns[0:4]].copy()
    df_top.columns = df_top.columns.droplevel([0,2,3])
    df_top.columns = ['SNP', 'Sex', 'CHR', 'POS']

    mask_rsid = df_top['SNP'].str.startswith('rs')

    df_top = df_top[mask_rsid].reset_index(drop=True)

    df_top.to_csv(output_path.with_suffix('.csv'), index=False, sep='\t')
    logger.info(f"Saved top hits to: {output_path.with_suffix('.csv')}")

    output_path.unlink()

    return output_csv


def get_bmi_japanese_gwas() -> tuple:

    library_path = Path(__file__).resolve().parent.parent

    URL = r"https://humandbs.dbcls.jp/files/hum0014/hum0014.v6.158k.v1.zip"
    FILENAME = "hum0014.v6.158k.v1.zip"
    LOCAL_PATH = library_path / "data" / "sumstats"

    LOCAL_PATH.mkdir(parents=True, exist_ok=True)

    local_filename = LOCAL_PATH / FILENAME

    logger.info(f"Downloading file: {URL} to {local_filename}")
    download_file(URL, local_filename=local_filename)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_f = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Female_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=False)

    logger.info(f"Extracting file: {local_filename}")
    extracted_gz_m = unzip_file_flat(local_filename, "hum0014.v6.158k.v1/Male_2017_BMI_BBJ_autosome.txt.gz", LOCAL_PATH, remove_zip=True)

    logger.info(f"Decompressing file: {extracted_gz_f}")
    uncompressed_file_f = extract_gz_file(extracted_gz_f, LOCAL_PATH, remove_gz=True)

    logger.info(f"Decompressing file: {extracted_gz_m}")
    uncompressed_file_m = extract_gz_file(extracted_gz_m, LOCAL_PATH, remove_gz=True)

    return uncompressed_file_f, uncompressed_file_m