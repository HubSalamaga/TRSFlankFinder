# TRSFlankFinder : Detection of unique sequences flanked by trinucleotide repeats

## Overview

> [!NOTE]
> 1. **`TRS_part.py`**: This script is used to obtain initial results for subsequent BLAST analysis.
>
> 2. **`TRS_BLAST.sh`** BLASTS sequences from a specific directory against nt database at 100% identity threshold.
>
> 3. **`Blast_part.py`**: This script is used to obtain final results from BLAST files. 
>
> 4. **Combined pipeline using `combined.py`**: This script executes all the steps above including the BLASTing step and automatically detects which version of the script to use depending on slurm availability (as of now threads and memory parameters can be changed only in the scripts themselves)

## Requirements:

1.Linux operating system preferably Ubuntu

2.Packages listed in the environment.yml file

## Installation:

> [!IMPORTANT]
> 1. **Environment Setup with Conda**: The necessary package list for script operation is in `environment.yml`.
>    
>    Quick installation command in terminal: `conda create env -f environment.yml -n TRS`
>    
> 2. Running `python/trsomix_wrapper/setup.py` to compile `_trs_omix`

## Run

> [!IMPORTANT]
> 1. **Activate Environment**: Use `conda activate TRS`. **All further operations should be performed in this environment**.
>
> 2. **Use -h with various scripts to see which arguments are optional and which are required for execution as well as descriptions of arugments**

```
usage: TRS_part.py [-h] --input_fasta_folder_path INPUT_FASTA_FOLDER_PATH
                   --tmin TMIN --tmax TMAX --mode MODE [--redo] [--cont]
                   --email EMAIL --threshold THRESHOLD --length_to_extract
                   LENGTH_TO_EXTRACT [--cd_hit_path CD_HIT_PATH]

This program extracts TRS sequences from a series of input genomes, allows for
length selection of extracted fragments, and prepares sequences for further
analysis.

optional arguments:
  -h, --help            show this help message and exit
  --input_fasta_folder_path INPUT_FASTA_FOLDER_PATH
                        Path to a folder containing genomes in fasta format
                        from which TRS sequences will be extracted
  --tmin TMIN           Minimum length of TRS sequences
  --tmax TMAX           Maximum length of TRS sequences
  --mode MODE           Mode of operation, must be 0 or 1
  --redo                Redo the analysis if results directory already exists
  --cont                Continue the analysis from saved TRS results file
  --email EMAIL         Address e-mail to be used for connection with NCBI
                        databases
  --threshold THRESHOLD
                        Identity threshold for clustering using cdhit has to
                        be between 0.8 and 1.0
  --length_to_extract LENGTH_TO_EXTRACT
                        Length of flanking sequences to be extracted from the
                        full TRS sequence
  --cd_hit_path CD_HIT_PATH
                        Path to the cd-hit-est executable

```

```
usage: Blast_part.py [-h] --blast_output_path BLAST_OUTPUT_PATH --taxonomy_db
                     TAXONOMY_DB --email EMAIL
                     [--ids_to_add_to_dictionary IDS_TO_ADD_TO_DICTIONARY | --ids_file IDS_FILE]
                     [--taxids_to_add_to_exceptions TAXIDS_TO_ADD_TO_EXCEPTIONS | --taxids_file TAXIDS_FILE]

This program analyzes the blast output files to find species specific
sequences coming from genomes analyzed in the previous step

optional arguments:
  -h, --help            show this help message and exit
  --blast_output_path BLAST_OUTPUT_PATH
                        Path to a folder containing blast results in a valid
                        format
  --taxonomy_db TAXONOMY_DB
                        Path to the taxonomy - accession database
  --email EMAIL         Addres email to be used for connection with NCBI
                        servers
  --ids_to_add_to_dictionary IDS_TO_ADD_TO_DICTIONARY
                        Comma separated list or single value of NCBI IDs
  --ids_file IDS_FILE   Path to the file containing NCBI IDs
  --taxids_to_add_to_exceptions TAXIDS_TO_ADD_TO_EXCEPTIONS
                        Comma separated list or single taxid to add to
                        filtering exceptions
  --taxids_file TAXIDS_FILE
                        Path to the file containing the taxids that should be
                        added to exceptions
```

```
usage: TRS_full_pipeline.py [-h] --input_fasta_folder_path
                            INPUT_FASTA_FOLDER_PATH --tmin TMIN --tmax TMAX
                            --mode MODE [--redo] [--cont] --email EMAIL
                            --threshold THRESHOLD --length_to_extract
                            LENGTH_TO_EXTRACT [--cd_hit_path CD_HIT_PATH]
                            [--taxonomy_db TAXONOMY_DB] [--update_taxonomy]
                            [--ids_to_add_to_dictionary IDS_TO_ADD_TO_DICTIONARY | --ids_file IDS_FILE]
                            [--taxids_to_add_to_exceptions TAXIDS_TO_ADD_TO_EXCEPTIONS | --taxids_file TAXIDS_FILE]
                            --blast_db BLAST_DB

This program extracts TRS sequences from a series of input genomes, allows for
length selection of extracted fragments, and prepares sequences for further
analysis.

optional arguments:
  -h, --help            show this help message and exit
  --input_fasta_folder_path INPUT_FASTA_FOLDER_PATH
                        Path to a folder containing genomes in fasta format
                        from which TRS sequences will be extracted
  --tmin TMIN           Minimum length of TRS sequences
  --tmax TMAX           Maximum length of TRS sequences
  --mode MODE           Mode of operation, must be 0 or 1
  --redo                Redo the analysis if results directory already exists
  --cont                Continue the analysis from saved TRS results file (currently disabled)
  --email EMAIL         Address e-mail to be used for connection with NCBI
                        databases
  --threshold THRESHOLD
                        Identity threshold for clustering using cdhit has to
                        be between 0.8 and 1.0
  --length_to_extract LENGTH_TO_EXTRACT
                        Length of flanking sequences to be extracted from the
                        full TRS sequence
  --cd_hit_path CD_HIT_PATH
                        Path to the cd-hit-est executable
  --taxonomy_db TAXONOMY_DB
                        Path to the taxonomy - accession database
  --update_taxonomy     Download and overwrite outdated taxonomy file
  --ids_to_add_to_dictionary IDS_TO_ADD_TO_DICTIONARY
                        Comma separated list or single value of NCBI IDs
  --ids_file IDS_FILE   Path to the file containing NCBI IDs
  --taxids_to_add_to_exceptions TAXIDS_TO_ADD_TO_EXCEPTIONS
                        Comma separated list or single taxid to add to
                        filtering exceptions
  --taxids_file TAXIDS_FILE
                        Path to the file containing the taxids that should be
                        added to exceptions
  --blast_db BLAST_DB   Path to your blast_db

```
# Acknowledgements

The content of the repository is derived from its predecessor, available at:
[https://github.com/TRS-omix/software](https://github.com/TRS-omix/software) described in more detail in **Sebastian Sakowski, Marta Majchrzak, Jacek Waldmajer, Pawel Parniewski**: *TRS-omix: a new search engine for trinucleotide flanked sequences*. 2021.

