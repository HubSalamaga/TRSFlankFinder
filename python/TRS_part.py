import cffi
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
import time
import argparse
import sys
from src.SequenceProcessor import SequenceProcessor
from src.FileHandler import FileHandler
from src.pytrsomix import SeqAnalyzer,TRScalculator
from src.BlastProcessor import BLASTProcessor
from src.stats import Stats
import subprocess
import pathlib
import requests
from datetime import datetime
from urllib.parse import urljoin
#comment

import pathlib
def check_mode(value):
    ivalue = int(value)
    if ivalue < 0 or ivalue > 1:
        raise argparse.ArgumentTypeError(f"Mode must be 0 or 1, got {value}")
    return ivalue

def check_threshold(value):
    fvalue = float(value)
    if fvalue < 0.8 or fvalue > 1.0:
        raise argparse.ArgumentTypeError(f"Identity threshold for cdhit must be between 0.8 and 1.0, got {value}")
    return fvalue

def find_latest_results_directory(base_dir, base_name):
    """
    Find the latest results directory that matches the base name pattern.
    """
    matching_dirs = [d for d in os.listdir(base_dir) if d.startswith(base_name)]
    if not matching_dirs:
        return None
    latest_dir = max(matching_dirs, key=lambda d: os.path.getmtime(os.path.join(base_dir, d)))
    return os.path.join(base_dir, latest_dir)

def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='''This program extracts TRS sequences from a series of input genomes,
                                     allows for length selection of extracted fragments, and prepares sequences for further analysis.''')
    parser.add_argument('--input_fasta_folder_path', help="Path to a folder containing genomes in fasta format from which TRS sequences will be extracted", 
                        required=True)
    parser.add_argument('--tmin', help="Minimum length of TRS sequences", required=True, type=int)
    parser.add_argument('--tmax', help="Maximum length of TRS sequences", required=True, type=int)
    parser.add_argument('--mode', help="Mode of operation, must be 0 or 1", required=True, type=check_mode)
    parser.add_argument('--redo', help="Redo the analysis if results directory already exists", required=False, action='store_true')
    parser.add_argument('--cont', help="Continue the analysis from saved TRS results file", required=False, action='store_true')
    parser.add_argument('--email', help="Address e-mail to be used for connection with NCBI databases", required=True, type=str)
    parser.add_argument('--threshold', help="Identity threshold for clustering using cdhit has to be between 0.8 and 1.0", required=True, type=check_threshold)
    parser.add_argument('--length_to_extract', help="Length of flanking sequences to be extracted from the full TRS sequence", required=True, type=int)
    parser.add_argument('--cd_hit_path', help="Path to the cd-hit-est executable", required=False, type=str)
    args = parser.parse_args()


    #Validate and set email used for querying ncbi services
    Entrez.email = SequenceProcessor.validate_and_set_email(args.email)
    print(f"Email adress is currently set to {Entrez.email}")

    # Define the base results directory name and ensure it exists
    input_fasta_folder_path_name = os.path.basename(args.input_fasta_folder_path)
    base_results_directory = os.path.join(os.getcwd(), f"{input_fasta_folder_path_name}_results")
    results_directory = base_results_directory
    FileHandler.ensure_directory_exists(base_results_directory)

    # Define the name of the CSV file that will store the results of the analysis 
    name_of_csv_file_storing_TRS_analysis_results = input_fasta_folder_path_name + "_results.csv"
    path_of_folder_storing_TRS_analysis_results = os.path.join(base_results_directory, "TRS_output")
    FileHandler.ensure_directory_exists(path_of_folder_storing_TRS_analysis_results)
    path_of_csv_file_storing_TRS_analysis_results = os.path.join(path_of_folder_storing_TRS_analysis_results, 
                                                                 name_of_csv_file_storing_TRS_analysis_results)

    if args.cont:
        # Check for the initial results file
        if os.path.exists(path_of_csv_file_storing_TRS_analysis_results):
            combined_trs_results = pd.read_csv(path_of_csv_file_storing_TRS_analysis_results) # no idea how to do that correctly
            print("DataFrame loaded from initial results. Resuming analysis...")
        else:
            # Check for the latest renamed results directory
            latest_results_directory = find_latest_results_directory(os.getcwd(), input_fasta_folder_path_name)
            if latest_results_directory:
                path_of_csv_file_storing_TRS_analysis_results = os.path.join(latest_results_directory, "TRS_output", name_of_csv_file_storing_TRS_analysis_results)
                if os.path.exists(path_of_csv_file_storing_TRS_analysis_results):
                    combined_trs_results = pd.read_csv(path_of_csv_file_storing_TRS_analysis_results)
                    print(f"DataFrame loaded from {latest_results_directory}. Resuming analysis...")
                    results_directory = latest_results_directory
                else:
                    print("No existing results file found to continue from. Exiting program.")
                    sys.exit(1)
            else:
                print("No existing results directory found to continue from. Exiting program.")
                sys.exit(1)
    else:
        # Check if the provided input directory exists
        if not os.path.exists(args.input_fasta_folder_path):
            print("Provided directory path does not exist. Quitting...")
            sys.exit(1)

        print("Existing directory provided. Proceeding....")
        # Check for fasta files in the directory
        fasta_files = [f for f in os.listdir(args.input_fasta_folder_path) if f.endswith('.fasta') or f.endswith('.fa')]
        if not fasta_files:
            print("No fasta files found in the provided directory. Quitting...")
            sys.exit(1)

        # Check if results already exist and handle redo flag
        if os.path.exists(path_of_csv_file_storing_TRS_analysis_results) and not args.redo:
            print("Results already exist. Use --redo to overwrite or --cont to continue from existing results. Exiting program.")
            sys.exit(1)

        print("Starting analysis...")

        trs_calculators = []

        # Iterate over each fasta file and calculate TRS
        for fasta_file in fasta_files:
            path_to_input_fasta = os.path.join(args.input_fasta_folder_path, fasta_file)
            if not os.path.exists(path_to_input_fasta):
                print(f"File '{fasta_file}' does not exist! Skipping....")
                continue

            # Define the TRS file path
            trs_file = os.path.abspath(os.path.join(args.input_fasta_folder_path, 'trs.txt'))
            print(f"TRS file found at {trs_file}")
            if not os.path.exists(trs_file):
                script_dir = os.path.dirname(os.path.abspath(__file__))
                trs_file = os.path.abspath(os.path.join(script_dir, 'trs.txt'))
                if not os.path.isfile(trs_file) or not os.access(trs_file, os.R_OK):
                    print(f"Cannot read TRS file at: {trs_file}. Check permissions or file encoding.")
                    sys.exit(1)
            
            interiors_file = os.path.join(path_of_folder_storing_TRS_analysis_results, 'interiors.txt')
            print(f"Interiors file will be created at: {interiors_file}")

            try:
                # Initialize TRS calculator for each sequence and perform TRS search
                trs_calculator = TRScalculator(sequence=path_to_input_fasta, trs=trs_file.encode('utf-8') if isinstance(trs_file, str) else trs_file, 
                                               interiors=interiors_file.encode('utf-8') if isinstance(interiors_file, str) else interiors_file,
                                               tmin=args.tmin, tmax=args.tmax, mode=args.mode)
                trs_calculator.calculate()
                trs_calculators.append(trs_calculator)
            except Exception as e:
                print(f"An error occurred while processing '{fasta_file}': {e}")
                continue

        list_of_trs_results = []

        # Iterate over each TRScalculator instance
        for trs_calculator in trs_calculators:
            # Extract results from the calculator
            result = trs_calculator.Result
            # Append the result to the list
            list_of_trs_results.append(result)

        # Concatenate all results into a single DataFrame
        combined_trs_results = pd.concat(list_of_trs_results, ignore_index=True)

        # Remove ">" from >SEQ column
        combined_trs_results[">SEQ"] = combined_trs_results[">SEQ"].str.replace(">","")

        # Save the results of the first analysis step to the CSV file
        combined_trs_results.to_csv(path_of_csv_file_storing_TRS_analysis_results, index=False)
        print(f"Results saved to {path_of_csv_file_storing_TRS_analysis_results}")

    trs_time = time.time() #Record time taken for TRS search
    print(f"TRS took {trs_time - start_time} seconds")

    l_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    r_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    l_chars = int(l_chars)
    r_chars = int(r_chars)
    print(combined_trs_results)
    combined_trs_results = SequenceProcessor.extract_sequences(combined_trs_results, l_chars, r_chars)

    results_directory_after_flanks_extracted = f"{results_directory}_L{l_chars}_R{r_chars}"
    results_directory_after_flanks_extracted_path = os.path.join(os.path.dirname(results_directory), results_directory_after_flanks_extracted)

    print(f"Results directory is set to: {results_directory_after_flanks_extracted_path}")

    # Rename the results directory if necessary
    if not args.cont:
        if not os.path.exists(results_directory_after_flanks_extracted_path):
            # Rename the existing results directory
            os.rename(results_directory, results_directory_after_flanks_extracted_path)
            print(f"The results directory has been renamed to: {results_directory_after_flanks_extracted_path}")
        else:
            print(f"Directory {results_directory_after_flanks_extracted_path} already exists. Consider using a different name or removing the existing directory.")

    results_directory = results_directory_after_flanks_extracted_path

    ncbi_ids = combined_trs_results["GENOME"].unique().tolist()
    if Entrez.email and Entrez.email == args.email :
         print(f"Email adress is still set to {Entrez.email}")
    organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids,email = Entrez.email)

    # Map NCBI IDs to taxonomic names
    combined_trs_results['Taxonomic Name'] = None
    combined_trs_results['Taxonomic Name'] = combined_trs_results['GENOME'].map(organism_map)
    
    # Identify unmatched genomes
    unmatched_genomes = combined_trs_results[combined_trs_results['Taxonomic Name'].isnull()]["GENOME"].unique()

    if len(unmatched_genomes) > 0:
        print(f"Warning: Some genome IDs could not be matched with taxonomic names: {unmatched_genomes}")
    
    # Generate sequence IDs for left and right flanking sequences 
    combined_trs_results['L_id'] = combined_trs_results['Taxonomic Name'] + '_L' + combined_trs_results['L-No'].astype(str)
    combined_trs_results['R_id'] = combined_trs_results['Taxonomic Name'] + '_R' + combined_trs_results['R-No'].astype(str)
    
    #Prepare a dataframe containing flanking sequneces and their IDs.
    sequences_df = combined_trs_results[['SEQ_L', 'SEQ_R', 'L_id', 'R_id']]

    #Write the extracted sequences to a FASTA file
    path_of_folder_storing_TRS_analysis_results = os.path.join(results_directory, "TRS_output")
    fasta_files_with_flanks = os.path.join(path_of_folder_storing_TRS_analysis_results, "combined_sequences.fasta")
    with open(fasta_files_with_flanks, 'w') as fasta_file:
        for _, row in sequences_df.iterrows():
            # Write left sequence
            fasta_file.write(f'>{row["L_id"]}\n')
            fasta_file.write(f'{row["SEQ_L"]}\n')
            # Write right sequence
            fasta_file.write(f'>{row["R_id"]}\n')
            fasta_file.write(f'{row["SEQ_R"]}\n')
    
    #Create unique FASTA file with renamed sequences (including ids and L/R identifiers)
    fasta_files_with_flanks_unique = os.path.join(path_of_folder_storing_TRS_analysis_results, "combined_sequences_unique.fasta")
    SequenceProcessor.rename_sequences(fasta_files_with_flanks, fasta_files_with_flanks_unique)
    
    #Ensure that directory for CD-HIT results exists
    cd_hit_results_folder = os.path.join(results_directory, "cd-hit-results")
    FileHandler.ensure_directory_exists(cd_hit_results_folder)
    #Define the path for output files
    cd_hit_output_file = os.path.join(cd_hit_results_folder, "combined_sequences_unique_cdhit")
    cd_hit_path = args.cd_hit_path

    #Run CD-HIT clustering on the sequences to identify the similar ones
    results_directory = SequenceProcessor.run_cdhit(cd_hit_path, input_file=fasta_files_with_flanks_unique, output_file=cd_hit_output_file,
                                                    results_directory=results_directory, sc=1, c=args.threshold)
    #Prepare for further processing
    cd_hit_results_folder = os.path.join(results_directory, "cd-hit-results")
    cdhit_clusters_file = os.path.join(cd_hit_results_folder, "combined_sequences_unique_cdhit.clstr")
    
    #NEW WAY
    input_fasta = os.path.join(results_directory,"TRS_output","combined_sequences_unique.fasta")
    output_folder_cluster = os.path.join(cd_hit_results_folder,"fasta_clusters")
    FileHandler.create_fasta_for_all_clusters(cdhit_clusters_file, input_fasta, output_folder_cluster,create_individual_files=True)
    
    cluster_size_1_fasta = os.path.join(output_folder_cluster, "cluster_size_1.fasta")
    sequences_after_clusters_filtering_folder = os.path.join(results_directory, "filtered_sequences")
    FileHandler.ensure_directory_exists(sequences_after_clusters_filtering_folder)
    renamed_fasta_file = os.path.join(sequences_after_clusters_filtering_folder, "not_in_clusters_combined_sequences_unique.fasta")
    
    try:
        with open(cluster_size_1_fasta, 'rb') as src_file, open(renamed_fasta_file, 'wb') as dest_file:
            dest_file.write(src_file.read())
        print(f"Copied and renamed {cluster_size_1_fasta} to {renamed_fasta_file}")
    except FileNotFoundError:
        print(f"Error: {cluster_size_1_fasta} not found.")
    except PermissionError:
        print(f"Error: Permission denied when accessing {cluster_size_1_fasta} or {renamed_fasta_file}.")
    except Exception as e:
        print(f"An error occurred: {e}")

    fasta_files_with_flanks_unique = os.path.join(results_directory,"TRS_output","combined_sequences_unique.fasta")
    blast_results_folder = os.path.join(results_directory, "blast_output")
    FileHandler.ensure_directory_exists(blast_results_folder)


    #NEW WAY
    
if __name__ == "__main__":
    main()
