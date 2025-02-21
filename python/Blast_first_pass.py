# %%
import argparse
import cffi
import fnmatch
import glob
import os
import pathlib
import re
import subprocess
import sys
import time
from datetime import datetime
from enum import Enum
from io import StringIO
from pathlib import Path
from urllib.parse import urljoin
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from Bio import Entrez, SeqIO
from src.BlastProcessor import BLASTProcessor
from src.FileHandler import FileHandler
from src.SequenceProcessor import SequenceProcessor
from src.pytrsomix import SeqAnalyzer, TRScalculator
from src.stats import Stats

#Utility function to check remote file date
def get_remote_file_date(url):
    response = requests.head(url)
    if 'Last-Modified' in response.headers:
        return datetime.strptime(response.headers['Last-Modified'], '%a, %d %b %Y %H:%M:%S %Z')

#Utility function to check local file date:
def get_local_file_date(filepath):
    if os.path.exists(filepath):
        return datetime.fromtimestamp(os.path.getmtime(filepath))
    
#Utility function to download the file
def download_file(url,dest):
    print("Downloading taxonomy database this might take a while...")
    response = requests.get(url,stream=True)
    with open(dest,"wb") as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)
    print(f"Downloaded and replaced {dest}")

def check_and_update_taxonomy_file(taxonomy_db_path, force_download=False):
    if not os.path.exists(taxonomy_db_path):
        os.makedirs(taxonomy_db_path)
        print(f"Created directory for taxonomy database at {taxonomy_db_path}")
    
    filename = 'nucl_gb.accession2taxid.gz'
    local_filepath = os.path.join(taxonomy_db_path, filename)
    remote_url = urljoin("https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/", filename)
    local_date = get_local_file_date(local_filepath)
    remote_date = get_remote_file_date(remote_url)

    if local_date and remote_date:
        if local_date >= remote_date:
            print("Local taxonomy file is up to date.")
            return
        else:
            print("Local taxonomy file is outdated.")
            if not force_download:
                print("Run the program with --update_taxonomy to download the latest version.")
                sys.exit(1)

    if force_download or not local_date:
        download_file(remote_url,local_filepath)


#
def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description='''This program analyzes the blast output files to find potential full sequence candidates 
                                     that should be BLASTED once again and the results passed to Blast_sceond_pass.py''')
    parser.add_argument('--blast_output_path', help="Path to a folder containing blast results in a valid format", type=str, required=True)
    parser.add_argument('--taxonomy_db',help='Path to the taxonomy - accession database', required=True, type= str)
    parser.add_argument('--update_taxonomy', help="Download and overwrite outdated taxonomy file", required=False, action='store_true')
    parser.add_argument('--email', help= "Addres email to be used for connection with NCBI servers", required= True, type= str)
    add_ids = parser.add_mutually_exclusive_group()
    add_ids.add_argument('--ids_to_add_to_dictionary', help= "Comma separated list or single value of NCBI IDs", required= False, type=str)
    add_ids.add_argument('--ids_file', help="Path to the file containing NCBI IDs", required= False, type= str)
    exceptions = parser.add_mutually_exclusive_group()
    exceptions.add_argument('--taxids_to_add_to_exceptions', help="Comma separated list or single taxid to add to filtering exceptions", required= False,type=str)
    exceptions.add_argument('--taxids_file', help= "Path to the file containing the taxids that should be added to exceptions", required= False, type=str)
    args = parser.parse_args()

    blast_output = args.blast_output_path
    taxonomy_db = args.taxonomy_db
    email = args.email
    ids_to_add_to_dictionary = args.ids_to_add_to_dictionary
    ids_file = args.ids_file
    taxids_to_add_to_exceptions = args.taxids_to_add_to_exceptions
    taxids_file = args.taxids_file

    if not args.taxonomy_db:
        default_taxonomy_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxonomy")
        taxonomy_db_path = default_taxonomy_path
    else:
        taxonomy_db_path = args.taxonomy_db

    print("Checking taxonomy file....")
    check_and_update_taxonomy_file(taxonomy_db_path, args.update_taxonomy)

    #Process BLAST analysis results
    FileHandler.convert_to_txt(blast_output)
    FileHandler.filter_and_overwrite_blast_file(blast_output)
    blast_output_path = Path(blast_output)
    results_directory = blast_output_path.parent
    modified_blast_path = os.path.join(blast_output_path, "modified_blast")
    FileHandler.ensure_directory_exists(modified_blast_path)


    #Collect accessions from BLAST files
    accessions = BLASTProcessor.collect_accessions_from_blast_files(blast_output_path)
    
    # Filter taxonomy file based on collected accessions
    taxonomy_db = os.path.join(taxonomy_db,"nucl_gb.accession2taxid.gz")
    tax_df = SequenceProcessor.filter_taxonomy_file(taxonomy_db, accessions)
    # print("Collected Accessions:", accessions)
    # print("Filtered Taxonomy DataFrame:\n", tax_df) Debug

    # Create a mapping between accession numbers and taxids
    print(f"Creating a mapping between accessions and taxids....")
    taxid_accessions_dict = {}
    for index, row in tax_df.iterrows():
        accession = row[tax_df.columns[0]]
        taxid = row[tax_df.columns[1]]

        if taxid in taxid_accessions_dict:
            taxid_accessions_dict[taxid].append(accession)
        else:
            taxid_accessions_dict[taxid] = [accession]

    accession_to_taxid = {accession: taxid for taxid, accessions in taxid_accessions_dict.items() for accession in accessions}

    # Process BLAST files to map accessions to TaxIDs and add them as a new column
    print(f"Appending a new column with taxids to the blast files.....")
    BLASTProcessor.match_accessions_to_taxids_and_add_column(
        blast_output_path,
        modified_blast_path,
        lambda accession: BLASTProcessor.map_accession_to_taxid(accession, accession_to_taxid)
    )

    end_time = time.time()
    print(f"Processing completed in {end_time - start_time:.2f} seconds")
    
    #Search for the results CSV file in the results directory
    print(f"Searching for results file...")
    results_file = FileHandler.find_file_by_name('*_results.csv', folder= results_directory)
    if results_file:
        result_file = results_file[0]
        print(f"Results file found at : {results_file}") #MAYBE MAKE IT INTO ARGUMENT ?
    else:
        print("No '*_results.csv' file selected or found.")
    
    #Load the combined results into a DataFrame
    combined_results = pd.read_csv(result_file)

    #Set the Entrez e-mail for further processing
    Entrez.email = SequenceProcessor.validate_and_set_email(email)
    
    #Fetch taxonomic information for the collected genome IDs
    ncbi_ids = combined_results["GENOME"].unique().tolist()
    tax_map = SequenceProcessor.fetch_organism_taxids(ncbi_ids)

    #Filter and clean the taxonomic mapping
    filtered_organism_taxid_map = {SequenceProcessor.filter_key_parts(key): value for key, value in tax_map.items()}
    print(f"Species - taxid pairs detected in dataset : {filtered_organism_taxid_map}")

    #Append taxid information to the filtered map
    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)

    #Interact and update the taxid-accession map using the provided NCBI IDs
    BLASTProcessor.interact_and_update_dict(filtered_organism_taxid_map,ncbi_ids_list=ids_to_add_to_dictionary,file_path=ids_file)

    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
    species_info = {SequenceProcessor.filter_key_parts(key): value for key, value in species_info.items()}
    
    # Prepare NaN file path for storing the results with missing accessions
    nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")

    # Construct dictionary of all taxid - acessions pairs in our data
    results_dict = BLASTProcessor.construct_dict_from_files(modified_blast_path,nan_file)

    #Load exceptions to filtering based on taxids 
    exceptions = BLASTProcessor.ask_for_exception(exception_ids=taxids_to_add_to_exceptions,file_path=taxids_file)
    if exceptions:
        print(f"Exceptions to filtering added: {exceptions}")
    else:
        print(f"No exceptions provided!")
    
    #Convert all values in the results and species info to integers
    for key, value_set in results_dict.items():
        results_dict[key] = {int(val) for val in value_set}    
        
    for key, value_set in species_info.items():
        species_info[key] = {int(val) for val in value_set}  

    '''This filtering step does the following:
    1. Looks through all key-value set pairs in results_dict
    2. If one of the values in the set is in exception removes it from the set(but remembers it)
    3. If leftover values match atleast one present in species_info 
    AND it's the only one associated with a given sequence the record is kept
    4. That means that the records which had more than one value associated with it but the other ones were exceptions
    are preserved
    5.Unfortunetely this didnt help with the Klebsiella_pneumoniae_subsp._pneumoniae completely disappearing from the dataset
    6. Still this is not a bug but expected behaviour as KP sequences match to A LOT of taxids'''
    
    #Perform filtering based on the exceptions and matching taxid criteria
    print("Filtering keys....")
    filtered_keys = BLASTProcessor.filter_with_exceptions(results_dict,species_info,exceptions)
    filtered_keys_final = BLASTProcessor.unpack_single_element_sets(filtered_keys)

    #Process the filtered BLAST result files
    print("Filtering files...")
    FileHandler.process_files_with_filter(modified_blast_path,filtered_keys_final)
    
    #NEW

    #Find the filtered FASTA file that contains all sequences(but unique)
    filtered_fasta_file = FileHandler.find_file_by_name('unique_taxids_not_in_clusters_combined_sequences_unique_blastn_out.txt',folder= results_directory)
    filtered_fasta_file = pathlib.Path(filtered_fasta_file[0])

    #Ensure the output directory exists for storing the filtered results
    intermediate_results = os.path.join(results_directory,"intermediate_results")
    FileHandler.ensure_directory_exists(intermediate_results)

    BLASTProcessor.extract_full_TRS_sequences(combined_results,filtered_fasta_file,results_directory,state=1)

if __name__ == "__main__":
    main()



