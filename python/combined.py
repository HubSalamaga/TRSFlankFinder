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

#Check and update taxonomy_file
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

# Utility function to check if the mode argument was provided and is within 0 to 1 range 
def check_mode(value):
    ivalue = int(value)
    if ivalue < 0 or ivalue > 1:
        raise argparse.ArgumentTypeError(f"Mode must be 0 or 1, got {value}")
    return ivalue
# Utility function to check if the threshold argument is within allowed range
def check_threshold(value):
    fvalue = float(value)
    if fvalue < 0.8 or fvalue > 1.0:
        raise argparse.ArgumentTypeError(f"Identity threshold for cdhit must be between 0.8 and 1.0, got {value}")
    return fvalue
#Function to find the latest results directory based on modification time - currently not used for anything
def find_latest_results_directory(base_dir, base_name):
    """
    Find the latest results directory that matches the base name pattern.
    """
    #Get all directories in base_dir that start with base_name
    matching_dirs = [d for d in os.listdir(base_dir) if d.startswith(base_name)]
    if not matching_dirs:
        return None
    #Find the latest directory based on modification time
    latest_dir = max(matching_dirs, key=lambda d: os.path.getmtime(os.path.join(base_dir, d)))
    return os.path.join(base_dir, latest_dir)

#Main function of the script, parses arguments and runs the TRS extraction pipeline 
def main():
    start_time = time.time() # Start the timer for benchmarking execution time
    #Argument parser for various command line arguments
    parser = argparse.ArgumentParser(description='''This program extracts TRS sequences from a series of input genomes,
                                     allows for length selection of extracted fragments, and prepares sequences for further analysis.''')
    # Adding various arguments some currently unused
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
    parser.add_argument('--taxonomy_db',help='Path to the taxonomy - accession database', required=False, type= str)
    parser.add_argument('--update_taxonomy', help="Download and overwrite outdated taxonomy file", required=False, action='store_true')
    #Mutually exclusive arguments for NCBI IDs or file containing the,
    add_ids = parser.add_mutually_exclusive_group()
    add_ids.add_argument('--ids_to_add_to_dictionary', help= "Comma separated list or single value of NCBI IDs", required= False, type=str)
    add_ids.add_argument('--ids_file', help="Path to the file containing NCBI IDs", required= False, type= str)
    #Same for exceptions (taxid is used here)
    exceptions = parser.add_mutually_exclusive_group()
    exceptions.add_argument('--taxids_to_add_to_exceptions', help="Comma separated list or single taxid to add to filtering exceptions", required= False,type=str)
    exceptions.add_argument('--taxids_file', help= "Path to the file containing the taxids that should be added to exceptions", required= False, type=str)
    #Path to blast db, has to be nt 
    parser.add_argument('--blast_db', help= "Path to your blast_db", required=True, type=str)
    args = parser.parse_args()

    #Determine the taxonomy path
    if not args.taxonomy_db:
        default_taxonomy_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxonomy")
        taxonomy_db_path = default_taxonomy_path
    else:
        taxonomy_db_path = args.taxonomy_db

    print("Checking taxonomy file...")
    check_and_update_taxonomy_file(taxonomy_db_path, args.update_taxonomy)


    #Assing parsed arguments to variables for use downstream
    taxonomy_db = args.taxonomy_db
    ids_to_add_to_dictionary = args.ids_to_add_to_dictionary
    ids_file = args.ids_file
    taxids_to_add_to_exceptions = args.taxids_to_add_to_exceptions
    taxids_file = args.taxids_file


    #Validate and set the email used for querying NCBI services 
    Entrez.email = SequenceProcessor.validate_and_set_email(args.email)
    print(f"Email adress is currently set to {Entrez.email}")

    # Setup the base directory for storing results
    input_fasta_folder_path_name = os.path.basename(args.input_fasta_folder_path)
    base_results_directory = os.path.join(os.getcwd(), f"{input_fasta_folder_path_name}_results")
    results_directory = base_results_directory
    FileHandler.ensure_directory_exists(base_results_directory) #Create the directory if it does not exist inform user if it does

    # Define the name of the CSV file that will store the results of the analysis 
    name_of_csv_file_storing_TRS_analysis_results = input_fasta_folder_path_name + "_results.csv"
    path_of_folder_storing_TRS_analysis_results = os.path.join(base_results_directory, "TRS_output")
    FileHandler.ensure_directory_exists(path_of_folder_storing_TRS_analysis_results)
    path_of_csv_file_storing_TRS_analysis_results = os.path.join(path_of_folder_storing_TRS_analysis_results, 
                                                                 name_of_csv_file_storing_TRS_analysis_results)
    '''Cont argument use is currently not implemented'''
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
        # Check if the provided input directory exists if it does the program will quit
        if not os.path.exists(args.input_fasta_folder_path):
            print("Provided directory path does not exist. Quitting...")
            sys.exit(1)
        print("Existing directory provided. Proceeding....")
        # Find all FASTA files in the input directory 
        fasta_files = [f for f in os.listdir(args.input_fasta_folder_path) if f.endswith('.fasta') or f.endswith('.fa')]
        if not fasta_files:
            print("No fasta files found in the provided directory. Quitting...")
            sys.exit(1)

        # If results already exist and --redo is not specified, prompt to exit
        if os.path.exists(path_of_csv_file_storing_TRS_analysis_results) and not args.redo:
            print("Results already exist. Use --redo to overwrite or --cont to continue from existing results. Exiting program.")
            sys.exit(1)

        print("Starting analysis...")   

        trs_calculators = []

        # Iterate over each FASTA file to calculate TRS sequences
        for fasta_file in fasta_files:
            path_to_input_fasta = os.path.join(args.input_fasta_folder_path, fasta_file)
            print(f"{path_to_input_fasta}")
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

        list_of_trs_results = [] #Store results

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

    #Adjust the length of flanking sequences to extract so that they don't overlap
    l_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    r_chars = SequenceProcessor.adjust_input_to_range(args.length_to_extract, args.tmin)
    l_chars = int(l_chars)
    r_chars = int(r_chars)
    print(combined_trs_results)
    
    #Extract flanking sequences and update the results DataFrame
    combined_trs_results = SequenceProcessor.extract_sequences(combined_trs_results, l_chars, r_chars)

    #Update the results directory after extraction of flanks
    results_directory_after_flanks_extracted = f"{results_directory}_L{l_chars}_R{r_chars}"
    results_directory_after_flanks_extracted_path = os.path.join(os.path.dirname(results_directory), results_directory_after_flanks_extracted)

    print(f"Results directory is set to: {results_directory_after_flanks_extracted_path}") #WE ARE HERE

    # Rename the results directory if necessary
    if not args.cont:
        if not os.path.exists(results_directory_after_flanks_extracted_path):
            # Rename the existing results directory
            os.rename(results_directory, results_directory_after_flanks_extracted_path)
            print(f"The results directory has been renamed to: {results_directory_after_flanks_extracted_path}")
        else:
            print(f"Directory {results_directory_after_flanks_extracted_path} already exists. Consider using a different name or removing the existing directory.")

    results_directory = results_directory_after_flanks_extracted_path

    #Fetch organism names from NCBI based on genome IDs
    ncbi_ids = combined_trs_results["GENOME"].unique().tolist()
    if Entrez.email and Entrez.email == args.email :
         print(f"Email adress is still set to {Entrez.email}")
    organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids,email = Entrez.email)

    # Map genome IDs to taxonomic names
    combined_trs_results['Taxonomic Name'] = None
    combined_trs_results['Taxonomic Name'] = combined_trs_results['GENOME'].map(organism_map)
    
    # Identify and warn if unmatched genomes were found in the dataset
    unmatched_genomes = combined_trs_results[combined_trs_results['Taxonomic Name'].isnull()]["GENOME"].unique()
    if len(unmatched_genomes) > 0:
        print(f"Warning: Some genome IDs could not be matched with taxonomic names: {unmatched_genomes}")
    
    # Generate sequence IDs for left and right flanking sequences 
    combined_trs_results['L_id'] = combined_trs_results['Taxonomic Name'] + '_L' + combined_trs_results['L-No'].astype(str)
    combined_trs_results['R_id'] = combined_trs_results['Taxonomic Name'] + '_R' + combined_trs_results['R-No'].astype(str)
    
    #Prepare a dataframe containing flanking sequneces and their IDs.
    sequences_df = combined_trs_results[['SEQ_L', 'SEQ_R', 'L_id', 'R_id']]
    print(sequences_df)
    
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

    #New way to filter clusters
    input_fasta = os.path.join(results_directory,"TRS_output","combined_sequences_unique.fasta")
    output_folder_cluster = os.path.join(cd_hit_results_folder,"fasta_clusters")
    FileHandler.create_fasta_for_all_clusters(cdhit_clusters_file, input_fasta, output_folder_cluster,create_individual_files=False)

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
    
    # Paths to SLURM and fallback script
    slurm_script = 'TRS_BLAST.sh'
    fallback_script = 'TRS_BLAST_NOSLURM.sh'

    # Set the paths to the query directory and BLAST database
    query_dir = sequences_after_clusters_filtering_folder  
    blast_db = args.blast_db
    perc_identity = 100 #should probably lower that to 90 or less - first pass let's keep it at 100

    # Set the environment variables before running the SLURM script
    env = os.environ.copy()
    env['QUERY_DIR'] = query_dir
    env['BLAST_DB'] = blast_db
    env['OUTPUT_DIR'] = blast_results_folder
    env['PERC_IDENTITY'] = str(perc_identity)


    try:
        # Attempt to submit the SLURM job using sbatch with the --wait option
        subprocess.run(['sbatch', '--wait', slurm_script], env=env, check=True
                                      ,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except:
        print(f"SLURM job submission failed")

    
            # If SLURM job submission fails, run the fallback script and wait for its completion
        try:
            print(f"Running fallback script from: {fallback_script}")
            result = subprocess.run([fallback_script, query_dir, blast_db], check=True, env=env)

            if result.returncode == 0:
                print("Fallback script completed successfully. Continuing with the Python script.")
            else:
                print("Fallback script failed.", file=sys.stderr)
                sys.exit(1)
        except FileNotFoundError:
            print(f"Fallback script {fallback_script} not found. Please check the path.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Fallback script failed with error: {e.stderr.decode()}")
            sys.exit(1)

    
    #Process BLAST analysis results
    blast_output_path = blast_results_folder
    FileHandler.convert_to_txt(blast_output_path) #convert output to .txt format
    FileHandler.filter_and_overwrite_blast_file(blast_output_path) # Filter and overwrite BLAST results file
    modified_blast_path = os.path.join(blast_output_path, "modified_blast")
    FileHandler.ensure_directory_exists(modified_blast_path)
    

    #Collect accessions from BLAST result files
    accessions = BLASTProcessor.collect_accessions_from_blast_files(blast_output_path)
    
    # Filter taxonomy file based on collected accessions
    taxonomy_db = os.path.join(taxonomy_db,"nucl_gb.accession2taxid.gz")
    tax_df = SequenceProcessor.filter_taxonomy_file(taxonomy_db, accessions)
    print("Collected Accessions:", accessions)
    print("Filtered Taxonomy DataFrame:\n", tax_df)

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

    # Append taxids as a new column to BLAST result files
    print(f"Appending a new column with taxids to the blast files.....")
    BLASTProcessor.match_accessions_to_taxids_and_add_column(
        blast_output_path,
        modified_blast_path,
        lambda accession: BLASTProcessor.map_accession_to_taxid(accession, accession_to_taxid)
    )

    end_time = time.time() # End the timer for benchmarking execution time
    print(f"Processing completed in {end_time - start_time:.2f} seconds")

    #Search for the results CSV file in the results directory (created during TRS step)
    print(f"Searching for results file...")
    results_file = FileHandler.find_file_by_name('*_results.csv', folder= results_directory)
    if results_file:
        result_file = results_file[0]
        print(f"Results file found at : {results_file}")
    else:
        print("No '*_results.csv' file selected or found.")

    #Load the combined results into a DataFrame
    combined_results = pd.read_csv(result_file)

    #Set the Entrez e-mail for further processing - seems like for some reason it refuses to remember it a this step and does not pass it into future ncbi requests
    Entrez.email = SequenceProcessor.validate_and_set_email(args.email)
    print(f"Current email is set to {Entrez.email}")

    #Fetch taxonomic information for the collected genome IDs
    ncbi_ids = combined_results["GENOME"].unique().tolist()
    tax_map = SequenceProcessor.fetch_organism_taxids(ncbi_ids)

    #Filter and clean the taxonomic mapping
    filtered_organism_taxid_map = {SequenceProcessor.filter_key_parts(key): value for key, value in tax_map.items()}
    print(f"Species - taxid pairs detected in dataset : {filtered_organism_taxid_map}")

    #Append taxid information to the filtered map
    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)


    #Interact and update the taxid-accession map using the provided NCBI IDs or file
    BLASTProcessor.interact_and_update_dict(filtered_organism_taxid_map,ncbi_ids_list=ids_to_add_to_dictionary,file_path=ids_file)

    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
    species_info = {SequenceProcessor.filter_key_parts(key): value for key, value in species_info.items()}

    # Prepare NaN file path for storing the results with missing acessions
    '''The blast db or taxonomy database is out of date if there are entries here, the problem is that the internal ncbi database gets updated much more
    frequently that the ones available for download, to remedy that atleast partially some way to user database to available online could be introduced'''
    nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")

    # Construct a dictionary of all taxid - acessions pairs in our data
    results_dict = BLASTProcessor.construct_dict_from_files(modified_blast_path,nan_file)


    #Ask for exceptions to filtering based on taxids - exceptions will not be consider as disquilifying when they are found as matches for a sequence
    '''Large pool of exceptions my lead to untrustworthy results, however they should get rectified in the 2nd pass'''
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

    # '''Separate reads in and outside of clusters into two categories 1) reads that have a matching pair (twins) 2) reads that do not (singles)'''
    # processed_fasta_file_path = FileHandler.find_file_by_name(file_name="unique_taxids_in_clusters_combined_sequences_unique_blastn_out.txt",folder=results_directory)
    # print(f"{processed_fasta_file_path}")
    # processed_fasta_file_path = pathlib.Path(processed_fasta_file_path[0]).parent

    # #Separate the sequences into "singles" and "twins" based on clustering results
    # BLASTProcessor.separate_into_singles_and_twins(processed_fasta_file_path,file_pattern="*.txt")

    # #Read the sequence IDs from the processed FASTA file and classify them into a valid dictionary
    # '''Create a special_dict that contains all sequence IDs that occur either in clusters or have no pair, the sequence_ids_dict contains the remaining IDs'''
    # sequence_ids_dict, special_dict = BLASTProcessor.read_sequence_ids(processed_fasta_file_path)

    #Find the filtered FASTA file that contains all sequences(but unique)
    filtered_fasta_file = FileHandler.find_file_by_name('unique_taxids_not_in_clusters_combined_sequences_unique_blastn_out.txt',folder= results_directory)
    filtered_fasta_file = pathlib.Path(filtered_fasta_file[0])

    #Ensure the output directory exists for storing of filtred results 
    intermediate_results = os.path.join(results_directory,"intermediate_results")
    FileHandler.ensure_directory_exists(intermediate_results)
    #fasta_to_parse = BLASTProcessor.filter_fasta_file_dict(filtered_fasta_file,sequence_ids_dict,special_dict,output_directory)
    BLASTProcessor.extract_full_TRS_sequences(combined_results,filtered_fasta_file,results_directory,state=1)

    #Stat modules
    import matplotlib.pyplot as plt
    import re
    import numpy as np
    #Change to instead of filtering output the desired names (edit file_paths)
    fasta_files_for_plotting = FileHandler.search_for_files(results_directory,'*.fasta')
    fasta_files_for_plotting_names = FileHandler.extract_file_names(fasta_files_for_plotting)
    print(f"Following fasta files were found in {results_directory} : {fasta_files_for_plotting_names}")
    statistics = Stats()
    for fasta_path in fasta_files_for_plotting:
        statistics.count_L_R(fasta_path)
    file_paths = [file_path for file_path in fasta_files_for_plotting if file_path in statistics.file_info]
    names_to_filter = ["combined_sequences_unique.fasta","full_sequences.fasta"]
    file_paths = [file_path for file_path in file_paths if os.path.basename(file_path) in names_to_filter]
    print(f"{file_paths}")
    # Dictionary to store titles for each file path
    file_titles = {}

    # Populate file_titles dictionary using file_paths list
    for file_path in file_paths:
    # Extract file name from file path
        file_name = os.path.basename(file_path)

    # Construct default title based on file name
    default_title = f"Title for {file_name}"

    # Add entry to file_titles dictionary
    file_titles[file_path] = default_title

    statistics.plot_lr_counts(file_paths=file_paths,file_titles=file_titles,results_directory=results_directory)

    txt_files_for_plotting = FileHandler.search_for_files(results_directory,'*.txt')
    txt_files_for_plotting_name = FileHandler.extract_file_names(txt_files_for_plotting)
    print(f"Following fasta files were found in {results_directory} : {txt_files_for_plotting_name}")
    for txt_path in txt_files_for_plotting:
        statistics.count_L_R(txt_path)
    file_paths = [file_path for file_path in txt_files_for_plotting if file_path in statistics.file_info]
    names_to_filter = ["unique_taxids_not_in_clusters_combined_sequences_unique_blastn_out.txt"]
    file_paths = [file_path for file_path in file_paths if os.path.basename(file_path) in names_to_filter]
    #Dictionary to store titles for each file path
    file_titles = {}

    for file_path in file_paths:
    # Extract file name from file path
        file_name = os.path.basename(file_path)

    # Construct default title based on file name
    default_title = f"Title for {file_name}"

    # Add entry to file_titles dictionary
    file_titles[file_path] = default_title



    statistics.plot_lr_counts(file_paths=file_paths,file_titles=file_titles,results_directory=results_directory)

    # Paths to SLURM and fallback script
    slurm_script = 'TRS_BLAST.sh'
    fallback_script = 'TRS_BLAST_NOSLURM.sh'

    # Set the paths to the query directory and BLAST database
    query_dir = intermediate_results  
    blast_db = args.blast_db
    perc_identity = 90 
    blast_results_folder_2nd_pass = os.path.join(results_directory,"blast_output_2nd_pass")
    

    # Set the environment variables before running the SLURM script
    env = os.environ.copy()
    env['QUERY_DIR'] = query_dir
    env['BLAST_DB'] = blast_db
    env['OUTPUT_DIR'] = blast_results_folder_2nd_pass
    env['PERC_IDENTITY'] = str(perc_identity)
    FileHandler.ensure_directory_exists(blast_results_folder_2nd_pass)

    try:
        # Attempt to submit the SLURM job using sbatch with the --wait option
        subprocess.run(['sbatch', '--wait', slurm_script], env=env, check=True
                                      ,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    except:
        print(f"SLURM job submission failed")

    
            # If SLURM job submission fails, run the fallback script and wait for its completion
        try:
            print(f"Running fallback script from: {fallback_script}")
            result = subprocess.run([fallback_script, query_dir, blast_db], check=True, env=env)

            if result.returncode == 0:
                print("Fallback script completed successfully. Continuing with the Python script.")
            else:
                print("Fallback script failed.", file=sys.stderr)
                sys.exit(1)
        except FileNotFoundError:
            print(f"Fallback script {fallback_script} not found. Please check the path.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"Fallback script failed with error: {e.stderr.decode()}")
            sys.exit(1)

    #Process BLAST analysis results
    blast_output_path = blast_results_folder_2nd_pass
    FileHandler.convert_to_txt(blast_output_path) #convert output to .txt format
    FileHandler.filter_and_overwrite_blast_file(blast_output_path) # Filter and overwrite BLAST results file
    modified_blast_path = os.path.join(blast_output_path, "modified_blast")
    FileHandler.ensure_directory_exists(modified_blast_path)
    

    #Collect accessions from BLAST result files
    accessions = BLASTProcessor.collect_accessions_from_blast_files(blast_output_path)
    
    # Filter taxonomy file based on collected accessions
    #taxonomy_db = os.path.join(taxonomy_db,"nucl_gb.accession2taxid.gz") we already did that 
    tax_df = SequenceProcessor.filter_taxonomy_file(taxonomy_db, accessions)
    print("Collected Accessions:", accessions)
    print("Filtered Taxonomy DataFrame:\n", tax_df)

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

    # Append taxids as a new column to BLAST result files
    print(f"Appending a new column with taxids to the blast files.....")
    BLASTProcessor.match_accessions_to_taxids_and_add_column(
        blast_output_path,
        modified_blast_path,
        lambda accession: BLASTProcessor.map_accession_to_taxid(accession, accession_to_taxid)
    )


    #Filter and clean the taxonomic mapping
    filtered_organism_taxid_map = {SequenceProcessor.filter_key_parts(key): value for key, value in tax_map.items()}
    print(f"Species - taxid pairs detected in dataset : {filtered_organism_taxid_map}")

    #Append taxid information to the filtered map
    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)


    #Interact and update the taxid-accession map using the provided NCBI IDs or file
    BLASTProcessor.interact_and_update_dict(filtered_organism_taxid_map,ncbi_ids_list=ids_to_add_to_dictionary,file_path=ids_file)

    species_info = BLASTProcessor.append_taxids_to_filtered_map(filtered_organism_taxid_map)
    species_info = {SequenceProcessor.filter_key_parts(key): value for key, value in species_info.items()}

    # Prepare NaN file path for storing the results with missing acessions
    '''The blast db or taxonomy database is out of date if there are entries here, the problem is that the internal ncbi database gets updated much more
    frequently that the ones available for download, to remedy that atleast partially some way to user database to available online could be introduced'''
    nan_file = os.path.join(modified_blast_path,"NaN acessions.csv")

    # Construct a dictionary of all taxid - acessions pairs in our data
    results_dict = BLASTProcessor.construct_dict_from_files(modified_blast_path,nan_file)


    #Ask for exceptions to filtering based on taxids - exceptions will not be consider as disquilifying when they are found as matches for a sequence
    '''Large pool of exceptions my lead to untrustworthy results, however they should get rectified in the 2nd pass'''
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

    
    #Perform filtering based on the exceptions and matching taxid criteria
    print("Filtering keys....")
    filtered_keys = BLASTProcessor.filter_with_exceptions(results_dict,species_info,exceptions)
    filtered_keys_final = BLASTProcessor.unpack_single_element_sets(filtered_keys)

    #Process the filtered BLAST result files 
    print("Filtering files...")
    FileHandler.process_files_with_filter(modified_blast_path,filtered_keys_final)

    #Save the file from this into final_results or smth

    #Find the filtered FASTA file that contains all sequences(but unique)
    filtered_fasta_file = FileHandler.find_file_by_name('unique_taxids_full_sequences_blastn_out.txt',folder= blast_results_folder_2nd_pass)
    filtered_fasta_file = pathlib.Path(filtered_fasta_file[0])

    # #Ensure the output directory exists for storing of filtred results 
    output_directory = os.path.join(results_directory,"final_results")
    FileHandler.ensure_directory_exists(output_directory)
    BLASTProcessor.extract_full_TRS_sequences(combined_results,filtered_fasta_file,results_directory,state=2)


    fasta_files_for_plotting = FileHandler.search_for_files(results_directory,'*.fasta')
    fasta_files_for_plotting_names = FileHandler.extract_file_names(fasta_files_for_plotting)
    print(f"Following fasta files were found in {results_directory} : {fasta_files_for_plotting_names}")
    statistics = Stats()
    for fasta_path in fasta_files_for_plotting:
        statistics.count_L_R(fasta_path)
    file_paths = [file_path for file_path in fasta_files_for_plotting if file_path in statistics.file_info]
    names_to_filter = ["full_sequences_final.fasta"]
    file_paths = [file_path for file_path in file_paths if os.path.basename(file_path) in names_to_filter]
   
    print(f"{file_paths}")
    # Dictionary to store titles for each file path
    file_titles = {}

    # Populate file_titles dictionary using file_paths list
    for file_path in file_paths:
    # Extract file name from file path
        file_name = os.path.basename(file_path)

    # Construct default title based on file name
    default_title = f"Title for {file_name}"

    # Add entry to file_titles dictionary
    file_titles[file_path] = default_title

    statistics.plot_lr_counts(file_paths=file_paths,file_titles=file_titles,results_directory=results_directory)

    txt_files_for_plotting = FileHandler.search_for_files(results_directory,'*.txt')
    txt_files_for_plotting_name = FileHandler.extract_file_names(txt_files_for_plotting)
    print(f"Following fasta files were found in {results_directory} : {txt_files_for_plotting_name}")
    for txt_path in txt_files_for_plotting:
        statistics.count_L_R(txt_path)
    file_paths = [file_path for file_path in txt_files_for_plotting if file_path in statistics.file_info]
    names_to_filter = ["unique_taxids_full_sequences_blastn_out.txt"]
    file_paths = [file_path for file_path in file_paths if os.path.basename(file_path) in names_to_filter]
    #Dictionary to store titles for each file path
    file_titles = {}

    for file_path in file_paths:
    # Extract file name from file path
        file_name = os.path.basename(file_path)

    # Construct default title based on file name
    default_title = f"Title for {file_name}"

    # Add entry to file_titles dictionary
    file_titles[file_path] = default_title

    statistics.plot_lr_counts(file_paths=file_paths,file_titles=file_titles,results_directory=results_directory)

    #BINOMIAL TEST
    full_seq_final = os.path.join(output_directory,"full_sequences_final.fasta")
    final_results_folder = output_directory
    df = Stats.create_trs_class_dataframe(output_folder_cluster)
    df = Stats.get_trs_class_totals(df) 
    df_full_fasta = statistics.create_lr_counts_dataframe(file_path=full_seq_final)
    merged_df = Stats.merge_trs_classes_results(df_full_fasta,df)
    result_df = Stats.test_binomial_for_classes(merged_df)
    result_df.index.name = "Class"
    Stats.save_binom_results(result_df,final_results_folder)


if __name__ == "__main__":
    main()
