import os
import fnmatch
import pandas as pd
from src.SequenceProcessor import SequenceProcessor
from Bio import Entrez
import pathlib

class BLASTProcessor:

    @staticmethod
    def process_blast_files_in_directory(directory_path, file_pattern="*.txt"):
        """
        Processes all BLAST output files in a specified directory matching a give pattern.
        Each file is filtered to remove duplicate pairs of sequence IDs and accession numbers.
        
        Args:
            directory_path (str): Path to the directory containing BLAST files.
            file_pattern (str): Pattern to match BLAST files (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                BLASTProcessor.filter_and_overwrite_blast_file(file_path)

    @staticmethod
    def filter_and_overwrite_blast_file(file_path):
        """
        Filters a BLAST output file by removing duplicate pairs of sequence IDs and accession numbers.
        Only the first occurrence of each unique pair is kept, and the file is overwritten.
        
        Args:
            file_path (str): Path to the BLAST file to be processed.
        """
        unique_pairs = set() # Store unique pairs of sequence IDs and accession IDs
        lines_to_keep = [] # Store lines to be retained

        try:
            with open(file_path, "r") as file:
                for line in file:
                    try:
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:
                            pair = (columns[0], columns[1]) # Get the pair
                            if pair not in unique_pairs:
                                unique_pairs.add(pair)
                                lines_to_keep.append(line)
                    except Exception as e:
                        print(f"An unexpected error has occurred while processing a line in {file_path}: {e}")
                print(f"Searching for unique sequence ID - accession pairs in {file_path}....")
        except PermissionError as e:
            print(f"Permission denied while trying to access {file_path}: {e}")
            return
        except IOError as e:
            print(f"An error has occurred while opening or reading {file_path}: {e}")
            return
        
        #Overwrite the file with filtered content
        try:
            with open(file_path, "w") as file:
                file.writelines(lines_to_keep)
        except PermissionError as e:
            print(f"Permission denied while accessing {file_path} for writing: {e}")
        except IOError as e:
            print(f"An error occurred while writing to {file_path}: {e}")

    @staticmethod
    def collect_accessions_from_blast_files(directory_path, file_pattern="*.txt"):
        """
        Collects accession numbers from BLAST files in a specified directory.
        
        Args:
            directory_path (str): Path to the directory containing BLAST files.
            file_pattern (str): Pattern to match BLAST files (default is "*.txt").
            
        Returns:
            set: A set of unique accession numbers extracted from BLAST files.
        """
        accessions = set()
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                with open(file_path, "r") as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        if len(columns) >= 2:
                            accessions.add(columns[1])
        return accessions

    @staticmethod
    def filter_and_overwrite_files_in_directory(directory_path, file_pattern="*.txt"): # i don't remember why i pretty much duplicated this function and at this point im too afraid to remove it
        """
         Filters and overwrites BLAST files in a directory by removing duplicates, similar to process_blast_files_in_directory.
        
        Args:
            directory_path (str): Path to the directory containing files to filter.
            file_pattern (str): Pattern to match files (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                BLASTProcessor.filter_and_overwrite_blast_file(file_path)
    @staticmethod
    def map_accession_to_taxid(accession, taxid_dict):
        """
        Maps an accession to a taxid using provided dictionary (in practice it iterates throguh multiple accessions).

        Args:
            accession(str): The accession to be mapped.
            taxid_dict(dict): A dictionary mapping accessions to taxids

        Returns:
            str: The taxid corresponding to acession or empty string if no match was found
        """
        return taxid_dict.get(accession, '')
    @staticmethod
    def match_accessions_to_taxids_and_add_column(blast_output_path, modified_blast_path, map_func):
        """
        Maps accessions to taxids in BLAST output files and adds them as additional column.

        Args:
            blast_output_path(str): Path to the directory containing the blast files that will be modified.
            modified_blast_path(str): Path to the directory in which results of this operation will be stored.
            map_func(function): A function that maps acession numbers to taxids.
        """
        for filename in os.listdir(blast_output_path):
            if filename.endswith(".txt"):
                input_file_path = os.path.join(blast_output_path, filename)
                
                # Read the file into a DataFrame
                df = pd.read_csv(input_file_path, sep='\t', header=None)  # No column names specified
                
                # Map accessions to TaxIDs and add them as a new column after the last existing column
                df[df.shape[1]] = df[1].map(lambda x: map_func(x))  # Apply the mapping function to column 1
                
                # Construct output file path
                modified_file_path = os.path.join(modified_blast_path, f"taxids_{filename}")
                
                # Save the modified DataFrame to a new file in the output directory
                df.to_csv(modified_file_path, sep='\t', index=False, header=False)
    @staticmethod
    def append_taxids_to_filtered_map(filtered_map):
        """
        Adds taxids to a dictionary of organism names by querying NCBI taxonomy.

        Args:
            filtered_map(dict): A dictionary mapping organism names to taxids.

        
        Returns:
            dict: Update dictionary where organism names are mapped to sets of taxids 
        """
        updated_map = {}
        for organism_name, original_taxid in filtered_map.items():
            # Makes is so the values are stored in sets 
            if organism_name not in updated_map:
                updated_map[organism_name] = set()
            updated_map[organism_name].add(original_taxid)

            try:
                handle = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
                search_results = Entrez.read(handle)
                handle.close()
                taxids = search_results.get("IdList",[])
                if taxids:
                    #Assume the first result is the most relevant
                    new_taxid = taxids[0]
                    if new_taxid not in updated_map[organism_name]:
                        updated_map[organism_name].add(new_taxid)
                        print(f"Updated {organism_name} with additional TaxID: {taxids[0]}")
                    else:
                        print(f"TaxID {new_taxid} already present for {organism_name}, skipping....")
                else:
                    print(f"No additional TaxID found for: {organism_name}")
            except Exception as e:
                print(f"Error fetching additional TaxID for {organism_name}: {e}")
        return updated_map
    @staticmethod
    def interact_and_update_dict(dictionary, ncbi_ids_list = None, file_path = None):
        """
        Updates a dictionary with a new NCBI ID specified by the user in the arguments

        Args:
            dictionary(dict): The existing dictionary mapping organism names to SETS of TaxIDs.
            user_input(str): Comma-separated list or single value of NCBI IDs.
            file_path(str): Path to file containing NCBI IDs.
        """
        if ncbi_ids_list:
            ids = [id.strip() for id in ncbi_ids_list.split(",")]
            if ids:
                new_organism_ids = SequenceProcessor.fetch_organism_taxids(ids)
                for organism_name, taxid_set in new_organism_ids.items():
                    if organism_name in dictionary:
                        print(f"Organism {organism_name} already present in dictionary,updating dictionary")
                        dictionary[organism_name].update(taxid_set)
                    else:
                        print(f"New organism detected : {organism_name}")
                        dictionary[organism_name] = taxid_set
        elif file_path:
            try:
                with open(file_path,'r') as file:
                    ids = [line.strip() for line in file if line.strip()]
                    new_organism_ids = SequenceProcessor.fetch_organism_taxids(ids)
                    for organism_name, taxid_set in new_organism_ids.items():
                        if organism_name in dictionary:
                            print(f"Organism {organism_name} already present in dictionary,updating dictionary")
                            dictionary[organism_name].update(taxid_set)
                        else:
                            print(f"New organism detected : {organism_name}")
                            dictionary[organism_name] = taxid_set
            except FileNotFoundError:
                print("The specified file was not found, check path and try again")
            except Exception as e:
                print("An unexpected error has occured: {e}")
        else:
            print("No input provided nothing will be added")
    @staticmethod
    def ask_for_exception(exception_ids = None, file_path = None):
        """
        Adds exceptions to further filtering steps based on provided arguments(taxids are used here).

        Args: 
            exception_ids(str): Comma separated list or single value of taxids.
            file_path(str): Path to file containing the exceptions.
        
        Returns:
            set: A set of taxid exceptions
        """
        exceptions = set()

        if exception_ids:
            exceptions = SequenceProcessor.process_taxid_input(exception_ids)
        elif file_path:
            try:
                with open(file_path, 'r') as file:
                    file_exceptions = set(line.strip() for line in file if line.strip().isdigit())
                    exceptions.update(file_exceptions)
            except FileNotFoundError:
                print("The specified file was not found, check path and try again")
            except Exception as e:
                print(f"An unexpected error has occurred: {e}")
        else:
            print("No input provided. No exceptions will be added.")

        return exceptions
    
    @staticmethod
    def filter_with_exceptions(results_dict, species_info, exceptions):
        """
        Filters the results dictionary based on the presence of taxids in the species info dictionary, taking into account exceptions.

        Args:
            results_dict(dict): Dictionary of results to filter.
            species_info(dict): Dictionary of species taxids.
            exceptions(set): Set of taxids to exclude.

        Returns:
            dict: Filtered results dictionary.
        """
        filtered_keys = {}
        all_species_taxids = set.union(*species_info.values())  # Combine all taxids from species info for easier lookup

        for key, values in results_dict.items():
            non_exception_values = values - exceptions  # Remove exceptions 
            if non_exception_values:  # If there are none left
                if len(non_exception_values) == 1 and non_exception_values & all_species_taxids:
                    filtered_keys[key] = non_exception_values
            else:  # If all values are exceptions 
                if values & all_species_taxids:
                    filtered_keys[key] = values
        return filtered_keys
    @staticmethod
    def unpack_single_element_sets(input_dict):
        """
        Unpacks single-element sets in the input dictionary replacing them with the element itself.

        Args:
            input_dict(dict): The dictionary whose single-element sets are to be unpacked.

        Returns:
            dict: Dictionary with single-element sets unpacked.
        """
        processed_dict = {}

        for key, value in input_dict.items():
            # Check if the value is a set with exactly one entry
            if isinstance(value, set) and len(value) == 1:
                # "Unpack" the set, storing its single element as the value
                processed_dict[key] = next(iter(value))
            else:
                # If the value is not a single-element set, retain it as is
                processed_dict[key] = value

        return processed_dict
    @staticmethod
    def separate_into_singles_and_twins(directory_path, file_pattern="*.txt"):
        """
        Separates sequences in files into those with paired identifiers (twins) and single identifiers
        Writes single entries to a new file named '<original_filename>_singles.txt'.

        Args:
            directory_path (str): Path to the directory containing the files to process.
            file_pattern (str): Pattern of the file names to process (default is "*.txt").
        """
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename, file_pattern):
                file_path = os.path.join(directory_path, filename)
                single_entries_path = os.path.join(directory_path, f"{os.path.splitext(filename)[0]}_singles.txt")

                print(f"Processing file: {filename}")  # Debugging

                with open(file_path, 'r') as file:
                    lines = file.readlines()

                print(f"Total lines read from {filename}: {len(lines)}")  # Debugging

                number_count = {}
                filtered_entries = []
                non_unique_entries = []

                # First pass to count occurrences of each identifier
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        number_count[number] = number_count.get(number, 0) + 1

                print(f"Number count: {number_count}")  # Debugging

                # Second pass to separate unique and non-unique entries
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        if number_count[number] > 1:
                            filtered_entries.append(line.strip())
                        else:
                            non_unique_entries.append(line.strip())

                print(f"Filtered entries for {filename}: {len(filtered_entries)}")  # Debugging
                print(f"Non-unique entries for {filename}: {len(non_unique_entries)}")  # Debugging

                # Overwrite the original file with paired identifiers
                with open(file_path, 'w') as file:
                    for entry in filtered_entries:
                        file.write(entry + '\n')

                print(f"Overwritten {filename} with filtered entries.")  # Debugging

                # Write out single entries to a separate file
                with open(single_entries_path, 'w') as file:
                    for entry in non_unique_entries:
                        file.write(entry + '\n')

                print(f"Written non-unique entries to {single_entries_path}.")  # Debugging

    @staticmethod
    def construct_dict_from_files(directory, nan_file):
        '''
        Construct dictionaries containing the sequence ids as keys and taxids as values from the files in the directory

        Args:
            directory(str): Path to files from which to construct dictionary 
            nan_file(str): Path to the file where keys without values will be stored
        
        Returns:
            dict: Dictionary mapping sequence IDs to taxids
        '''
        data_dict = {}
        nan_keys = set()
        for filename in os.listdir(directory):
            if filename.endswith('.txt'):
                filepath = os.path.join(directory,filename)
                # determine the number of columns
                with open(filepath,"r") as f:
                    first_line = f.readline()
                    num_columns = len(first_line.split('\t'))
                
                #Read only necessary columns (1st and last)
                df = pd.read_csv(filepath, sep='\t',header=None, usecols=[0,num_columns-1])

                #Group by the first column and convert the last one to sets
                grouped = df.groupby(0)[num_columns-1].apply(set).to_dict()

                #Merge the current file's dictionary with the main dictionary
                for key, value_set in grouped.items():
                    '''
                    This is a bandaid fix for a problem that can be encountered when either taxonomy or blast database is out of date 
                    '''
                    value_set = {val for val in value_set if pd.notna(val)}
                    if len(value_set) < len(grouped[key]):
                        nan_keys.add(key) # Add key which had NaN
                    if key in data_dict:
                        data_dict[key].update(value_set)
                    else:
                        data_dict[key] = value_set
        with open(nan_file,"w") as nf:
            for key in nan_keys:
                nf.write(key + '\n')
        if nan_keys:
            print("NaN values were found while constructing dictionary you blast db or taxonomy db might be out of date." + "\n"
                  + f"Keys associated with NaN values were saved to: {nan_file} analysis can proceed without those values.")
        else:
            print("No NaN values detected")
        return data_dict
    @staticmethod
    def read_sequence_ids(directory_path):
        """
        Reads sequence IDs from files in the specified directory and organizes them into two dictionaries.
        
        Arguments:
        directory_path : str or Path : Path to the directory containing the sequence ID files.
        
        Returns:
        tuple : (sequence_ids_dict, special_dict)
            - sequence_ids_dict: A dictionary containing the remaining sequence IDs (not in clusters paired).
            - special_dict: A dictionary containing the sequence IDs from sequences in clusters or those that are not in clusters but don't have a pair.
        """
        sequence_ids_dict = {}
        special_dict = {}

        directory = pathlib.Path(directory_path)
        for file_path in directory.glob('unique_*'):
            with file_path.open('r') as file:
                sequence_ids = [line.split()[0] for line in file if line.strip()]
            
            if sequence_ids:
                if file_path.name == 'unique_taxids_in_clusters_combined_sequences_unique_blastn_out.txt' or file_path.name =="unique_taxids_in_clusters_combined_sequences_unique_blastn_out_singles.txt":
                    special_dict['special'] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
                else:
                    sequence_ids_dict[file_path.stem] = {seq_id.split('_')[-1] for seq_id in sequence_ids}
            else:
                print(f'The file {file_path.name} is empty.')
        return sequence_ids_dict,special_dict
    @staticmethod
    def filter_fasta_file_dict(fasta_file_path, sequence_ids_dict, special_dict, output_directory):
        """
        Filters sequences from a FASTA file based on the provided sequence IDs and writes the filtered sequences to an output file.
        If a sequence is part of the special IDs, its header is modified with a 'CLUSTER_' prefix.
        
        Arguments:
            fasta_file_path : str or Path : Path to the input FASTA file.
            sequence_ids_dict : dict : Dictionary of sequence IDs to filter (from read_sequence_ids function).
            special_dict : dict : Dictionary containing special sequence IDs (from read_sequence_ids function).
            output_directory : str or Path : Directory to save the filtered output FASTA file.
        
        Returns:
            Path : Path to the filtered output FASTA file.
        """
        # Combine all sequence IDs from sequence_ids_dict into a single set
        combined_ids = set()
        for ids in sequence_ids_dict.values():
            combined_ids.update(ids)

        # Get the special IDs from special_dict
        special_ids = special_dict.get('special', set())

        # Prepare file paths
        fasta_file_path = pathlib.Path(fasta_file_path)
        output_path = pathlib.Path(output_directory) / f"{fasta_file_path.stem}_filtered.fasta"

        with fasta_file_path.open('r') as fasta, output_path.open('w') as output:
            write_sequence = False
            output_line = ""

            for line in fasta:
                if line.startswith('>'):  # Header line
                    # Extract the sequence ID (part after '>' and the last number after '_')
                    seqID_last_number = line.split('>')[1].split('_')[-1].strip()

                    # Check if the ID is in special_ids or combined_ids
                    prepend_cluster = seqID_last_number in special_ids
                    if seqID_last_number in combined_ids or prepend_cluster:
                        output_line = f">CLUSTER_{line[1:]}" if prepend_cluster else line
                        write_sequence = True
                    else:
                        write_sequence = False
                if write_sequence:
                    # Write the header or sequence to the output file
                    output.write(output_line if line.startswith('>') else line)

        print(f"Filtered FASTA file has been written to {output_path}")
        return output_path
    @staticmethod
    def create_fasta_file_with_full_TRS_sequences(df, filtered_df, results_directory, incremented_indices,state):
        """
        Creates a FASTA file containing full sequences using the filtered DataFrame and the original one.
        Assigns new unique IDs to each sequence that contain information about flanking TRS regions.
        
        Args:
            df (DataFrame): The original DataFrame with all sequences.
            filtered_df (DataFrame): The filtered DataFrame with sequences to include in the FASTA file.
            results_directory (str): Directory to save the resulting FASTA file.
            incremented_indices (list): A list of incremented indices used for assigning sequence identifiers.
        
        Returns:
            str: Path to the resulting FASTA file.
        """
        # Fetch unique genome IDs from the DataFrame
        ncbi_ids = df["GENOME"].unique().tolist()
        
        # Map genome IDs to organism names
        organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids)
        
        # Add a new column with taxonomic names
        filtered_df['Taxonomic Name'] = filtered_df["GENOME"].map(organism_map)
        
        # Generate the L_R_id column for unique IDs containing information about flanking TRS
        filtered_df['L_R_id'] = (
            filtered_df['Taxonomic Name'] + '_L' + filtered_df['L-No'].astype(str) + '_R' + filtered_df['R-No'].astype(str)
        )
        
        # Create the directory if it does not exist
        if state == 1:
            full_sequences_path = os.path.join(results_directory, "intermediate_results")
            os.makedirs(full_sequences_path, exist_ok=True)
        elif state == 2:
            full_sequences_path = os.path.join(results_directory, "final_results")
            os.makedirs(full_sequences_path, exist_ok=True)
        
        # Define the output FASTA file path
        if state == 1:
            full_sequences_fasta = os.path.join(full_sequences_path, "full_sequences.fasta")
        if state == 2:
            full_sequences_fasta = os.path.join(full_sequences_path, "full_sequences_final.fasta")
        
        # Write the sequences to the FASTA file
        with open(full_sequences_fasta, 'w') as file:
            for i, (index, row) in enumerate(filtered_df.iterrows()):
                # Write the sequence header with the L_R_id and incremented index
                file.write(f'>{row["L_R_id"]}_{incremented_indices[i]}\n')
                
                # Write the actual sequence
                file.write(f'{row[">SEQ"]}\n')
        
        print(f"Full sequences FASTA file has been created: {full_sequences_fasta}")
        
        return full_sequences_fasta
    @staticmethod
    def parse_fasta_indices(fasta_file):
        """
        Parse a FASTA file and extract the part after the last '_' from each sequence ID.
        Increment those values by 1 and return both the original and incremented indices as unique, sorted lists.

        These incremented values represent the real positions of the sequences in the CSV file containing all sequences.

        Args:
            fasta_file (str): Path to the input FASTA file.
        
        Returns:
            tuple: Two sorted lists, one of unique original indices and one of unique incremented indices.
        """
        incremented_indices = []
        original_indices = []

        with open(fasta_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    index = line.strip().split('_')[-1]
                    try:
                        original_index = int(index) - 1
                        original_indices.append(original_index)
                        incremented_index = original_index + 2
                        incremented_indices.append(incremented_index)
                    except ValueError:
                        print(f"Skipping invalid sequence ID: {line.strip()}")

        return sorted(set(original_indices)), sorted(set(incremented_indices))
    @staticmethod
    def extract_rows_from_csv(df,row_indices):
        """Filters the DataFrame using indices corresponding to rows in dataframe and sorts them using those indexes"""
        filtered_df = df.iloc[row_indices].sort_index()
        return filtered_df
    @staticmethod
    def extract_full_TRS_sequences(df,fasta_file,results_directory,state):
        """Executes the below functions"""
        original_indices,incremented_indices = BLASTProcessor.parse_txt_indices(fasta_file,state)
        filtered_df = BLASTProcessor.extract_rows_from_csv(df=df,row_indices=original_indices)
        BLASTProcessor.create_fasta_file_with_full_TRS_sequences(df,filtered_df,results_directory,incremented_indices,state)

    @staticmethod
    def parse_txt_indices(txt_file, state):
        """
        Parse a text file and extract the part after the last '_' from each sequence name in the first column.
        Increment the indices by 1 if state is 1; otherwise, return them as is. Return unique, sorted lists of indices.

        Args:
            txt_file (str): Path to the input text file.
            state (int): If 1, increment the indices; if 2, do not increment.

        Returns:
            tuple: Two sorted lists, one of unique original indices and one of unique processed indices.
        """
        original_indices = []
        processed_indices = []

        with open(txt_file, 'r') as file:
            for line in file:
                # Split the line into columns and get the first column
                first_column = line.strip().split('\t')[0]
                # Extract the part after the last '_'
                index = first_column.split('_')[-1]
                try:
                    original_index = int(index) - 1
                    original_indices.append(original_index)

                    # Process index based on state
                    if state == 1:
                        processed_index = original_index + 2
                    elif state == 2:
                        processed_index = original_index
                    else:
                        raise ValueError("Invalid state value. State must be 1 or 2.")

                    processed_indices.append(processed_index)
                except ValueError:
                    print(f"Skipping invalid sequence ID: {first_column}")

        return sorted(set(original_indices)), sorted(set(processed_indices))

