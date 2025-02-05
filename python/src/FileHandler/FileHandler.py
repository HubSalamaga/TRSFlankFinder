import os 
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import fnmatch
import sys
from src.SequenceProcessor import SequenceProcessor

class FileHandler:
    
    @staticmethod
    def ensure_directory_exists(directory_path):
        '''
        Ensures that the specified directory exists if it does not creates it. 
        '''
        try:
            if not os.path.exists(directory_path):
                os.makedirs(directory_path)
            else:
                print(f"Directory {directory_path} already exists")
        except PermissionError:
            print(f"Premission denied: Unable to create {directory_path}")
        except FileNotFoundError:
            print(f"One or more intermediate directories in the path do not exists")
        except Exception as e:
            print(f"Unidentified error has occured while creating {directory_path}: {e}")
    
    @staticmethod
    def convert_to_txt(directory_path):
        '''
        Identifies non-txt files in the directory and converts them to tab-separated .txt files
        '''
        try:
            non_txt_files = [f for f in os.listdir(directory_path) if not f.endswith(".txt")]
        except FileNotFoundError:
            print("The specified directory does not exist.")
            return
        except PermissionError:
            print(f"Permission denied to {directory_path}.")
            return
        except Exception as e:
            print(f"An unexpected error has occurred: {e}")
        
        if non_txt_files:
            print(f"No .txt files found in {directory_path}")
            for file in non_txt_files:
                print(f"Found file: {file} lacks extension")

            print("Directory should contain only blast results!")
            for file in non_txt_files:
                original_path = os.path.join(directory_path,file)
                new_path = f"{original_path}.txt"
                try:
                    os.rename(original_path,new_path)
                    print(f"Renamed {file} to {file}.txt")
                except FileNotFoundError:
                    print(f"Failed to rename {file}: File not found.")
                except PermissionError:
                    print(f"Failed to rename {file}: Permission denied.")
                except Exception as e:
                    print(f"Failed to rename {file}: {e}")
        else:
            print("No non-txt files were found")

    @staticmethod
    def filter_and_overwrite_files_in_directory(directory_path, file_pattern= "*.txt"):
        """
        Filters each file in a specified directory based on the last number in the identifiers of its lines,
        then overwrites each file with its filtered content for entries with paired identifiers.
        Writes entries with unique identifiers to a separate file named '<original_filename>_singles.txt'.
        Targets files matching a given file pattern (default is "*.txt").
        """       
        for filename in os.listdir(directory_path):
            if fnmatch.fnmatch(filename,file_pattern):
                file_path = os.path.join(directory_path,filename)
                single_entries_path = os.path.join(directory_path,f"{os.path.splitext(filename)[0]}_singles.txt")

                with open(file_path,'r') as file:
                    lines = file.readlines()

                number_count = {}
                filtered_entries = []
                non_unique_entries = []

                #First pass to count the occurences of each identifier
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        number_count[number] = number_count.get(number,0) + 1
                

                #Second pass to separate unique and non-unique identifiers
                for line in lines:
                    identifier = line.strip().split('\t')[0]
                    number = identifier.split('_')[-1]
                    if number.isdigit():
                        if number_count[number] > 1:
                            filtered_entries.append(line.strip())
                        else:
                            non_unique_entries.append(line.strip())
                
                # Overwrite the original file with paired identifiers
                with open(file_path, 'w') as file:
                    for entry in filtered_entries:
                        file.write(entry + '\n')

                # Write out single entries to a separate file
                with open(single_entries_path, 'w') as file:
                    for entry in non_unique_entries:
                        file.write(entry + '\n')
    
    @staticmethod
    def find_directory_by_name(directory_name, search_paths = None,auto=False):
        """
        Searches for directories with a specified name in provided locations and allows the user to select one if multiple are found.

        Params: 
        - directory_name (str): The name of the directory to search for.
        - search_paths (list of str): Optional. List of paths to search in. Defaults to common locations.
        - auto (bool): If True, automatically select the first found directory.
        
        Returns:
        - str : The selected path to the directory.
        """

        if search_paths is None:
            user_dir = os.path.expanduser('~')
            search_paths = [
                user_dir,
                os.path.join(user_dir,'Documents'),
                os.path.join(user_dir,'Desktop')
            ]
        
        found_directories = []

        for location in search_paths:
            try:
                for root,dirs, _ in os.walk(location):
                    if directory_name in dirs:
                        found_directories.append(os.path.join(root,directory_name))
            except PermissionError as e:
                print(f"Permission denied when accessing {location}. Error: {e}")
            except Exception as e:
                print(f"An unexpected error has occurred: {e}")

        if not found_directories:
            print("No directories found check your input.")
            return None
        
        if len(found_directories) == 1 or auto:
            return found_directories[0]
        
        print(f"Multiple {directory_name} directories found. Please select one:")
        for i, directory in enumerate(found_directories,start=1):
            print(f"{i}. {directory}")

        while True:
            selection = input("Enter the number of the directory you want to use(or type 'exit' to exit): ")
            if selection.lower() == 'exit':
                print('Exiting')
                return None
            
            try:
                selected_index = int(selection) - 1
                if 0 <= selected_index < len(found_directories):
                    return found_directories[selected_index]
                else:
                    print("Invalid selection. Try again")
            except ValueError:
                print("Invalid input. Please enter a number.")
    
    @staticmethod
    def read_fasta_ids(filename):
        """
        Read and return a set of FASTA IDs for a given file.
        """
        fasta_ids = set()

        def process_line(line):
            """
            Process a single line from the FASTA file, adding the ID to the set if the line is a header.
            """
            if line.startswith('>'):
                fasta_id = line.strip()[1:]
                fasta_ids.add(fasta_id)

        try:
            with open(filename, 'r') as f:
                for line in f:
                    process_line(line)
        except FileNotFoundError:
            print(f"Error: File {filename} not found. Do not move files during execution!")
        except PermissionError:
            print(f"Error: Permission denied when accessing {filename}")
        print("Fasta IDs processed successfully")
        return fasta_ids
    @staticmethod
    def parse_fasta_indices(fasta_file):
        """
        Parse fasta file and extract the part after the last '_' from each sequence id. 
        Increment those values by 1 and return them.
        """
        incremented_indices = []
        indices = []
        with open(fasta_file, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    index = line.strip().split('_')[-1]
                    indices.append(index)
                    incremented_index = int(index) - 1
                    incremented_indices.append(incremented_index)
        incremented_indices = list(set(incremented_indices))
        return incremented_indices
    @staticmethod
    def extract_rows_from_csv(df,row_indices):
        """
        Filters the DataFrame using indices corresponding to rows in dataframe and sorts them using those indexes.
        """
        filtered_df = df.iloc[row_indices].sort_index()
        return filtered_df
    @staticmethod
    def create_fasta_file_with_full_TRS_sequences(df,filtered_df, results_directory,incrementent_indices):
        """
        Creates fasta file containing full sequences using filtered dataframe and the original one.
        """
        ncbi_ids = df["GENOME"].unique().tolist()
        organism_map = SequenceProcessor.fetch_organism_names(ncbi_ids)
        filtered_df['Taxonomic Name'] = df["GENOME"].map(organism_map)
        indices = incrementent_indices
        filtered_df["L_R_id"] = filtered_df['Taxonomic Name'] + '_L' + filtered_df['L-No'].astype(str) + '_R' + filtered_df['R-No'].astype(str)
        full_sequences_path = os.path.join(results_directory, "final_output")
        full_sequences_fasta = os.path.join(full_sequences_path,"full_sequences.fasta")
        with open(full_sequences_fasta, 'w') as file:
            for index, row in filtered_df.iterrows():
                file.write(f'>{row["L_R_id"]}_{index}\n')
                file.write(f'{row[">SEQ"]}\n')
    @staticmethod
    def filter_and_overwrite_blast_file(directory_path):
        """
        Processes BLAST output files in the given directory to remove duplicate pairs of sequence IDs and accession numbers,
        keeping the first occurrence of each unique pair along with the entire line of data.
        """
        for filename in os.listdir(directory_path):
            if filename.endswith(".txt"):  # Process only .txt files
                file_path = os.path.join(directory_path, filename)
                unique_pairs = set()
                lines_to_keep = []

                try:
                    with open(file_path, 'r') as file:
                        for line in file:
                            try:
                                columns = line.strip().split('\t')
                                if len(columns) >= 2:
                                    pair = (columns[0], columns[1])
                                    if pair not in unique_pairs:
                                        unique_pairs.add(pair)
                                        lines_to_keep.append(line)
                            except Exception as e:
                                print(f"An unexpected error occurred while processing the line in {file_path}: {e}")
                    
                    # Write the filtered lines back to the file or a new file
                    with open(file_path, 'w') as file:
                        file.writelines(lines_to_keep)

                    print(f"Filtered unique sequence ID - accession pairs in {file_path}.")
                
                except PermissionError as e:
                    print(f"Permission denied while trying to access {file_path}: {e}")
                except IOError as e:
                    print(f"An error occurred while opening or reading {file_path}: {e}")
    @staticmethod
    def search_for_files(search_path,file_name):
        found = []
        for root, _, files in os.walk(search_path):
            for file in files:
                if fnmatch.fnmatch(file,file_name):
                    found.append(os.path.join(root,file))
        return found
    @staticmethod
    def find_file_by_name(file_name, folder):
        """
        Searches for files with a specified name pattern within a given folder.

        Args:
            file_name (str): The name of the file to search for, supports wildcards.
            folder (str): The directory to search within.

        Returns:
            list: A list of file paths that match the search pattern, or an empty list if no files are found.
        """
        if not file_name or not isinstance(file_name, str):
            print("Invalid or empty file name provided")
            return []
        
        if not folder or not os.path.isdir(folder):
            print("Invalid or empty folder provided")
            return []
        
        found_files = []
        for root, dirs, files in os.walk(folder):
            for name in fnmatch.filter(files, file_name):
                found_files.append(os.path.join(root, name))
        
        return found_files
    @staticmethod
    def process_files_with_filter(directory,filtered_dict):
        output_directory = os.path.join(directory,"unique_sequences")
        FileHandler.ensure_directory_exists(output_directory)

        for filename in os.listdir(directory):
            if filename.endswith('.txt'):
                input_filepath = os.path.join(directory,filename)
                output_filepath = os.path.join(output_directory,f"unique_{filename}")

                df = pd.read_csv(input_filepath,sep='\t',header=None)

                #Create a boolean mask for rows where the first column matches a key in filtered_dict
                mask = df[0].isin(filtered_dict.keys())

                #Instead of creating a filtered dataframe we will directly alter the last column
                df.loc[mask,df.columns[-1]] = df.loc[mask,0].map(filtered_dict)

                #Now, filter the DataFrame to keep only the rows that match the mask
                filtered_df = df[mask]
                unique_df = filtered_df.drop_duplicates(subset=[0])
                unique_df.to_csv(output_filepath, sep='\t', index= False, header=False)
                print(f"Final output saved at {output_filepath}")
    @staticmethod
    def check_job_status(job_id):
        """Check the job status using squeue."""
        try:
            result = subprocess.run(['squeue', '--job', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = result.stdout.decode().strip()
            return output
        except subprocess.CalledProcessError as e:
            print(f"Failed to query job status: {e.stderr.decode()}")
            return None
    @staticmethod
    def extract_file_names(file_paths):
        """
        Extracts file names from a list of file paths or a single file path.

        Parameters:
            file_paths(str or list): A single file path or a list of them.

        Returns:
            list: List of file names if input is a list or a single file name when it's a single file
        """
        if isinstance(file_paths,list):
            return [os.path.basename(file_path) for file_path in file_paths]
        elif isinstance(file_paths,str):
            return os.path.basename(file_paths)
        else:
            raise ValueError("Input should be either a string or a list of directory paths")
   
    @staticmethod
    def find_largest_cluster_size(cdhit_clusters_file):
        """
        Finds the size of the largest cluster in the CD-HIT `.clstr` file.

        Parameters:
            cdhit_clusters_file (str): Path to the CD-HIT `.clstr` file.

        Returns:
            int: Size of the largest cluster.
        """
        largest_cluster_size = 0
        current_cluster_size = 0

        try:
            with open(cdhit_clusters_file, "r") as file:
                for line in file:
                    if line.startswith('>Cluster'):
                        if current_cluster_size > largest_cluster_size:
                            largest_cluster_size = current_cluster_size
                        current_cluster_size = 0
                    else:
                        current_cluster_size += 1

                # Handle the last cluster in the file
                if current_cluster_size > largest_cluster_size:
                    largest_cluster_size = current_cluster_size

        except FileNotFoundError:
            print(f"CD-HIT cluster file {cdhit_clusters_file} not found.")
        except Exception as e:
            print(f"Error finding largest cluster size: {e}")

        return largest_cluster_size
    
    @staticmethod
    def get_sequence_ids_from_clusters(cdhit_clusters_file, cluster_size):
        """
        Extracts sequence IDs from clusters of a specific size in the CD-HIT `.clstr` file.

        Parameters:
            cdhit_clusters_file (str): Path to the CD-HIT `.clstr` file.
            cluster_size (int): Size of the cluster to extract IDs for.

        Returns:
            list: List of sequence IDs in clusters of the specified size.
        """
        sequence_ids = []
        current_cluster = []
        current_cluster_size = 0

        try:
            with open(cdhit_clusters_file, "r") as file:
                for line in file:
                    if line.startswith('>Cluster'):
                        # Process the previous cluster
                        if current_cluster_size == cluster_size:
                            sequence_ids.extend(current_cluster)
                        # Reset for the next cluster
                        current_cluster = []
                        current_cluster_size = 0
                    elif line.strip():
                        # Increment size and add sequence ID
                        current_cluster_size += 1
                        seq_id = line.split(">")[1].split("...")[0]  # Extract sequence ID
                        current_cluster.append(seq_id)

                # Handle the last cluster
                if current_cluster_size == cluster_size:
                    sequence_ids.extend(current_cluster)

        except FileNotFoundError:
            print(f"CD-HIT cluster file {cdhit_clusters_file} not found.")
        except Exception as e:
            print(f"Error processing clusters: {e}")

        return sequence_ids
    @staticmethod
    def filter_fasta_by_exact_cluster_size(input_fasta, output_fasta, cluster_ids):
        """
        Filters a FASTA file to include only sequences belonging to a specific cluster size.

        Parameters:
        - input_fasta (str): Path to the input FASTA file.
        - output_fasta (str): Path where the filtered FASTA file will be saved.
        - cluster_ids (list): List of FASTA IDs belonging to the specific cluster size.
        """
        try:
            ids_to_include = set(cluster_ids)  # Convert list to set for faster lookups
            write_sequence = False  # Flag to track whether to write the current sequence

            with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
                for line in infile:
                    if line.startswith(">"):  # Header line
                        current_id = line.strip()[1:]  # Extract sequence ID
                        write_sequence = current_id in ids_to_include
                    if write_sequence:
                        outfile.write(line)

        except FileNotFoundError:
            print(f"Error: File {input_fasta} not found.")
        except PermissionError:
            print(f"Error: Permission denied when accessing {input_fasta} or {output_fasta}.")
        except Exception as e:
            print(f"An error occurred: {e}")

    @staticmethod
    def extract_clusters_by_size(cdhit_clusters_file, cluster_size):
        """
        Extracts all clusters of a given size from the CD-HIT `.clstr` file.

        Parameters:
            cdhit_clusters_file (str): Path to the CD-HIT `.clstr` file.
            cluster_size (int): Size of the clusters to extract.

        Returns:
            list: A list of lists, where each sublist contains the sequence IDs of a cluster.
        """
        clusters = []
        current_cluster = []
        current_cluster_size = 0

        try:
            with open(cdhit_clusters_file, "r") as file:
                for line in file:
                    if line.startswith('>Cluster'):
                        # Process the previous cluster
                        if current_cluster_size == cluster_size:
                            clusters.append(current_cluster)
                        # Reset for the next cluster
                        current_cluster = []
                        current_cluster_size = 0
                    elif line.strip():
                        current_cluster_size += 1
                        seq_id = line.split(">")[1].split("...")[0]  # Extract sequence ID
                        current_cluster.append(seq_id)

                # Handle the last cluster
                if current_cluster_size == cluster_size:
                    clusters.append(current_cluster)

        except FileNotFoundError:
            print(f"CD-HIT cluster file {cdhit_clusters_file} not found.")
        except Exception as e:
            print(f"Error processing clusters: {e}")

        return clusters

    @staticmethod
    def write_fasta_for_individual_clusters(cdhit_clusters_file, input_fasta, cluster_size, output_folder):
        """
        Creates individual FASTA files for each cluster in a given cluster size.

        Parameters:
            cdhit_clusters_file (str): Path to the CD-HIT `.clstr` file.
            input_fasta (str): Path to the input FASTA file containing sequences.
            cluster_size (int): Cluster size to process.
            output_folder (str): Folder to save the FASTA files for individual clusters.
        """
        # Extract all clusters of the specified size
        clusters = FileHandler.extract_clusters_by_size(cdhit_clusters_file, cluster_size)

        for cluster_idx, cluster_ids in enumerate(clusters, start=1):
            # Define the output file for the current cluster
            cluster_fasta = os.path.join(output_folder, f"cluster_{cluster_idx}.fasta")

            # Write the sequences for the current cluster
            FileHandler.filter_fasta_by_exact_cluster_size(input_fasta, cluster_fasta, cluster_ids)
            print(f"Cluster {cluster_idx} in size {cluster_size}: {len(cluster_ids)} sequences written to {cluster_fasta}.")
    
    @staticmethod
    def create_fasta_for_all_clusters(cdhit_clusters_file, input_fasta, output_folder, create_individual_files=True):
        """
        Creates FASTA files for each cluster size in the CD-HIT `.clstr` file, skipping empty clusters.
        Each file contains only sequences from clusters of the exact specified size.
        Optionally, creates individual FASTA files for each cluster in sizes greater than 1.

        Parameters:
            cdhit_clusters_file (str): Path to the CD-HIT `.clstr` file.
            input_fasta (str): Path to the input FASTA file containing sequences.
            output_folder (str): Folder to save the FASTA files for each cluster size.
            create_individual_files (bool): Whether to create individual FASTA files for each cluster.
        """
        # Ensure the output folder exists
        FileHandler.ensure_directory_exists(output_folder)

        # Find the largest cluster size
        largest_cluster_size = FileHandler.find_largest_cluster_size(cdhit_clusters_file)
        print(f"Largest cluster size: {largest_cluster_size}")

        for cluster_size in range(1, largest_cluster_size + 1):
            # Get sequence IDs only from clusters of the exact size
            sequence_ids = FileHandler.get_sequence_ids_from_clusters(cdhit_clusters_file, cluster_size)

            # Skip if no sequences found for the cluster size
            if not sequence_ids:
                print(f"Cluster size {cluster_size}: No sequences found. Skipping.")
                continue

            # Define the output folder for the current cluster size
            cluster_size_folder = os.path.join(output_folder, f"cluster_size_{cluster_size}")
            FileHandler.ensure_directory_exists(cluster_size_folder)

            # Write sequences belonging to the current cluster size into a single FASTA file
            output_fasta = os.path.join(output_folder, f"cluster_size_{cluster_size}.fasta")
            FileHandler.filter_fasta_by_exact_cluster_size(input_fasta, output_fasta, sequence_ids)
            print(f"Cluster size {cluster_size}: {len(sequence_ids)} sequences written to {output_fasta}.")

            # Skip creating individual cluster FASTA files if disabled or for cluster size 1
            if not create_individual_files or cluster_size == 1:
                if cluster_size == 1:
                    print(f"Skipping individual FASTA file creation for clusters of size 1.")
                continue
            # Separate sequences for each cluster within the cluster size 
            '''Maybe something to read what clusters were already calculated because for larger datasets it's very time consuming mayber something like
            
            if not create_individual_files or cluster_size in range(1, clusters_to_skip):
                if cluster_size in range(1,clusters_to_skip):
                print(f"Skipping individual FASTA file creation for clusters of cluster_size)
            '''
            FileHandler.write_fasta_for_individual_clusters(cdhit_clusters_file, input_fasta, cluster_size, cluster_size_folder)
