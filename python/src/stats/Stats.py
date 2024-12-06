import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
import glob
from scipy.stats import binom_test

class Stats:
    def __init__(self):
        self.file_info = {}  # Stores file details: {file_path: {'type': 'txt'/'fasta', 'lines'/ 'sequences': count, 'folder': folder_name}}
        self.lr_counts = {}  # Stores L/R counts for files: {file_path: {'Lxx': count, 'Rxx': count}}
    def count_L_R(self, file_path):
        """
        Reads a .txt file to count lines and analyze L(number)/R(number) patterns.
        Parameters:
            file_path (str): Path to the .txt file.
        """
        line_count = 0
        lr_pattern_counts = {}
        folder_name = os.path.dirname(file_path)
        with open(file_path, 'r') as file:
            for line in file:
                line_count += 1
                first_column = line.split()[0] if line.split() else ""
                for part in first_column.split('_'):
                    if part.startswith(('L', 'R')) and part[1:].isdigit():
                        lr_pattern_counts[part] = lr_pattern_counts.get(part, 0) + 1

        self.file_info[file_path] = {'type': 'txt', 'lines': line_count, 'folder': folder_name}
        self.lr_counts[file_path] = lr_pattern_counts
    def count_sequences_and_l_r_patterns(self, file_path):
        """
        Reads a .fasta file to count the number of sequences and analyze L(number)/R(number) patterns in seqID.
        Parameters:
            file_path (str): Path to the .fasta file.
        """
        sequence_count = 0
        lr_pattern_counts = {}
        folder_name = os.path.dirname(file_path)
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    sequence_count += 1
                    seqID = line[1:].strip()  # Extract seqID, removing '>' and whitespace
                    for part in seqID.split('_'):
                        if part.startswith(('L', 'R')) and part[1:].isdigit():
                            lr_pattern_counts[part] = lr_pattern_counts.get(part, 0) + 1

        self.file_info[file_path] = {'type': 'fasta', 'sequences': sequence_count, 'folder': folder_name}
        self.lr_counts[file_path] = lr_pattern_counts
    def get_lr_counts(self, file_path):
        """
        Returns the L(number)/R(number) pattern counts for a file.
        Parameters:
            file_path (str): Path to the file.
        """
        lr_counts = self.lr_counts.get(file_path,{})
        sorted_lr_counts = {k: lr_counts[k] for k in sorted(lr_counts.keys())}
        return sorted_lr_counts
    def plot_lr_counts(self,file_paths,file_titles,results_directory):
        """
        Plots L/R counts for each file and saves the plots to specified results folder
        """
        for file_path in file_paths:
            lr_counts = self.get_lr_counts(file_path)
            max_count = max(lr_counts.values()) if lr_counts else 0

            file_name = os.path.basename(file_path)
            if max_count == 0:
                print(f"No L/R counts found for file: {file_name}")
            
            #Get the title for the current filepath from title dictionary
            title = file_titles.get(file_path,os.path.basename(file_path))

            #Separate counts
            l_counts = {k: v for k, v in lr_counts.items() if 'L' in k}
            r_counts = {k: v for k, v in lr_counts.items() if 'R' in k}

            #Total L 
            total_l_sequences = sum(l_counts.values())
            #Total R
            total_r_sequences = sum(r_counts.values())
            #Total L+R 
            total_sequences = sum(lr_counts.values())

            #Calculate average counts
            avg_l_count = np.mean(list(l_counts.values())) if l_counts else 0
            avg_r_count = np.mean(list(r_counts.values())) if r_counts else 0

            fig,axs = plt.subplots(1,2, figsize=(20,6)) #adjust as needed


            #L counts plot
            sorted_l_keys = sorted(l_counts.keys(), key= lambda x: int(re.search(r'\d+', x).group()))
            axs[0].bar(sorted_l_keys, [l_counts[key] for key in sorted_l_keys], color='green')
            axs[0].axhline(y=avg_l_count, linestyle='--', color='red', label=f'Average: {avg_l_count:.2f}')
            axs[0].set_title(f'{title} - L Counts ({total_l_sequences} L, Total: {total_sequences})')
            axs[0].set_xlabel('L Combinations')
            axs[0].set_ylabel('Counts')
            axs[0].set_ylim(0, max_count + (max_count / 10))  # Set upper bound for Y-axis based on file
            axs[0].tick_params(axis='x', rotation=45)
            axs[0].set_xticks(range(len(sorted_l_keys)))  # Set the tick positions
            axs[0].set_xticklabels(sorted_l_keys, fontsize=7)  # Set the tick labels
            #R counts plot
            sorted_r_keys = sorted(r_counts.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
            axs[1].bar(sorted_r_keys, [r_counts[key] for key in sorted_r_keys], color='blue')
            axs[1].axhline(y=avg_r_count, linestyle='--', color='red', label=f'Average: {avg_r_count:.2f}')
            axs[1].set_title(f'{title} - R Counts ({total_r_sequences} R, Total: {total_sequences})')
            axs[1].set_xlabel('R Combinations')
            axs[1].set_ylabel('Counts')
            axs[1].set_ylim(0, max_count + (max_count / 10))  # Set upper bound for Y-axis based on file
            axs[1].tick_params(axis='x', rotation=45)
            axs[1].set_xticks(range(len(sorted_r_keys)))  # Set the tick positions
            axs[1].set_xticklabels(sorted_r_keys, fontsize=7)  # Set the tick labels
            #Legends
            axs[0].legend()
            axs[1].legend()

            # Adjust layout
            plt.tight_layout()

            save_path = os.path.join(results_directory,"summary")
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            save_path = os.path.join(save_path,f"{file_name}_summary.png")
            plt.savefig(save_path, bbox_inches='tight', facecolor="white", dpi=300)
    @staticmethod
    def create_trs_class_dataframe(cluster_folder):
        """
        Generates a DataFrame with TRS class counts across cluster sizes.

        Parameters:
            cluster_folder (str): Path to the folder containing `cluster_size_*.fasta` files.

        Returns:
            pd.DataFrame: A DataFrame with cluster sizes as columns and TRS classes as rows.
        """


        # Dictionary to store TRS counts for each cluster size
        trs_class_counts = {}

        # Get all `cluster_size_*.fasta` files
        fasta_files = glob.glob(os.path.join(cluster_folder, "cluster_size_*.fasta"))

        for fasta_file in fasta_files:
            # Extract cluster size from file name
            cluster_size = int(os.path.basename(fasta_file).split("_")[2].split(".")[0])

            # Initialize cluster size entry if not present
            if cluster_size not in trs_class_counts:
                trs_class_counts[cluster_size] = {}

            # Process the FASTA file
            with open(fasta_file, "r") as f:
                for line in f:
                    if line.startswith(">"):  # Sequence header
                        seq_name = line.strip()[1:]  # Remove '>'
                        trs_class = seq_name.split("_")[-2]  # Extract only the TRS class prefix
                        trs_class_counts[cluster_size][trs_class] = trs_class_counts[cluster_size].get(trs_class, 0) + 1

        # Build the DataFrame from the dictionary
        # Extract all unique TRS classes
        trs_classes = set()
        for class_counts in trs_class_counts.values():
            trs_classes.update(class_counts.keys())

        # Custom sorting for TRS classes (e.g., L1, L2, ..., L60, R1, R2, ..., R60)
        def natural_sort_key(trs_class):
            prefix, number = trs_class[0], int(trs_class[1:])
            return prefix, number

        sorted_trs_classes = sorted(trs_classes, key=natural_sort_key)

        # Create a DataFrame
        df = pd.DataFrame(columns=sorted(trs_class_counts.keys()), index=sorted_trs_classes).fillna(0)

        # Populate the DataFrame
        for cluster_size, class_counts in trs_class_counts.items():
            for trs_class, count in class_counts.items():
                df.at[trs_class, cluster_size] = count

        # Rename columns for better readability
        df.columns = [f"Cluster Size {size}" for size in df.columns]
        df.index.name = "TRS Class"
        return df
    @staticmethod
    def get_trs_class_totals(df):
        """
        Creates a new DataFrame with total counts of each TRS class across the dataset.

        Parameters:
            df (pd.DataFrame): Original DataFrame with TRS classes as rows and cluster sizes as columns.

        Returns:
            pd.DataFrame: A new DataFrame with TRS classes and their total counts.
        """
        # Calculate the total count for each TRS class
        total_counts = df.sum(axis=1)

        # Create a new DataFrame
        total_counts_df = pd.DataFrame(total_counts, columns=["Total Count"])
        total_counts_df.index.name = "TRS Class"

        # Sort by total counts in descending order (optional)
        total_counts_df = total_counts_df.sort_values(by="Total Count", ascending=False)

        return total_counts_df
    def create_lr_counts_dataframe(self, file_path):
        """
        Reads a file, counts L and R occurrences, and creates a DataFrame of TRS class counts.

        Parameters:
            file_path (str): Path to the file containing TRS data.

        Returns:
            pd.DataFrame: A DataFrame with TRS classes and their counts, sorted in descending order.
        """
        # Perform the counting and get LR counts
        testing_fasta_counting = self.count_L_R(file_path)
        sorted_lr_counts = self.get_lr_counts(file_path)
        
        # Convert to DataFrame and process
        sorted_lr_counts = pd.DataFrame([sorted_lr_counts]).transpose()
        sorted_lr_counts = sorted_lr_counts.rename(columns={0: "Counts"})
        sorted_lr_counts = sorted_lr_counts.sort_values(by="Counts", ascending=False)
        sorted_lr_counts.index.name = "TRS class"
        
        return sorted_lr_counts
    @staticmethod
    def test_binomial_for_classes(df):
        """
        Performs a Binomial Test for each class to determine if observed counts
        significantly differ from expected proportions.

        Parameters:
            df (pd.DataFrame): DataFrame with 'Counts', 'Total Count', 'Proportion_Found'.

        Returns:
            pd.DataFrame: A DataFrame with observed counts, expected proportions,
                        and p-values for each class.
        """
        # Total number of sequences that passed the filters
        total_passed = df["Counts"].sum()/2

        # Initialize lists to store results
        p_values = []
        # Calculate expected counts for each class
        df["Expected Count"] = df["Proportion_Found"] * total_passed

        # Perform Binomial Test for each class
        for _, row in df.iterrows():
            observed = row["Counts"] if not pd.isna(row["Counts"]) else 0
            expected_prob = row["Proportion_Found"]  # Proportion of this class in the dataset

            # Perform binomial test
            if total_passed > 0 and expected_prob > 0:  # Ensure valid test conditions
                p_val = binom_test(observed, n=total_passed, p=expected_prob, alternative="two-sided")
            else:
                p_val = float("nan")  # Assign NaN if test conditions are invalid
            p_values.append(p_val)

        # Add p-values to DataFrame
        df["P-Value"] = p_values

        #Rename and reformat columns
        df = df.rename(columns={
        "Proportion_Found": "% Found",
        "Proportion_Passed": "% Passed"
        })
        df["% Found"] = df["% Found"] * 100
        df["% Passed"] = df["% Passed"] * 100

        return df[["Counts", "Total Count", "% Found", "% Passed", "Expected Count", "P-Value"]]            
    @staticmethod
    def merge_trs_classes_results(df1,df2):
        merged_df = pd.merge(df1,df2, left_index=True,right_index=True,how="outer")
        total_sequences_found = merged_df['Total Count'].sum()/2
        total_sequences_passed = merged_df['Counts'].sum()
        merged_df["Proportion_Found"] = merged_df["Total Count"]/ total_sequences_found
        merged_df["Proportion_Passed"] = merged_df["Counts"]/ total_sequences_passed
        return merged_df
    @staticmethod
    def save_binom_results(df,folder_path):
        full_results = os.path.join(folder_path,"binomial_results_full.csv")
        df.to_csv(full_results,index=True)
        with_counts = df[(df["Counts"] >= 1) & (df["Expected Count"] <= df["Counts"])]
        results_with_counts = os.path.join(folder_path,"binomial_results_with_counts.csv")
        with_counts.to_csv(results_with_counts,index=True)
        significant = df[df["P-Value"] <= 0.05]
        results_significant = os.path.join(folder_path,"binomial_results_significant.csv")
        significant.to_csv(results_significant,index=True)
        return with_counts,significant
