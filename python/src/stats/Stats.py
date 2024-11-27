import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re

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

            
