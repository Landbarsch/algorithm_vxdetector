#!/usr/bin/python

import argparse
import glob
import os
import sys
import warnings
import pandas as pd
import Output_counter as Output_counter
import files_manager as files_manager
from interact_bowtie2 import mapbowtie2, buildbowtie2
from interact_bedtools import overlap
import random

def sample_fastq(file_path, sample_size, sampled_indices=None):
    '''Get random parts from the FASTQ file based on shared indices for paired-end reads.'''
    sampled_reads = []  # List to hold sampled reads
    with open(file_path, 'r') as f:
        lines = f.readlines()
        total_reads = len(lines) // 4  # FASTQ files store 4 lines per read
        if sampled_indices is None:
            # If no indices are provided, sample randomly from the total reads
            sampled_indices = sorted(random.sample(range(total_reads), sample_size))
        for idx in sampled_indices:
            # Get the read (4 lines per read) and append to the sampled_reads list
            read = lines[idx*4:(idx+1)*4]
            sampled_reads.extend(read)
    return sampled_reads, sampled_indices

def save_sampled_fastq(sampled_reads, output_path):
    '''Save sampled reads to a temporary FASTQ file.'''
    with open(output_path, 'w') as f:
        f.writelines(sampled_reads)

def do_statistic(result):
    '''Compute statistics (mean and standard deviation) from alignment results.'''
    average = result.mean(numeric_only=True).to_frame().T  # Compute the mean of numeric columns
    region = (result['Sequenced variable region'].mode().values)  # Get most frequent region
    region = ' / '.join(str(r) for r in region)  # Join the modes into a single string
    region = region.replace('\'', '').replace('[', '').replace(']', '')  # Clean up string format
    average['Sequenced variable region'] = region  # Add region info to the average dataframe
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'
    std_dev = result.std(numeric_only=True).to_frame().T  # Compute standard deviation of numeric columns
    statistic = pd.concat([average, std_dev], axis=0)  # Combine mean and standard deviation
    # Select specific columns for the final output
    statistic = statistic[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                           'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    statistic['row_descriptor'] = ['Average', 'Standard deviation']  # Add row descriptors
    statistic = statistic.set_index('row_descriptor')  # Set descriptor as index
    result = pd.concat([statistic, result], axis=0)  # Append stats to original results
    return result

def do_output(result, new_file, single_file):
    '''Process and output results to a CSV file.'''
    warnings.simplefilter(action='ignore', category=FutureWarning)
    result = pd.DataFrame(result).T.sort_index()  # Convert result to DataFrame and sort

    # Convert columns to numeric where applicable
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')
    
    if single_file is False:
        result = do_statistic(result)  # Add statistics if processing multiple files
    else:
        result = result[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                         'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    
    # Write the result DataFrame to a CSV file
    result.to_csv(new_file, index=True)

def workflow(file_dir, new_file, write_csv, sample_size):
    '''Main workflow for processing FASTQ files and performing sequence alignment.'''
    path = files_manager.get_lib()  # Get the main working directory
    temp_path = files_manager.tmp_dir(path, temp_path=None)  # Create a temporary directory
    paired = False  # To check if the input is paired-end
    single_file = False  # To check if processing a single file
    result = dict()  # Dictionary to store results from different files
    buildbowtie2(path)  # Build the Bowtie2 index

    # Check if the input path exists
    if not os.path.exists(file_dir):
        raise ValueError(f"The provided path {file_dir} does not exist.")
    
    # Check if there are FASTQ files in the directory
    if glob.glob(f'{file_dir}**/*.fastq*', recursive=True) == [] and os.path.isdir(file_dir):
        files_manager.tmp_dir(path, temp_path)
        raise ValueError('There were no FASTQ files in this directory')
    
    # Process single file
    if os.path.isfile(file_dir):
        single_file = True
        file_name = os.path.basename(file_dir)
        read2_file = os.path.join(os.path.dirname(file_dir), file_name.replace('_R1_', '_R2_'))  # Get the paired file if exists
        if '_R1_' in file_name and os.path.exists(read2_file):
            paired = True  # Mark as paired if both R1 and R2 files exist

        # Sample reads from R1 and R2 with same indices
        sampled_reads_R1, sampled_indices = sample_fastq(file_dir, sample_size)
        sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)

        # Save sampled reads to temporary files
        temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
        temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
        save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
        save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

        # Run Bowtie2 alignment
        aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2, path, temp_path, paired)

        if Error is True:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file does not look like a fastq file')
        if paired is True and Output_counter.rawincount(f'{temp_path}paired.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')
        
        # Perform overlap analysis using Bedtools
        overlap(path, temp_path, aligned_path)
        if paired is False and Output_counter.rawincount(f'{temp_path}BED.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')

        # Generate results from output
        file_name = file_name.rsplit('.f', 1)[0]
        file_name = file_name.replace('_R1_001', '')
        result[file_name] = Output_counter.create_row(temp_path, paired)
        
    # Process directory of files
    elif os.path.isdir(file_dir):
        single_file = False
        total_files = len(glob.glob(f'{file_dir}**/*.fastq*', recursive=True))  # Count number of FASTQ files
        processed_files = 0
        for fq_file in glob.glob(f'{file_dir}**/*.fastq*', recursive=True):
            processed_files += 1
            paired = False
            if '_R2_' in fq_file:
                continue  # Skip R2 files to avoid processing duplicates
            
            file_name = os.path.basename(fq_file)
            read2_file = os.path.join(os.path.dirname(fq_file), file_name.replace('_R1_', '_R2_'))
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True  # Mark as paired if both R1 and R2 files exist

            # Sample reads from R1 and R2 with same indices
            sampled_reads_R1, sampled_indices = sample_fastq(fq_file, sample_size)
            sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)
            
            # Save sampled reads to temporary files
            temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
            temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
            save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
            save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

            # Run Bowtie2 alignment
            aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2, path, temp_path, paired)
            if Error is True:
                continue  # Skip if there's an error in processing the file

            # Perform overlap analysis using Bedtools
            overlap(path, temp_path, aligned_path)

            # Generate results from output
            file_name = file_name.rsplit('.f', 1)[0]
            file_name = file_name.replace('_R1_001', '')
            result[file_name] = Output_counter.create_row(temp_path, paired)
                
    # Clean up temporary directory
    files_manager.tmp_dir(path, temp_path)
    
    # Output the result as CSV or display
    do_output(result, new_file, single_file)
    if write_csv is True:
        new_file = (f'{path}Output/{os.path.basename(os.path.dirname(file_dir))}.csv')
        do_output(result, new_file, single_file)

def main():
    '''Main function to parse user arguments and initiate the workflow.'''
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This program tries to find which variable region of the 16S sequence was sequenced'))
    parser.add_argument('dir_path', help=('Directory path of the directory containing multiple fastq or fasta files.'))
    parser.add_argument('-o', '--output', dest='output_file', default=sys.stdout, 
                        help='User can specify a file format in which the output is written in the Output folder.')
    parser.add_argument('-c', '--csv', dest='write_csv', action='store_true', 
                        help='If set the output will be written in a .csv file in the Output folder')
    parser.add_argument('-s', '--sample_size', dest='sample_size', type=int, default=1000, 
                        help='Number of reads to sample from each FASTQ file')
    
    # Parse the user input arguments
    args = parser.parse_args()
    
    # Start the workflow with provided arguments
    workflow(args.dir_path, args.output_file, args.write_csv, args.sample_size)

if __name__ == '__main__':
    main()  # Entry point of the script
