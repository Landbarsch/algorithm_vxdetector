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
import tempfile
from multiprocessing import Pool


def do_statistic(result):
    '''Performs statistical analysis on the DataFrame

    Calculates the mean and standard deviation for all numeric columns
    in the DataFrame, and determines the most common sequenced variable 
    region. Adds these statistics to the DataFrame and returns the updated 
    DataFrame.
    '''
    average = result.mean(numeric_only=True).to_frame().T  # Calculate mean for numeric columns
    region = (result['Sequenced variable region'].mode().values)  # Find the most common variable region
    region = ' / '.join(str(r) for r in region)  # Format the result for better readability
    region = region.replace('\'', '').replace('[', '').replace(']', '')  # Clean up formatting
    average['Sequenced variable region'] = region
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'  # Handle cases where reads are not paired
    std_dev = result.std(numeric_only=True).to_frame().T  # Calculate standard deviation for numeric columns
    statistic = pd.concat([average, std_dev], axis=0)  # Combine mean and standard deviation dataframes
    # Select relevant columns and reorder them
    statistic = statistic[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                           'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    statistic['row_descriptor'] = ['Average', 'Standard deviation']  # Add descriptors for the statistics
    statistic = statistic.set_index('row_descriptor')  # Set the row descriptors as the index
    result = pd.concat([statistic, result], axis=0)  # Combine the statistical results with the original data
    return result


def do_output(result, new_file, single_file):
    '''Writes the results into a CSV file

    Converts the dictionary of results into a DataFrame, calculates statistics 
    if working with multiple files, and writes the data into a CSV file.
    '''
    warnings.simplefilter(action='ignore', category=FutureWarning)  # Ignore warnings about future behavior
    result = pd.DataFrame(result).T.sort_index()  # Convert the dictionary to a DataFrame and sort by index
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')  # Convert columns to numeric where possible
    if single_file is False:
        result = do_statistic(result)  # If multiple files, calculate statistics
    else:
        # Select relevant columns for single file processing
        result = result[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                         'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    result.to_csv(new_file, index=True)  # Write the results to a CSV file


def process_file(fq_file, path, bowtie2_params):
    '''Processes a single fastq file for alignment and analysis
    
    Aligns the sequences using Bowtie2, analyzes overlaps using Bedtools, 
    and collects the output in a temporary directory. Handles paired-end 
    reads if applicable.
    '''
    paired = False  # Initialize as unpaired
    result = {}
    try:
        with tempfile.TemporaryDirectory() as temp_path:  # Create a temporary directory for intermediate files
            file_name = os.path.basename(fq_file)  # Get the base filename
            read2_file = os.path.join(os.path.dirname(fq_file),
                                      file_name.replace('_R1_', '_R2_'))  # Look for the corresponding paired file
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True  # Set as paired if the second read file exists
            # Call Bowtie2 to align the sequences
            aligned_path, Error = mapbowtie2(fq_file, read2_file, path, temp_path, paired)
            if Error:
                return None  # Skip the file if there is an error
            # Call Bedtools to analyze overlaps
            overlap(path, temp_path, aligned_path)
            # Prepare the result dictionary
            file_name = file_name.rsplit('.f', 1)[0]
            file_name = file_name.replace('_R1_001', '')
            result[file_name] = Output_counter.create_row(temp_path, paired)  # Collect output data
    except Exception as e:
        print(f"Error processing file {fq_file}: {e}")  # Print error message for failed processing
    return result


def workflow(file_dir, new_file, write_csv,bowtie2_params):
    '''Main workflow for handling file directories or single files
    
    This function builds the Bowtie2 index, checks if a single file or a directory
    is provided, and processes files using multiprocessing. Writes the output 
    to CSV if specified.
    '''
    path = files_manager.get_lib()  # Get the path to the reference library
    buildbowtie2(path)  # Build the Bowtie2 index

    result = {}
    single_file = False

    # Check if input is a single file or directory
    if os.path.isfile(file_dir):
        single_file = True  # Set to single file mode
        result = process_file(file_dir, path, bowtie2_params)  # Process the single file
    elif os.path.isdir(file_dir):
        # Get all fastq files in the directory (ignoring read 2 files for paired-end)
        fastq_files = glob.glob(f'{file_dir}/**/*.fastq*', recursive=True)
        # Process files in parallel using multiprocessing
        with Pool() as pool:
            results = pool.starmap(process_file, [(fq_file, path,bowtie2_params) for fq_file in fastq_files if '_R2_' not in fq_file])
        
        # Collect the results from the parallel processing
        for res in results:
            if res:
                result.update(res)
    
    do_output(result, new_file, single_file)  # Write the results to output file
    if write_csv:
        output_csv = (f'{path}Output/{os.path.basename(os.path.dirname(file_dir))}.csv')  # Prepare CSV path
        do_output(result, output_csv, single_file)  # Write CSV output if specified


def main():
    '''Main function for parsing user input and starting the workflow
    
    Uses argparse to collect command-line arguments for directory path, 
    output file, and CSV option, then passes them to the workflow function.
    '''
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This program tries to determine which variable region of the 16S '
        'sequence was sequenced.'))  
    parser.add_argument('dir_path',
                        help=('Directory path of the directory containing '
                              'multiple fastq or fasta files.'))  
    parser.add_argument('-o', '--output', dest='output_file',
                        default=sys.stdout,
                        help='User can specify a file format in which the \
                        output is written in the Output folder.')  
    parser.add_argument('-c', '--csv', dest='write_csv', action='store_true',
                        help='If set the output will be written in a \
                       .csv file in the Output folder')  
    parser.add_argument('-b','--bowtie2-params', dest='bowtie2_params', default="--threads 4,--fast",
                        help='Additional parameters to pass to Bowtie2')
    
    args = parser.parse_args()  
    workflow(args.dir_path, args.output_file, args.write_csv, args.bowtie2_params)  # Call the workflow with parsed arguments


if __name__ == '__main__':
    main()  
