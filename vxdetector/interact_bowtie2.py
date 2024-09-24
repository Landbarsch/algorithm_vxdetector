#!/usr/bin/python

import os
import shutil
import subprocess

bowtie2_path = shutil.which('bowtie2')
if bowtie2_path is None:
    bowtie2_path = '$CONDA/bin/bowtie2'
samtools_path = shutil.which('samtools')
if samtools_path is None:
    samtools_path = '$CONDA/bin/samtools'
bedtools_path = shutil.which('bedtools')
if bedtools_path is None:
    bedtools_path = '$CONDA/bin/bedtools'
# finds where bowtie2, samtools and bedtools
# are installed


def buildbowtie2(path):
    r'''Builds bowtie2 index

    This function builds an index for bowtie2 it is the equivalent
    of building it manually with the "" command.

    Parameters
    ----------
    path : str
        The program filepath. The directory where it needs to look
        for the reference genome and where the index should be saved.

    '''
    input_ref = f'{path}Indexed_bt2/85_otus.fasta'
    # reference genome greengenes is used
    output_path = f'{path}Indexed_bt2/bowtie2'
    if os.path.exists(f'{output_path}.1.bt2'):
        # only builds an index for the alignment if there isn't already one
        pass
    else:
        os.system(f'{bowtie2_path}-build -f {input_ref} {output_path} '
                  '> /dev/null')
        # if a indexfile is missing a new index is build


def mapbowtie2(fasta_file, read2_file, path, temp_path, paired, bowtie2_params=None):    
    r'''Maps reads against index

    This function maps read files against a previously build index.
    The bowtie standard output .sam is piped into samtools view and
    converted to .bam. In this step all unmapped reads and those with a
    lower mapquality than requested are discarded.
    If the reads are paired the .bam file is converted to .bed with bedtools
    bamtobed that way they are properly paired and can be intersected.

    Parameters
    ----------
    fasta_file : str
        Path to a .fastq file (or .fastq.gz). This file contains the (forward)
        reads which are mapped against the reference.
    read2_file : str
        Path to a .fastq file (or .fastq.gz) containing backwards reads.
        Is disregarded if paired = False. If not needed can be ''.
    path : str
        The program filepath. The directory where it needs to look
        for the index files.
    temp_path : str
        Path to a (temporary) folder where the resulting .log files and output
        files are saved
    paired : Bool
        Dictates wether bowtie2 alignes paired or unpaired reads.
        If paired is set to False bowtie2 will do an unpaired alignment
        otherwise it will do a paired one.

    Returns
    -------
    aligned_path : str
        Path to the converted bowtie2 Output file in which either paired or
        unpaired reads have been mapped against the reference.
    Error : Bool
        Wether or not bowtie2 ended with an error.
        In the context of the whole program Error = True leads to either
        a raised Exception (if only a single file was given) or to the current
        file being skipped (if a directory was given).

    '''
    
# Define the path to the Bowtie2 index
    index_path = f'{path}Indexed_bt2/bowtie2'
    
    # Check if the index files exist
    if not os.path.exists(f'{index_path}.1.bt2'):
        raise FileNotFoundError(f'No index files found under "{index_path}"')

    # Define paths for log files
    log_path = f'{temp_path}bowtie2.log'
    bed_logpath = f'{temp_path}bed.log'
    Error = False
    
    # Create the base Bowtie2 command
    bowtie2_base_cmd = [bowtie2_path, '-x', index_path, '--fast']
    
    # Add additional Bowtie2 parameters if provided
    if bowtie2_params:
        bowtie2_base_cmd += bowtie2_params

    # If the reads are paired, use paired-end processing
    if paired:
        aligned_path = f'{temp_path}paired.bed'
        bowtie2_base_cmd += ['-1', fasta_file, '-2', read2_file]
        
        # Create the pipeline for Bowtie2 -> Samtools -> Bedtools
        cmd_bowtie = bowtie2_base_cmd
        cmd_samtools = [samtools_path, 'view', '-b', '-q', '30', '-S', '-F', '4']
        cmd_bedtools = [bedtools_path, 'bamtobed', '-bedpe', '-i', 'stdin']

        # Run the subprocess pipeline, writing output directly to files
        with open(log_path, 'w') as log_file, open(bed_logpath, 'w') as bed_log, open(aligned_path, 'w') as output_file:
            # Start Bowtie2 process
            bowtie_process = subprocess.Popen(cmd_bowtie, stderr=log_file, stdout=subprocess.PIPE)
            # Pass output of Bowtie2 to Samtools
            samtools_process = subprocess.Popen(cmd_samtools, stdin=bowtie_process.stdout, stdout=subprocess.PIPE)
            # Pass output of Samtools to Bedtools and write to file
            subprocess.run(cmd_bedtools, stdin=samtools_process.stdout, stdout=output_file, stderr=bed_log)
        
    # If the reads are unpaired, process them without Bedtools
    else:
        aligned_path = f'{temp_path}unpaired.bam'
        bowtie2_base_cmd += ['-U', fasta_file]
        
        # Create the pipeline for Bowtie2 -> Samtools
        cmd_bowtie = bowtie2_base_cmd
        cmd_samtools = [samtools_path, 'view', '-b', '-q', '30', '-S', '-F', '4', '-o', aligned_path]

        # Run the Bowtie2 and Samtools subprocess pipeline
        with open(log_path, 'w') as log_file:
            # Start Bowtie2 process and pipe output to Samtools
            bowtie_process = subprocess.Popen(cmd_bowtie, stderr=log_file, stdout=subprocess.PIPE)
            subprocess.run(cmd_samtools, stdin=bowtie_process.stdout)

    # Check for errors in the Bowtie2 log file
    with open(log_path, 'r') as log:
        if any("error" in line.lower() for line in log):
            Error = True

    return aligned_path, Error
