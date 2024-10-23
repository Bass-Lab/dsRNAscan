#!/usr/bin/env python3
import os
import locale
import glob
from Bio import SeqIO
from Bio.Emboss.Applications import EInvertedCommandline
import argparse
import subprocess
import re
import RNA
import sys
import numpy as np
import pandas as pd
import multiprocessing

# Set environment variables for locale
os.environ['LC_ALL'] = 'C.UTF-8'
os.environ['LANG'] = 'C.UTF-8'
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

def is_valid_fragment(fragment):
    # Validation logic for fragment
    return fragment != len(fragment) * "N"

def process_window(i, window_start, window_size, basename, algorithm, args, fasta_file, chromosome):
    pid = os.getpid()
    temp_filename = f"{basename}_{pid}.txt"

    if algorithm == "einverted":
        start_pos = i + 1  # EMBOSS uses 1-based indexing
        end_pos = i + window_size
        einverted_cmd = f"/usr/local/bin/einverted -sequence {fasta_file} -sbegin {start_pos} -send {end_pos} -gap {args.gaps} -threshold {args.score} -match {args.match} -mismatch {args.mismatch} -maxrepeat {args.max_span} -filter"
        process = subprocess.Popen(einverted_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, _ = process.communicate()
        ein_results = stdout.split("\n")
        parse_einverted_results(ein_results, window_start, window_size, basename, args, temp_filename, chromosome)

def parse_einverted_results(ein_results, window_start, window_size, basename, args, temp_filename, chromosome):
    with open(temp_filename, 'a') as temp_file:
        j = 0
        while j < len(ein_results) - 1:
            # Extract score and other details
            score_line = ein_results[j + 1].split()
            seq_i_full = ein_results[j + 2].split()
            #pairs = ein_results[j + 3]
            seq_j_full = ein_results[j + 4].split()
            #print(f"Score line: {score_line}")
            # print(f"Seq i full: {seq_i_full}")
            # print(f"Pairs: {pairs}")
            # print(f"Seq j full: {seq_j_full}")

            # Extracting score, raw match, percentage match, and gaps
            score = score_line[2]
            raw_match = score_line[3]
            matches, total = map(int, raw_match.split('/'))
            match_perc = round((matches / total) * 100, 2)    
            # find gaps one column from last column            
            gap_numb = score_line[-2]
            
            # Correctly calculate the genomic coordinates
            # Add the window_start to the relative positions
            i_start = int(seq_i_full[0]) + window_start
            i_end = int(seq_i_full[2]) + window_start
            j_start = int(seq_j_full[2]) + window_start
            j_end = int(seq_j_full[0]) + window_start

            # RNA folding and scoring
            i_seq = seq_i_full[1].replace("-", "").upper()
            j_seq = ''.join(reversed(seq_j_full[1].replace("-", ""))).upper()
            duplex_info = RNA.duplexfold(i_seq, j_seq)
            energy = round(float(duplex_info.energy), 2)
            structure = str(duplex_info.structure)
            pairs = int(structure.count('(') * 2)
            percent_paired = round(float(pairs / (len(structure) - 1)) * 100, 2)
            if match_perc < args.paired_cutoff:
                print(f"Skipping {i_start} to {j_end} due to low percentage of pairs: {percent_paired}")
                #print(f"{basename}\t{i_start}\t{i_end}\t{j_start}\t{j_end}\t{score}\t{raw_match}\t{match_perc}\t{gap_numb}\t{i_seq}\t{j_seq}\t{structure}\t{energy}\t{percent_paired}\n")
                j += 5
                continue
            
            # Writing results to the ein_results file
            temp_file.write(f"{chromosome}\t{i_start}\t{i_end}\t{j_start}\t{j_end}\t{score}\t{raw_match}\t{match_perc}\t{gap_numb}\t{i_seq}\t{j_seq}\t{structure}\t{energy}\t{percent_paired}\n")
            #print(f"{basename}\t{i_start}\t{i_end}\t{j_start}\t{j_end}\t{score}\t{raw_match}\t{match_perc}\t{gap_numb}\t{i_seq}\t{j_seq}\t{structure}\t{energy}\t{percent_paired}\n")
            
            # Increment j based on the structure of your einverted output
            j += 5

# Define the process_frame function
def process_frame(frame_start, frame_step_size, end_coordinate, window_size, basename, algorithm, args, fasta_file, chromosome, pool):
    for start in range(frame_start, end_coordinate, frame_step_size):
        window_start = start
        end = min(start + window_size, end_coordinate)
        pool.apply_async(process_window, (start, window_start, window_size, basename, algorithm, args, fasta_file, chromosome))


def main():
    ### Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',  type=str,
                        help='input filename')
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('-s', '--step', type=int, default=150,
                        help='Step size; default = 500')
    parser.add_argument('-w', type=int, default=10000,
                        help='Window size; default = 10000')
    parser.add_argument('--max_span', type=int, default=10000,
                        help='Max span of inverted repeat; default = 10000')
    parser.add_argument('--score', type=int, default=50,
                        help='Minimum score threshold for inverted repeat; Default = 50')
    parser.add_argument('--min', type=int, default=30,
                        help='Minimum length of inverted repeat; Default = 30')
    parser.add_argument('--max', type=int, default=10000,
                        help='Max length of inverted repeat; Default = 10000')
    parser.add_argument('--gaps', type=int, default=12,
                        help='Gap penalty')
    parser.add_argument('--start', type=int, default=0,
                        help='Starting coordinate for scan; Default = 0')
    parser.add_argument('--end', type=int, default=0,
                        help='Ending coordinate for scan; Default = 0')
    parser.add_argument('-x', '--mismatch', type=int, default=-4,
            help='Mismatch score')
    parser.add_argument('--match', type=int, default=3,
            help='Match score')
    parser.add_argument('--paired_cutoff', type=int, default=70,
                        help='Cutoff to ignore sturctures with low percentage of pairs; Default <70')
    parser.add_argument('--algorithm', type=str, default="einverted",
            help='Inverted repeat finding algorithm (einverted or iupacpal)')
    parser.add_argument('--reverse', action='store_true', default=False,
                        help='Use this option if running on reverse strand')
    parser.add_argument('--chr', type=str, required=True,
                        help='Chromosome name, if chromosome name in header type "header"')
    parser.add_argument('-c', '--cpus', type=int, default=4,
                        help='Number of cpus to use; Default = 4')
    
    args = parser.parse_args()
    chromosome = args.chr
    end_coordinate = int(args.end)
    fasta_file = args.filename
    cpu_count = args.cpus
    step_size = args.step

    with open(args.filename) as f:
        # Create a pool of workers for multiprocessing, starting with 2 workers
        pool = multiprocessing.Pool(cpu_count)
        
        # Process each sequence
        tasks = []

        for cur_record in SeqIO.parse(f, "fasta"):            
            # Check if chromosome is 'header'
            if args.chr == "header":
                chromosome = cur_record.name
            else:
                chromosome = args.chr

            # Determine strand and set up basename
            strand = "-" if args.reverse else "+"
            if args.reverse:
                # Reverse complement the sequence and print to new fasta file
                with open(f"{args.filename.split('.')[0]}.{chromosome}.reverse.fasta", 'w+') as reverse_file:
                    reverse_file.write(f">{cur_record.name}\n")
                    reverse_file.write(str(cur_record.seq.reverse_complement()))
                # Set filename to the new reverse complemented fasta file
                fasta_file = f"{args.filename.split('.')[0]}.{chromosome}.reverse.fasta"

            basename = f"{args.filename.split('.')[0]}.{chromosome}.{'reverse' if args.reverse else 'forward'}_win{args.w}_step{args.step}_start{args.start}_score{args.score}"

            # Set up result files
            with open(f"{basename}.ein_results.txt", 'w+') as results_file:
                results_file.write("Chromosome\tStrand\tScore\tRawMatch\tPercMatch\tGaps\ti_start\ti_end\tj_start\tj_end\ti_seq\tj_seq\tstructure\tdG(kcal/mol)\tpercent_paired\n")
            
            with open(f"{basename}.dsRNApredictions.bp", 'w+') as bp_file:
                # Example header - adjust based on your requirements
                bp_file.write("# Base Pair Predictions\n")
                bp_file.write("# Format: sequence_id\tstart\tend\n")

            # Process each sequence
            end_coordinate = args.end if args.end != 0 else len(cur_record.seq)
            
            frame_step_size = step_size * cpu_count
            for cpu_index in range(cpu_count):
                frame_start = cpu_index * step_size

                # Start processing at each frame and jump by frame_step_size
                for start in range(frame_start, end_coordinate, frame_step_size):
                    end = min(start + args.w, end_coordinate)
                    tasks.append(pool.apply_async(process_window, (start, args.w, basename, args.algorithm, args, fasta_file, chromosome)))
            
            # Close the pool
            pool.close()
            pool.join()

            # Merge, sort results, and clean up temporary files
            temp_files = glob.glob(f"{basename}_*.txt")
            merged_filename = f"{basename}_merged_results.txt"
            with open(merged_filename, 'w') as merged_file:
                merged_file.write("Chromosome\tStrand\tScore\tRawMatch\tPercMatch\tGaps\ti_start\ti_end\tj_start\tj_end\ti_seq\tj_seq\tstructure\tdG(kcal/mol)\tpercent_paired\n")
                all_data = []
                for temp_file in temp_files:
                    with open(temp_file, 'r') as file:
                        lines = file.readlines()
                        if lines:  # Check if file is not empty
                            all_data.extend(lines[1:])  # Skip header and add the rest

                    os.remove(temp_file)

                all_data.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1]), int(x.split("\t")[2])))
                for line in all_data:
                    merged_file.write(line)

# Run the main function
if __name__ == "__main__":
    main()


