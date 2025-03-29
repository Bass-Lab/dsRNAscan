#!/usr/bin/env python3
import os
import locale
import glob
from Bio import SeqIO
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

# Determine the directory of this script and set the local path for einverted
script_dir = os.path.dirname(os.path.abspath(__file__))
einverted_bin = os.path.join(script_dir, "./tools", "einverted")

def is_valid_fragment(fragment):
    # Validation logic for fragment
    return fragment != len(fragment) * "N"


def generate_bp_file(input_file, output_file):
    """
    Generate a BP file from the merged dsRNA results file using the correct format for IGV.
    Uses strand information and percent_paired for color coding.
    
    Args:
        input_file (str): Path to the merged results file
        output_file (str): Path to the output BP file
    """
    print(f"Reading data from {input_file}")
    
    try:
        # Check if the file exists and is not empty
        if not os.path.exists(input_file):
            print(f"Error: File {input_file} does not exist")
            return
            
        if os.path.getsize(input_file) == 0:
            print(f"Error: File {input_file} is empty")
            return
            
        # Read the merged results file with better error handling
        try:
            df = pd.read_csv(input_file, sep="\t")
        except pd.errors.EmptyDataError:
            print(f"Error: No data found in {input_file}")
            return
        except Exception as e:
            print(f"Error reading file {input_file}: {str(e)}")
            return
            
        # Check if DataFrame is empty
        if df.empty:
            print(f"Warning: No data found in {input_file}")
            return
        
        # Print the column names for debugging
        print(f"Columns in file: {', '.join(df.columns)}")
        
        # Verify required columns exist
        required_cols = ["Chromosome", "i_start", "i_end", "j_start", "j_end", "percent_paired"]
        
        # Handle inconsistent column naming
        column_mappings = {
            "Chromosome": ["chromosome", "chr", "chrom"],
            "i_start": ["start1", "start_1", "left_start"],
            "i_end": ["end1", "end_1", "left_end"],
            "j_start": ["start2", "start_2", "right_start"],
            "j_end": ["end2", "end_2", "right_end"],
            "percent_paired": ["percpaired", "perc_paired", "percentpaired", "percent_match", "PercMatch"]
        }
        
        # Handle strand column specifically - it might be missing but we can default it
        has_strand = "Strand" in df.columns
        if not has_strand:
            for alt_name in ["strand", "str"]:
                if alt_name in df.columns:
                    df = df.rename(columns={alt_name: "Strand"})
                    has_strand = True
                    break
                    
        # If still no strand column, add default
        if not has_strand:
            print("No strand column found, defaulting to '+' strand")
            df["Strand"] = "+"
        
        # Check for missing columns and try to use alternatives
        missing_cols = []
        for col in required_cols:
            if col not in df.columns:
                # Try alternative names
                found = False
                if col in column_mappings:
                    for alt_col in column_mappings[col]:
                        if alt_col in df.columns:
                            df = df.rename(columns={alt_col: col})
                            found = True
                            break
                
                if not found:
                    missing_cols.append(col)
        
        if missing_cols:
            print(f"Error: Missing required columns: {', '.join(missing_cols)}")
            print(f"Available columns: {', '.join(df.columns)}")
            return
        
        # Create a new BP file
        with open(output_file, 'w') as bp_file:
            # Write header with color definitions according to the correct BP format
            bp_file.write("color:\t100\t149\t237\t70-80% paired (forward strand)\n")
            bp_file.write("color:\t65\t105\t225\t80-90% paired (forward strand)\n")
            bp_file.write("color:\t0\t0\t139\t90-100% paired (forward strand)\n")
            bp_file.write("color:\t205\t92\t92\t70-80% paired (reverse strand)\n")
            bp_file.write("color:\t178\t34\t34\t80-90% paired (reverse strand)\n")
            bp_file.write("color:\t139\t0\t0\t90-100% paired (reverse strand)\n")
            
            # Process each row
            for idx, row in df.iterrows():
                # Get chromosome and positions
                chrom = row["Chromosome"]
                
                # Get the strand and convert to "+" or "-" if needed
                strand = row.get("Strand", "+")
                if strand not in ["+", "-"]:
                    # Handle numeric or other formats
                    if strand == "1" or str(strand).lower() == "forward":
                        strand = "+"
                    elif strand == "-1" or str(strand).lower() == "reverse":
                        strand = "-"
                    else:
                        strand = "+"  # Default
                
                # Get percent paired - try different column names if needed
                if "percent_paired" in row:
                    percent_paired = row["percent_paired"]
                elif "PercMatch" in row:
                    percent_paired = row["PercMatch"]
                else:
                    percent_paired = 75.0  # Default
                
                # Convert percent_paired to float if it's not already
                try:
                    if isinstance(percent_paired, str):
                        percent_paired = float(percent_paired.replace('%', ''))
                    else:
                        percent_paired = float(percent_paired)
                except ValueError:
                    print(f"Warning: Could not parse percent_paired value '{percent_paired}', defaulting to 75")
                    percent_paired = 75.0
                
                # Determine color index based on strand and percent_paired
                # Color indices match the header order (0-5)
                if strand == "+":
                    # Forward strand
                    if percent_paired >= 90:
                        color_idx = 2  # dark blue
                    elif percent_paired >= 80:
                        color_idx = 1  # royal blue
                    else:
                        color_idx = 0  # cornflower blue
                else:
                    # Reverse strand
                    if percent_paired >= 90:
                        color_idx = 5  # dark red
                    elif percent_paired >= 80:
                        color_idx = 4  # firebrick
                    else:
                        color_idx = 3  # indian red
                
                # Get the coordinates for both arms
                try:
                    i_start = int(row["i_start"])
                    i_end = int(row["i_end"])
                    j_start = int(row["j_start"])
                    j_end = int(row["j_end"])
                except (ValueError, TypeError) as e:
                    print(f"Warning: Could not parse coordinate values for row {idx}, skipping: {e}")
                    continue
                
                # Write the BP record with coordinates from both arms forming a pair
                # Format: <chrom> <left_start> <left_end> <right_start> <right_end> <color_idx>
                bp_file.write(f"{chrom}\t{i_start}\t{i_end}\t{j_start}\t{j_end}\t{color_idx}\n")
            
            print(f"Successfully wrote BP file to {output_file}")
    except Exception as e:
        print(f"Error generating BP file: {str(e)}")
        import traceback
        traceback.print_exc()

def fix_merge_temp_files(basename):
    """
    Fix the merging of temporary files into a single output file.
    
    Args:
        basename (str): Base name for the temporary files pattern
    
    Returns:
        str: Path to the merged results file
    """
    temp_files = glob.glob(f"{basename}_*.txt")
    merged_filename = f"{basename}_merged_results.txt"

    # Check if any temp files exist
    if not temp_files:
        print(f"Warning: No temporary files found matching pattern {basename}_*.txt")
        # Create an empty output file with headers
        with open(merged_filename, 'w') as merged_file:
            merged_file.write("Chromosome\ti_start\ti_end\tj_start\tj_end\teff_i_start\teff_i_end\teff_j_start\teff_j_end\tScore\tRawMatch\tPercMatch\tGaps\ti_seq\tj_seq\tstructure\tdG(kcal/mol)\tpercent_paired\n")
        return merged_filename
    else:
        print(f"Found {len(temp_files)} temporary files to merge")
        
        # Initialize an empty DataFrame with the expected columns
        column_names = ["Chromosome", "i_start", "i_end", "j_start", "j_end", 
                       "eff_i_start", "eff_i_end", "eff_j_start", "eff_j_end", 
                       "Score", "RawMatch", "PercMatch", "Gaps", 
                       "i_seq", "j_seq", "structure", "dG(kcal/mol)", "percent_paired"]
        
        all_dfs = []
        
        # Process each temp file individually to better handle errors
        for temp_file in temp_files:
            try:
                # Check if file is empty
                if os.path.getsize(temp_file) == 0:
                    print(f"Skipping empty file: {temp_file}")
                    continue
                    
                # Try to read the file with various approaches
                try:
                    # First try reading without header assumptions
                    df = pd.read_csv(temp_file, sep="\t", header=None)
                    
                    # If we got here, the file was read successfully
                    if len(df.columns) == len(column_names):
                        df.columns = column_names
                        all_dfs.append(df)
                    else:
                        print(f"Warning: File {temp_file} has {len(df.columns)} columns, expected {len(column_names)}")
                        print(f"First row: {df.iloc[0].tolist()}")
                        
                        # Try to handle common cases - first line might be header
                        if len(df.columns) == 1 and isinstance(df.iloc[0, 0], str) and "\t" in df.iloc[0, 0]:
                            print(f"Attempting to parse as TSV with embedded tabs")
                            # Re-read with pandas' flexible parsing
                            df = pd.read_csv(temp_file, sep=None, engine='python')
                            if len(df.columns) == len(column_names):
                                df.columns = column_names
                                all_dfs.append(df)
                except Exception as e:
                    print(f"Error reading file {temp_file}: {str(e)}")
                    
                    # Try a more basic approach - read line by line
                    print("Attempting manual parsing...")
                    manual_rows = []
                    with open(temp_file, 'r') as f:
                        for line in f:
                            if line.strip() and not line.startswith('#'):
                                fields = line.strip().split('\t')
                                if len(fields) == len(column_names):
                                    manual_rows.append(fields)
                    
                    if manual_rows:
                        print(f"Manually parsed {len(manual_rows)} rows from {temp_file}")
                        df = pd.DataFrame(manual_rows, columns=column_names)
                        all_dfs.append(df)
                    
            except Exception as e:
                print(f"Failed to process file {temp_file}: {str(e)}")
        
        if all_dfs:
            # Combine all successfully read DataFrames
            df = pd.concat(all_dfs, ignore_index=True)
            
            # Handle empty DataFrame case
            if df.empty:
                print("Warning: No data was successfully read from temp files")
                # Create empty file with headers
                with open(merged_filename, 'w') as merged_file:
                    merged_file.write("\t".join(column_names) + "\n")
            else:
                # Convert coordinate columns to numeric for proper sorting
                for col in ["i_start", "i_end", "j_start", "j_end", "eff_i_start", 
                            "eff_i_end", "eff_j_start", "eff_j_end"]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                
                # Drop duplicate rows
                df = df.drop_duplicates()
                
                # Convert all sequence columns to uppercase RNA
                df['i_seq'] = df['i_seq'].str.upper().str.replace("T", "U")
                df['j_seq'] = df['j_seq'].str.upper().str.replace("T", "U")
                
                # Sort the DataFrame by chromosome and numeric coordinates
                df = df.sort_values(by=["Chromosome", "i_start", "i_end"])
                
                df['i_eff_seq'] = df.apply(lambda row: row['i_seq'][row['eff_i_start']-1:row['eff_i_end']], axis=1)
                df['j_eff_seq'] = df.apply(lambda row: row['j_seq'][row['eff_j_start']-1:row['eff_j_end']], axis=1)

                # Sort the DataFrame by chromosome and numeric coordinates
                df = df.sort_values(by=["Chromosome", "i_start", "i_end"])
                
                # Write out the merged results
                df.to_csv(merged_filename, sep="\t", index=False)
                print(f"Successfully wrote {len(df)} records to {merged_filename}")
        else:
            print("Warning: No data could be read from any temp files")
            # Create empty file with headers
            with open(merged_filename, 'w') as merged_file:
                merged_file.write("\t".join(column_names) + "\n")
        
        return merged_filename
    
def predict_hybridization(seq1, seq2):
    """
    Predict RNA-RNA interactions using RNAduplex.
    
    Args:
        seq1 (str): First RNA sequence
        seq2 (str): Second RNA sequence
    
    Returns:
        str: Raw output from RNAduplex containing structure and energy
    """
    process = subprocess.Popen(
        ["RNAduplex"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    input_data = f"{seq1}\n{seq2}\n".encode("utf-8")
    stdout, stderr = process.communicate(input_data)
    #print(f"[DEBUG] RNAduplex stdout: {stdout.decode('utf-8')}")
    if stderr:
         print(f"[DEBUG] RNAduplex stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8").strip()

def parse_rnaduplex_output(output):
    """
    Parse the output from RNAduplex.
    
    Example output: "(((...)))&(((...)))  1,9 : 3,11  (-10.40)"
    
    Args:
        output (str): Output string from RNAduplex
    
    Returns:
        tuple: (structure, indices_seq1, indices_seq2, energy)
    """
    try:
        # print(f"[DEBUG] Parsing RNAduplex output: {output}")
        parts = output.split()
        # print(f"[DEBUG] RNAduplex parts: {parts}")
        
        # Handle empty or invalid output
        if not parts or len(parts) < 4:
            print(f"Warning: Invalid RNAduplex output: {output}")
            return "", [0, 0], [0, 0], 0.0
        
        structure = parts[0]
        
        # Parse indices more safely
        try:
            indices_seq1 = [int(x) for x in parts[1].split(',')]
            if len(indices_seq1) != 2:
                raise ValueError(f"Expected 2 indices for seq1, got {len(indices_seq1)}")
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse indices_seq1 from {parts[1] if len(parts) > 1 else 'missing'}: {e}")
            indices_seq1 = [1, 1]  # Default to position 1
            
        try:
            indices_seq2 = [int(x) for x in parts[3].split(',')]
            if len(indices_seq2) != 2:
                raise ValueError(f"Expected 2 indices for seq2, got {len(indices_seq2)}")
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse indices_seq2 from {parts[3] if len(parts) > 3 else 'missing'}: {e}")
            indices_seq2 = [1, 1]  # Default to position 1
        
        # Extract energy from the output
        energy = None
        if len(parts) > 4:
            # Energy is typically in the format (-10.40)
            energy_str = parts[4].strip('()')
            try:
                energy = float(energy_str)
            except ValueError:
                print(f"Warning: Could not parse energy value '{energy_str}' from RNAduplex output")
                energy = 0.0
        else:
            energy = 0.0
        
        return structure, indices_seq1, indices_seq2, energy
    except Exception as e:
        print(f"Error parsing RNAduplex output: {e}")
        return "", [1, 1], [1, 1], 0.0

def safe_extract_effective_seq(row, seq_col, start_col, end_col):
    """
    Safely extract a subsequence based on the effective indices.
    Handles type conversion, boundary checking, and exceptions.
    
    Args:
        row: DataFrame row
        seq_col: Column name for the sequence
        start_col: Column name for the start index
        end_col: Column name for the end index
        
    Returns:
        str: The extracted subsequence or the full sequence if extraction fails
    """
    try:
        # Make sure we have non-empty sequence
        if not row[seq_col] or pd.isna(row[seq_col]):
            return ""
            
        # Make sure we have valid numbers for indices
        if pd.isna(row[start_col]) or pd.isna(row[end_col]):
            return row[seq_col]
            
        # Make sure we have integers for slicing
        start_idx = int(float(row[start_col])) - 1  # Convert to 0-based index
        end_idx = int(float(row[end_col]))
        
        # Make sure indices are valid for the sequence
        if start_idx < 0:
            start_idx = 0
            
        seq = str(row[seq_col])
        if end_idx > len(seq):
            end_idx = len(seq)
        
        # Skip extraction if indices are invalid
        if start_idx >= end_idx or start_idx >= len(seq):
            return seq
            
        # Return the slice
        return seq[start_idx:end_idx]
    except (ValueError, TypeError, IndexError) as e:
        print(f"Warning: Could not extract effective sequence: {e}. Using full sequence instead.")
        # Return the original sequence as fallback
        return str(row[seq_col]) if row[seq_col] and not pd.isna(row[seq_col]) else ""

def process_window(i, window_start, window_size, basename, algorithm, args, fasta_file, chromosome, strand):
    """Process a genomic window to identify dsRNA structures"""
    pid = os.getpid()
    temp_filename = f"{basename}_{pid}.txt"

    if algorithm == "einverted":
        start_pos = i + 1  # EMBOSS uses 1-based indexing
        end_pos = i + window_size
        einverted_cmd = f"{einverted_bin} -sequence {fasta_file} -sbegin {start_pos} -send {end_pos} -gap {args.gaps} -threshold {args.score} -match {args.match} -mismatch {args.mismatch} -maxrepeat {args.max_span} -filter"
        
        process = subprocess.Popen(einverted_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, _ = process.communicate()
        
        ein_results = stdout.split("\n")
        parse_einverted_results(ein_results, window_start, window_size, basename, args, temp_filename, chromosome, strand)
        
def parse_einverted_results(ein_results, window_start, window_size, basename, args, temp_filename, chromosome, strand):
    """Parse results from einverted and save to temporary file, including strand information"""
    with open(temp_filename, 'a') as temp_file:
        j = 0
        while j < len(ein_results) - 1:
            # Skip if we don't have at least 5 lines in the current result block
            if j + 4 >= len(ein_results):
                break
                
            # Extract score and other details
            try:
                score_line = ein_results[j + 1].split()
                seq_i_full = ein_results[j + 2].split()
                seq_j_full = ein_results[j + 4].split()
                
                # Skip if we don't have enough data
                if len(score_line) < 4 or len(seq_i_full) < 3 or len(seq_j_full) < 3:
                    j += 5
                    continue
                    
                # Extracting score, raw match, percentage match, and gaps
                score = score_line[2]
                raw_match = score_line[3]
                matches, total = map(int, raw_match.split('/'))
                match_perc = round((matches / total) * 100, 2)    
                # find gaps one column from last column            
                gap_numb = score_line[-2]
                
                # Calculate the genomic coordinates directly from the einverted output
                i_start = int(seq_i_full[0])
                i_end = int(seq_i_full[2])
                j_start = int(seq_j_full[2])
                j_end = int(seq_j_full[0])
                
                # Double-check the coordinates are in the correct order
                if i_start > i_end or j_start > j_end or i_start > j_start or i_end > j_end:
                    print(f"Warning: Coordinates not in correct order for {i_start} to {j_end}. Sorting...")
                    coords = sorted([i_start, i_end, j_start, j_end])
                    i_start, i_end, j_start, j_end = coords[0], coords[1], coords[2], coords[3]
                
                # RNA folding and scoring
                i_seq = seq_i_full[1].replace("-", "").upper()
                j_seq = ''.join(reversed(seq_j_full[1].replace("-", ""))).upper()
                output = predict_hybridization(i_seq, j_seq)
                structure, indices_seq1, indices_seq2, energy = parse_rnaduplex_output(output)
                
                # Skip if we got empty results from RNAduplex
                if not structure:
                    j += 5
                    continue
                    
                eff_i = (indices_seq1[0], indices_seq1[1])
                eff_j = (indices_seq2[0], indices_seq2[1])
                pairs = int(structure.count('(') * 2)
                
                # Calculate percent_paired safely
                try:
                    percent_paired = round(float(pairs / (len(structure) - 1)) * 100, 2)
                except (ZeroDivisionError, ValueError):
                    percent_paired = 0
                
                if match_perc < args.paired_cutoff:
                    print(f"Skipping {i_start} to {j_end} due to low percentage of pairs: {percent_paired}")
                    j += 5
                    continue
                
                # Writing results to the ein_results file - include strand information
                temp_file.write(f"{chromosome}\t{strand}\t{score}\t{raw_match}\t{match_perc}\t{gap_numb}\t{i_start}\t{i_end}\t{j_start}\t{j_end}\t{eff_i[0]}\t{eff_i[1]}\t{eff_j[0]}\t{eff_j[1]}\t{i_seq}\t{j_seq}\t{structure}\t{energy}\t{percent_paired}\n")
            except Exception as e:
                print(f"Error processing result block at index {j}: {str(e)}")
            
            # Increment j based on the structure of your einverted output
            j += 5
            
# Define the process_frame function
def process_frame(frame_start, frame_step_size, end_coordinate, window_size, basename, algorithm, args, fasta_file, chromosome, pool):
    for start in range(frame_start, end_coordinate, frame_step_size):
        window_start = start
        end = min(start + window_size, end_coordinate)
        pool.apply_async(process_window, (start, window_start, window_size, basename, algorithm, args, fasta_file, chromosome))
        # For debugging, run the process_window function directly
        # process_window(start, start, args.w, basename, args.algorithm, args, fasta_file, chromosome)

def main():
    ### Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',  type=str,
                        help='input filename')
    parser.add_argument('-t', type=int, default=37,
                        help='Folding temperature in celsius; default = 37C')
    parser.add_argument('-s', '--step', type=int, default=150,
                        help='Step size; default = 150')
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
    parser.add_argument('--clean', action='store_false', default=True,
                    help='Clean up temporary files after processing')
    
    args = parser.parse_args()
    chromosome = args.chr
    end_coordinate = int(args.end)
    fasta_file = args.filename
    cpu_count = args.cpus
    step_size = args.step

    with open(args.filename) as f:
        # Create a pool of workers for multiprocessing, starting with 2 workers
        # pool = multiprocessing.Pool(cpu_count)
        
        # Process each sequence
        tasks = []

        for cur_record in SeqIO.parse(f, "fasta"): 
            # Correct for starting coordinate
            print(f"Processing sequence: {cur_record.name}")   
            # Print the sequence length
            #print(f"Sequence length: {len(cur_record.seq)}")
            
            # Convert to RNA uppercase
            #cur_record.seq = cur_record.seq.transcribe().upper()
            
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
                    reverse_file.write(str(cur_record.seq.reverse_complement().transcribe().upper()))
                # Set filename to the new reverse complemented fasta file
                fasta_file = f"{args.filename.split('.')[0]}.{chromosome}.reverse.fasta"
            else:
                # Convert forward strand to RNA and print to new fasta file
                with open(f"{args.filename.split('.')[0]}.{chromosome}.forward.fasta", 'w+') as forward_file:
                    forward_file.write(f">{cur_record.name}\n")
                    forward_file.write(str(cur_record.seq.transcribe().upper()))
                # Set filename to the new RNA fasta file
                fasta_file = f"{args.filename.split('.')[0]}.{chromosome}.forward.fasta"
            
            # Set up basename for output files
            basename = f"{args.filename.split('.')[0]}.{chromosome}.{'reverse' if args.reverse else 'forward'}_win{args.w}_step{args.step}_start{args.start}_score{args.score}"

            # Set up result files
            with open(f"{basename}.ein_results.txt", 'w+') as results_file:
                results_file.write("Chromosome\tStrand\tScore\tRawMatch\tPercMatch\tGaps\ti_start\ti_end\tj_start\tj_end\teff_i_start\teff_i_end\teff_j_start\teff_j_end\ti_seq\tj_seq\tstructure\tdG(kcal/mol)\tpercent_paired\n")
            
            # with open(f"{basename}.dsRNApredictions.bp", 'w+') as bp_file:
            #     # Example header - adjust based on your requirements
            #     bp_file.write("# Base Pair Predictions\n")
            #     bp_file.write("# Format: sequence_id\tstart\tend\n")

            # Process each sequence
            end_coordinate = args.end if args.end != 0 else len(cur_record.seq)
            seq_length = end_coordinate - args.start

            # Determine if the sequence is short (less than window size)
            is_short_sequence = seq_length < args.w

            # Print what we're scanning now
            if is_short_sequence:
                print(f"Short sequence detected: {cur_record.name} length {seq_length} bp")
                print(f"Using single window approach for the entire sequence")
                
                # Just process the entire sequence as one window
                process_window(args.start, args.start, seq_length, basename, args.algorithm, 
                            args, fasta_file, chromosome, strand)
            else:
                # Normal processing for longer sequences
                print(f"Scanning {cur_record.name} from {args.start} to {end_coordinate} with window size {args.w} and step size {args.step}")
                
                # Create a pool of workers for multiprocessing 
                pool = multiprocessing.Pool(cpu_count)
                tasks = []
                
                
                # Use multiprocessing for longer sequences
                frame_step_size = step_size * cpu_count
                for cpu_index in range(cpu_count):
                    # Start from the specified start coordinate plus the CPU's offset
                    frame_start = args.start + (cpu_index * step_size)

                    # Start processing at each frame and jump by frame_step_size
                    for start in range(frame_start, end_coordinate, frame_step_size):
                        window_end = min(start + args.w, end_coordinate)
                        window_size = window_end - start
                        
                        # Only process if we have a meaningful window
                        if window_size >= args.min:
                            tasks.append(pool.apply_async(process_window, 
                                        (start, start, window_size, basename, args.algorithm, 
                                        args, fasta_file, chromosome, strand)))
                # Close the pool
                pool.close()
                pool.join()

            # Merge, sort results, and clean up temporary files
            temp_files = glob.glob(f"{basename}_*.txt")
            merged_filename = f"{basename}_merged_results.txt"

            # Check if any temp files exist
            if not temp_files:
                print(f"Warning: No temporary files found matching pattern {basename}_*.txt")
                # Create an empty output file with headers to avoid downstream errors
                with open(merged_filename, 'w') as merged_file:
                    merged_file.write("Chromosome\ti_start\ti_end\tj_start\tj_end\teff_i_start\teff_i_end\teff_j_start\teff_j_end\tScore\tRawMatch\tPercMatch\tGaps\ti_seq\tj_seq\tstructure\tdG(kcal/mol)\tpercent_paired\n")
            else:
                print(f"Found {len(temp_files)} temporary files to merge")
                
                # Initialize an empty DataFrame with the expected columns
                column_names = ["Chromosome", "Strand", "Score", "RawMatch", "PercMatch", "Gaps", 
                            "i_start", "i_end", "j_start", "j_end", "eff_i_start", "eff_i_end", 
                            "eff_j_start", "eff_j_end", "i_seq", "j_seq", "structure", 
                            "dG(kcal/mol)", "percent_paired"]

                all_dfs = []
                
                # Process each temp file individually to better handle errors
                for temp_file in temp_files:
                    try:
                        # Check if file is empty
                        if os.path.getsize(temp_file) == 0:
                            print(f"Skipping empty file: {temp_file}")
                            continue
                            
                        # Try to read the file with various approaches
                        try:
                            # First try reading without header assumptions
                            df = pd.read_csv(temp_file, sep="\t", header=None)
                            
                            # If we got here, the file was read successfully
                            if len(df.columns) == len(column_names):
                                df.columns = column_names
                                all_dfs.append(df)
                            else:
                                print(f"Warning: File {temp_file} has {len(df.columns)} columns, expected {len(column_names)}")
                                print(f"First row: {df.iloc[0].tolist()}")
                                
                                # Try to handle common cases - first line might be header
                                if len(df.columns) == 1 and isinstance(df.iloc[0, 0], str) and "\t" in df.iloc[0, 0]:
                                    print(f"Attempting to parse as TSV with embedded tabs")
                                    # Re-read with pandas' flexible parsing
                                    df = pd.read_csv(temp_file, sep=None, engine='python')
                                    if len(df.columns) == len(column_names):
                                        df.columns = column_names
                                        all_dfs.append(df)
                        except Exception as e:
                            print(f"Error reading file {temp_file}: {str(e)}")
                            
                            # Try a more basic approach - read line by line
                            print("Attempting manual parsing...")
                            manual_rows = []
                            with open(temp_file, 'r') as f:
                                for line in f:
                                    if line.strip() and not line.startswith('#'):
                                        fields = line.strip().split('\t')
                                        if len(fields) == len(column_names):
                                            manual_rows.append(fields)
                            
                            if manual_rows:
                                print(f"Manually parsed {len(manual_rows)} rows from {temp_file}")
                                df = pd.DataFrame(manual_rows, columns=column_names)
                                all_dfs.append(df)
                            
                    except Exception as e:
                        print(f"Failed to process file {temp_file}: {str(e)}")
                
                
                if all_dfs:
                    # Combine all successfully read DataFrames
                    df = pd.concat(all_dfs, ignore_index=True)
                    
                    # Handle empty DataFrame case
                    if df.empty:
                        print("Warning: No data was successfully read from temp files")
                        # Create empty file with headers
                        with open(merged_filename, 'w') as merged_file:
                            merged_file.write("\t".join(column_names) + "\n")
                    else:
                        # Convert coordinate columns to numeric for proper sorting
                        for col in ["i_start", "i_end", "j_start", "j_end", "eff_i_start", 
                                    "eff_i_end", "eff_j_start", "eff_j_end"]:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                                        
                        # Drop duplicate rows
                        df = df.drop_duplicates()
                        
                        # Convert all sequence columns to uppercase RNA
                        df['i_seq'] = df['i_seq'].str.upper().str.replace("T", "U")
                        df['j_seq'] = df['j_seq'].str.upper().str.replace("T", "U")
                        
                        # Sort the DataFrame by chromosome and numeric coordinates
                        df = df.sort_values(by=["Chromosome", "i_start", "i_end"])
                        
                        # Safe extraction of effective sequences
                        def safe_extract_effective_seq(row, seq_col, start_col, end_col):
                            try:
                                # Make sure we have integers for slicing
                                start_idx = int(row[start_col]) - 1  # Convert to 0-based index
                                end_idx = int(row[end_col])
                                
                                # Make sure indices are valid for the sequence
                                if start_idx < 0:
                                    start_idx = 0
                                    
                                seq = str(row[seq_col])
                                if end_idx > len(seq):
                                    end_idx = len(seq)
                                    
                                # Return the slice
                                return seq[start_idx:end_idx]
                            except (ValueError, TypeError, IndexError) as e:
                                print(f"Warning: Could not extract effective sequence: {e}. Using full sequence instead.")
                                return row[seq_col]
                        

                        # With these lines:
                        df['i_eff_seq'] = df.apply(
                            lambda row: safe_extract_effective_seq(row, 'i_seq', 'eff_i_start', 'eff_i_end'), 
                            axis=1
                        )
                        df['j_eff_seq'] = df.apply(
                            lambda row: safe_extract_effective_seq(row, 'j_seq', 'eff_j_start', 'eff_j_end'), 
                            axis=1
                        )

                        # Add FORNA visualization links (with error checking)
                        df['FORNA_Link'] = (
                            "http://rna.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=" + 
                            df['i_eff_seq'].fillna('') + df['j_eff_seq'].fillna('') + 
                            "&structure=" + 
                            df['structure'].fillna('').str.replace("&", "")  
                        )

                        # Write out the merged results
                        df.to_csv(merged_filename, sep="\t", index=False)
                        print(f"Successfully wrote {len(df)} records to {merged_filename}")
                
                else:
                    print("Warning: No data could be read from any temp files")
                    # Create empty file with headers
                    with open(merged_filename, 'w') as merged_file:
                        merged_file.write("\t".join(column_names) + "\n")

            # Clean up temp files if desired (uncomment to enable)
            # for temp_file in temp_files:
            #     try:
            #         os.remove(temp_file)
            #     except Exception as e:
            #         print(f"Failed to remove temp file {temp_file}: {str(e)}")

            # Now generate the BP file
            try:
                # Use the function from the previous script to generate BP file
                generate_bp_file(merged_filename, f"{basename}.dsRNApredictions.bp")
            except NameError:
                print("BP file generation function not defined. Please add the generate_bp_file function to your script.")
            
            # Clean up temporary files if requested
            if args.clean:
                for temp_file in temp_files:
                    try:
                        os.remove(temp_file)
                        print(f"Removed temporary file: {temp_file}")
                    except Exception as e:
                        print(f"Failed to remove temp file {temp_file}: {str(e)}")
                        
            # Print file names and paths
            print(f"Results written to {merged_filename}")
            print(f"Base Pair predictions written to {basename}.dsRNApredictions.bp")
            
# Run the main function
if __name__ == "__main__":
    main()
