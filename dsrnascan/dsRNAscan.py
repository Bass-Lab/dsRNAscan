#!/usr/bin/env python3
"""
dsRNAscan - A tool for genome-wide prediction of double-stranded RNA structures
Copyright (C) 2024 Bass Lab
Version: 0.3.0
"""

__version__ = '0.3.0'
__author__ = 'Bass Lab'

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
import gzip
from datetime import datetime
from queue import Empty

# Set environment variables for locale
os.environ['LC_ALL'] = 'C.UTF-8'
os.environ['LANG'] = 'C.UTF-8'
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

# Determine the directory of this script and set the local path for einverted
script_dir = os.path.dirname(os.path.abspath(__file__))

# Try multiple locations for einverted binary
possible_paths = [
    os.path.join(script_dir, "tools", "einverted"),  # Development location
    os.path.join(script_dir, "..", "tools", "einverted"),  # Installed location
    os.path.join(os.path.dirname(script_dir), "tools", "einverted"),  # Alternative installed
    "/usr/local/bin/einverted",  # System installation
    "/usr/bin/einverted",  # System installation
]

einverted_bin = None
for path in possible_paths:
    if os.path.exists(path) and os.access(path, os.X_OK):
        einverted_bin = path
        break

if not einverted_bin:
    # Last resort: check if einverted is in PATH
    from shutil import which
    einverted_bin = which("einverted")
    
if not einverted_bin:
    einverted_bin = os.path.join(script_dir, "tools", "einverted")  # Default for error message

def smart_open(filename, mode='rt'):
    """
    Open a file, automatically detecting if it's gzipped based on extension.
    
    Args:
        filename: Path to file
        mode: Mode to open file (default 'rt' for reading text)
        
    Returns:
        File handle (either regular file or gzip file)
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

# Check if einverted binary exists
if not os.path.exists(einverted_bin):
    print(f"Error: einverted binary not found at {einverted_bin}")
    print("Please ensure the einverted tool is installed in the 'tools' subdirectory.")
    print("You can install it by running the installation script or downloading from EMBOSS.")
    sys.exit(1)

# Check if einverted is executable
if not os.access(einverted_bin, os.X_OK):
    print(f"Error: einverted binary at {einverted_bin} is not executable.")
    print("Please run: chmod +x {}".format(einverted_bin))
    sys.exit(1)

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
    
def predict_hybridization(seq1, seq2, temperature=37):
    """
    Predict RNA-RNA interactions using RNAduplex Python bindings.
    
    Args:
        seq1 (str): First RNA sequence
        seq2 (str): Second RNA sequence
        temperature (int): Folding temperature in Celsius (default 37)
    
    Returns:
        str: Formatted output string compatible with parse_rnaduplex_output
    """
    try:
        # Set temperature for ViennaRNA
        RNA.cvar.temperature = temperature
        
        # Calculate duplex using Python bindings
        result = RNA.duplexfold(seq1, seq2)
        
        # Parse structure to get lengths
        structure_parts = result.structure.split('&')
        if len(structure_parts) != 2:
            print(f"Warning: Invalid duplex structure: {result.structure}")
            return ""
        
        len1 = len(structure_parts[0])
        len2 = len(structure_parts[1])
        
        # Calculate indices using the convention we discovered:
        # i = 3' end of seq1 duplex region (1-based)
        # j = 5' start of seq2 duplex region (1-based)
        seq1_start = result.i - len1 + 1
        seq1_end = result.i
        seq2_start = result.j
        seq2_end = result.j + len2 - 1
        
        # Format output to match CLI format for compatibility with parse_rnaduplex_output
        # Format: "structure   seq1_start,seq1_end  :  seq2_start,seq2_end  (energy)"
        output = f"{result.structure}   {seq1_start},{seq1_end}  :  {seq2_start},{seq2_end}  ({result.energy:.2f})"
        
        return output
        
    except ImportError:
        print("Error: RNA module not found. Please install ViennaRNA Python bindings.")
        print("Install with: conda install -c bioconda viennarna")
        return ""
    except Exception as e:
        print(f"Error running RNAduplex Python bindings: {str(e)}")
        return ""

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

def result_writer(output_file, result_queue, num_workers):
    """
    Dedicated process that writes results to file as they arrive from worker processes.
    This runs in a separate process to avoid blocking workers.
    Deduplicates results based on coordinates.
    """
    with open(output_file, 'w') as f:
        # Write header - basic coordinates first, structural details, then effective coords and sequences
        f.write("Chromosome\tStrand\ti_start\ti_end\tj_start\tj_end\t"
                "Score\tRawMatch\tPercMatch\tGaps\t"
                "dG(kcal/mol)\tpercent_paired\tlongest_helix\t"
                "orig_arm_length\teff_arm_length\t"
                "eff_i_start\teff_i_end\teff_j_start\teff_j_end\t"
                "i_seq\tj_seq\tstructure\n")
        
        workers_done = 0
        results_written = 0
        seen_coordinates = set()  # Track unique dsRNA coordinates
        
        while workers_done < num_workers:
            try:
                result = result_queue.get(timeout=1)
                
                if result == "DONE":
                    workers_done += 1
                    continue
                
                # Create unique key based on coordinates
                coord_key = (result['chromosome'], result['strand'], 
                            result['i_start'], result['i_end'],
                            result['j_start'], result['j_end'])
                
                # Skip if we've already seen this dsRNA
                if coord_key in seen_coordinates:
                    continue
                
                seen_coordinates.add(coord_key)
                
                # Write result as TSV line - basic coords, structural details, eff coords, sequences
                f.write(f"{result['chromosome']}\t{result['strand']}\t"
                       f"{result['i_start']}\t{result['i_end']}\t"
                       f"{result['j_start']}\t{result['j_end']}\t"
                       f"{result['score']}\t{result['raw_match']}\t"
                       f"{result['match_perc']}\t{result['gap_numb']}\t"
                       f"{result['energy']}\t{result['percent_paired']}\t{result['longest_helix']}\t"
                       f"{result['orig_arm_length']}\t{result['eff_arm_length']}\t"
                       f"{result['eff_i_start']}\t{result['eff_i_end']}\t"
                       f"{result['eff_j_start']}\t{result['eff_j_end']}\t"
                       f"{result['i_seq']}\t{result['j_seq']}\t"
                       f"{result['structure']}\n")
                
                results_written += 1
                
                # Flush periodically for real-time output
                if result_queue.qsize() < 100:
                    f.flush()
                    
            except Empty:
                continue
            except Exception as e:
                print(f"Error writing result: {e}")
        
        print(f"Writer process finished. Wrote {results_written} results.")

def process_window(i, window_start, window_size, basename, algorithm, args, full_sequence, chromosome, strand, result_queue):
    """Process a genomic window to identify dsRNA structures and stream results to queue
    
    Args:
        i: Start position in the sequence
        window_start: Start position for coordinate calculations
        window_size: Size of the window to process
        basename: Base name for output files
        algorithm: Algorithm to use (einverted)
        args: Command line arguments
        full_sequence: The complete sequence string (already reverse complemented if needed)
        chromosome: Chromosome name
        strand: Strand (+ or -)
        result_queue: Queue for results
    """
    results = []  # Collect results for this window

    if algorithm == "einverted":
        # Extract the window sequence from the provided full sequence
        window_seq = full_sequence[i:i+window_size].upper()
        
        # Check if the sequence is all Ns
        if all(base == 'N' for base in window_seq):
            # Skip this window
            return
        
        if not window_seq:
            return
        
        # Use stdin with einverted - provide sequence directly
        # einverted can accept stdin with -sbegin and -send for coordinates within the stdin sequence
        einverted_cmd = [
            einverted_bin,
            "-sequence", "stdin",  # Read from stdin
            "-sbegin", "1",       # Start at position 1 of stdin sequence
            "-send", str(len(window_seq)),  # End at the length of the window
            "-gap", str(args.gaps),
            "-threshold", str(args.score),
            "-match", str(args.match),
            "-mismatch", str(args.mismatch),
            "-maxrepeat", str(args.max_span),
            "-outfile", "stdout",  # Write to stdout
            "-outseq", "/dev/null",  # Suppress sequence output file
            "-filter"
        ]
        
        # Create FASTA input for stdin
        stdin_input = f">{chromosome}:{i+1}-{i+window_size}\n{window_seq}\n"
        
        process = subprocess.Popen(einverted_cmd, 
                                 stdin=subprocess.PIPE, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, 
                                 text=True)
        stdout, stderr = process.communicate(input=stdin_input)
        
        ein_results = stdout.split("\n")
        results = parse_einverted_results(ein_results, window_start, window_size, basename, args, chromosome, strand)
        
        # Put results in queue for streaming to writer
        for result in results:
            result_queue.put(result)
    
    # Signal this worker is done with this window
    return len(results)
        
def parse_einverted_results(ein_results, window_start, window_size, basename, args, chromosome, strand):
    """Parse results from einverted and return as list of result dictionaries"""
    results = []
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
            score = score_line[2].rstrip(':')  # Remove trailing colon from score
            raw_match = score_line[3]
            matches, total = map(int, raw_match.split('/'))
            match_perc = round((matches / total) * 100, 2)    
            # find gaps one column from last column            
            gap_numb = score_line[-2]
            
            # Calculate the genomic coordinates from einverted output
            # Since we're using stdin, einverted returns coordinates relative to the stdin sequence
            # We need to add window_start to convert back to genomic coordinates
            
            # Calculate the genomic coordinates from einverted output
            # For both strands, we maintain the same coordinate system
            # The only difference is the sequence content (reverse complement for negative strand)
            i_start = int(seq_i_full[0]) + window_start
            i_end = int(seq_i_full[2]) + window_start
            j_start = int(seq_j_full[2]) + window_start
            j_end = int(seq_j_full[0]) + window_start
            
            # Double-check the coordinates are in the correct order
            if i_start > i_end or j_start > j_end or i_start > j_start or i_end > j_end:
                print(f"Warning: Coordinates not in correct order for {i_start} to {j_end}. Sorting...")
                coords = sorted([i_start, i_end, j_start, j_end])
                i_start, i_end, j_start, j_end = coords[0], coords[1], coords[2], coords[3]
            
            # RNA folding and scoring
            # Extract sequences from einverted output
            i_seq = seq_i_full[1].replace("-", "").upper()
            j_seq = ''.join(reversed(seq_j_full[1].replace("-", ""))).upper()
            output = predict_hybridization(i_seq, j_seq, temperature=args.t)
            structure, indices_seq1, indices_seq2, energy = parse_rnaduplex_output(output)
            
            # Skip if we got empty results from RNAduplex
            if not structure:
                j += 5
                continue
                
            # Convert 1-based RNAduplex indices to genomic coordinates
            # RNAduplex returns positions relative to input sequences (1-based)
            # We need to add these to the genomic start positions (0-based adjustment)
            
            # Calculate effective coordinates based on RNAduplex trimming
            # RNAduplex returns 1-based indices for the portions of sequences that form the duplex
            if strand == "-":
                # For reverse strand, the sequences are reverse complemented
                # RNAduplex indices are from the start of the RC sequences
                # We need to adjust for the fact that trimming from the start of RC seq
                # is actually trimming from the end in genomic coordinates
                
                # Get the lengths of the sequences
                i_seq_len = len(i_seq)
                j_seq_len = len(j_seq)
                
                # For reverse strand:
                # - Trimming from start of RC sequence = trimming from end in genomic coords
                # - Trimming from end of RC sequence = trimming from start in genomic coords
                
                # If RNAduplex says use positions 3-10 in a 15bp RC sequence:
                # That means it trimmed 2bp from start and 5bp from end of RC sequence
                # In genomic coords, that's trimming 5bp from start and 2bp from end
                
                eff_i_start = i_start + (i_seq_len - indices_seq1[1])
                eff_i_end = i_end - (indices_seq1[0] - 1)
                eff_j_start = j_start + (j_seq_len - indices_seq2[1])
                eff_j_end = j_end - (indices_seq2[0] - 1)
            else:
                # For forward strand, RNAduplex indices map directly to genomic positions
                # indices are 1-based, so we need to adjust
                eff_i_start = i_start + (indices_seq1[0] - 1)
                eff_i_end = i_start + indices_seq1[1] - 1
                eff_j_start = j_start + (indices_seq2[0] - 1)
                eff_j_end = j_start + indices_seq2[1] - 1
            
            # Debug coordinate conversion if needed
            # print(f"[DEBUG] Coordinate conversion: i_start={i_start}, indices_seq1={indices_seq1} -> eff_i=({eff_i_start}, {eff_i_end})")
            # print(f"[DEBUG] Coordinate conversion: j_start={j_start}, indices_seq2={indices_seq2} -> eff_j=({eff_j_start}, {eff_j_end})")
            
            # Store as tuples for compatibility with existing code
            eff_i = (eff_i_start, eff_i_end)
            eff_j = (eff_j_start, eff_j_end)
            
            # Validate that the effective sequences match the structure length
            i_arm_length = indices_seq1[1] - indices_seq1[0] + 1
            j_arm_length = indices_seq2[1] - indices_seq2[0] + 1
            structure_parts = structure.split('&')
            
            if len(structure_parts) == 2:
                i_structure_length = len(structure_parts[0])
                j_structure_length = len(structure_parts[1])
                
                if i_arm_length != i_structure_length or j_arm_length != j_structure_length:
                    print(f"Warning: Structure length mismatch - i_arm: {i_arm_length} vs {i_structure_length}, j_arm: {j_arm_length} vs {j_structure_length}")
            
            pairs = int(structure.count('(') * 2)
            
            # Calculate percent_paired safely
            try:
                percent_paired = round(float(pairs / (len(structure) - 1)) * 100, 2)
            except (ZeroDivisionError, ValueError):
                percent_paired = 0
            
            # Calculate longest continuous helix
            longest_helix = find_longest_helix(structure)
            
            # Calculate arm lengths
            orig_arm_length = (i_end - i_start + 1)  # Should be same as (j_end - j_start + 1)
            eff_arm_length = (eff_i[1] - eff_i[0] + 1)  # Should be same as (eff_j[1] - eff_j[0] + 1)
            
            if match_perc < args.paired_cutoff:
                print(f"Skipping {i_start} to {j_end} due to low percentage of pairs: {percent_paired}")
                j += 5
                continue
            
            # Use the structure as-is from RNAduplex
            # The structure should remain in the same format regardless of strand
            display_structure = structure
            
            # Create result dictionary instead of writing to file
            result = {
                'chromosome': chromosome,
                'strand': strand,
                'score': score,
                'raw_match': raw_match,
                'match_perc': match_perc,
                'gap_numb': gap_numb,
                'i_start': i_start,
                'i_end': i_end,
                'j_start': j_start,
                'j_end': j_end,
                'eff_i_start': eff_i[0],
                'eff_i_end': eff_i[1],
                'eff_j_start': eff_j[0],
                'eff_j_end': eff_j[1],
                'i_seq': i_seq.replace("T", "U"),  # Convert to RNA for output
                'j_seq': j_seq.replace("T", "U"),  # Convert to RNA for output
                'structure': display_structure,
                'energy': energy,
                'percent_paired': percent_paired,
                'longest_helix': longest_helix,
                'orig_arm_length': orig_arm_length,
                'eff_arm_length': eff_arm_length
            }
            results.append(result)
        except Exception as e:
            print(f"Error processing result block at index {j}: {str(e)}")
        
        # Increment j based on the structure of your einverted output
        j += 5
    
    return results

def find_longest_helix(structure):
    """
    Find the longest stretch of contiguous base pairs in an RNA structure.
    
    Args:
        structure (str): String representing RNA structure (e.g., "(((...)))&(((...)))")
        
    Returns:
        int: Length of the longest contiguous helix (minimum of both arms)
    """
    try:
        # Handle invalid or empty structures
        if not structure or "&" not in structure:
            return 0
            
        # Split structure into both arms
        parts = structure.split("&")
        if len(parts) != 2:
            return 0
            
        left_arm, right_arm = parts[0], parts[1]
        
        # Find longest stretch of opening brackets "(" in left arm
        left_max = 0
        current_left = 0
        for char in left_arm:
            if char == "(":
                current_left += 1
                left_max = max(left_max, current_left)
            else:
                current_left = 0
                
        # Find longest stretch of closing brackets ")" in right arm
        right_max = 0
        current_right = 0
        for char in right_arm:
            if char == ")":
                current_right += 1
                right_max = max(right_max, current_right)
            else:
                current_right = 0
                
        # Return the minimum of the two arms (since a helix requires both sides)
        return min(left_max, right_max)
    except Exception as e:
        print(f"Error calculating longest helix: {e}")
        return 0
            
# Define the process_frame function
def process_frame(frame_start, frame_step_size, end_coordinate, window_size, basename, algorithm, args, full_sequence, chromosome, strand, result_queue, pool):
    for start in range(frame_start, end_coordinate, frame_step_size):
        window_start = start
        end = min(start + window_size, end_coordinate)
        pool.apply_async(process_window, (start, window_start, window_size, basename, algorithm, args, full_sequence, chromosome, strand, result_queue))
        # For debugging, run the process_window function directly
        # process_window(start, start, args.w, basename, args.algorithm, args, full_sequence, chromosome, strand, result_queue)

def main():
    ### Arguments
    parser = argparse.ArgumentParser(
        description='dsRNAscan - A tool for genome-wide prediction of double-stranded RNA structures',
        epilog='Version: {} | Copyright (C) 2024 Bass Lab'.format(__version__)
    )
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
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
    parser.add_argument('--chr', type=str, default='header',
                        help='Chromosome name, if chromosome name in header type "header" (default: header)')
    parser.add_argument('--output_label', type=str, default='header',
                        help='Label for output files and results (default: use sequence header)')
    parser.add_argument('--only_seq', type=str, default=None,
                        help='Only scan this specific sequence')
    parser.add_argument('-c', '--cpus', type=int, default=4,
                        help='Number of cpus to use; Default = 4')
    parser.add_argument('--clean', action='store_false', default=True,
                    help='Clean up temporary files after processing')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: dsrnascan_YYYYMMDD_HHMMSS)')
    
    args = parser.parse_args()
    
    # Validate command line arguments
    if args.w <= 0:
        parser.error("Window size must be greater than 0")
    if args.step <= 0:
        parser.error("Step size must be greater than 0")
    if args.step > args.w:
        print("Warning: Step size is larger than window size. This may cause gaps in coverage.")
    if args.min <= 0:
        parser.error("Minimum inverted repeat length must be greater than 0")
    if args.max < args.min:
        parser.error("Maximum inverted repeat length must be greater than or equal to minimum length")
    if args.cpus <= 0:
        parser.error("Number of CPUs must be greater than 0")
    if args.paired_cutoff < 0 or args.paired_cutoff > 100:
        parser.error("Paired cutoff must be between 0 and 100")
    if args.start < 0:
        parser.error("Start coordinate must be non-negative")
    if args.end < 0:
        parser.error("End coordinate must be non-negative")
    if args.end != 0 and args.end <= args.start:
        parser.error("End coordinate must be greater than start coordinate")
        
    # Check if input file exists and is readable
    if not os.path.exists(args.filename):
        parser.error(f"Input file '{args.filename}' does not exist")
    if not os.access(args.filename, os.R_OK):
        parser.error(f"Input file '{args.filename}' is not readable")
    if os.path.getsize(args.filename) == 0:
        parser.error(f"Input file '{args.filename}' is empty")
        
    # Try to open the file to ensure it's a valid FASTA
    try:
        with smart_open(args.filename) as test_file:
            first_record = next(SeqIO.parse(test_file, "fasta"), None)
            if first_record is None:
                parser.error(f"Input file '{args.filename}' does not appear to be a valid FASTA file")
    except Exception as e:
        parser.error(f"Error reading input file '{args.filename}': {str(e)}")
    
    # Create output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        # Create timestamped directory name
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = f"dsrnascan_{timestamp}"
    
    # Create the directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")
    
    # Create a multiprocessing Manager for the result queue
    manager = multiprocessing.Manager()
    result_queue = manager.Queue(maxsize=10000)  # Buffer up to 10k results
    
    chromosome = args.chr
    end_coordinate = int(args.end)
    fasta_file = args.filename
    cpu_count = args.cpus
    step_size = args.step
    sequence_count = 0

    try:
        with smart_open(args.filename) as f:
            # Create a pool of workers for multiprocessing, starting with 2 workers
            # pool = multiprocessing.Pool(cpu_count)
            
            # Process each sequence
            tasks = []

            for cur_record in SeqIO.parse(f, "fasta"): 
                # Skip sequences if only_seq is specified and this isn't it
                if hasattr(args, 'only_seq') and args.only_seq and cur_record.name != args.only_seq:
                    continue
                    
                sequence_count += 1
                # Correct for starting coordinate
                print(f"Processing sequence: {cur_record.name}")
                
                # Validate sequence
                if len(cur_record.seq) == 0:
                    print(f"Warning: Sequence {cur_record.name} is empty, skipping...")
                    continue
                    
                if not cur_record.seq:
                    print(f"Warning: No sequence data for {cur_record.name}, skipping...")
                    continue   
                # Print the sequence length
                #print(f"Sequence length: {len(cur_record.seq)}")
                
                # Convert to RNA uppercase
                #cur_record.seq = cur_record.seq.transcribe().upper()
                
                # Check if chromosome is 'header' - use new parameter names with fallback
                if hasattr(args, 'output_label'):
                    if args.output_label == "header":
                        chromosome = cur_record.name
                    else:
                        chromosome = args.output_label
                elif args.chr == "header":
                    chromosome = cur_record.name
                else:
                    chromosome = args.chr

                # Determine strand and set up basename
                strand = "-" if args.reverse else "+"
                # Get base filename without extension(s)
                base_filename = args.filename
                if base_filename.endswith('.gz'):
                    base_filename = base_filename[:-3]  # Remove .gz
                if base_filename.endswith('.fa') or base_filename.endswith('.fasta'):
                    base_filename = os.path.splitext(base_filename)[0]
                
                # Prepare the sequence in memory (reverse complement if needed)
                if args.reverse:
                    # For reverse strand, reverse complement the entire sequence once
                    full_sequence = str(cur_record.seq.reverse_complement().upper())
                else:
                    # For forward strand, just use the sequence as-is
                    full_sequence = str(cur_record.seq.upper())
                
                # Set up basename for output files
                basename = f"{base_filename}.{chromosome}.{'reverse' if args.reverse else 'forward'}_win{args.w}_step{args.step}_start{args.start}_score{args.score}"

                # Result files are now written directly via streaming (merged_results.txt)
                
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
                    # For single window, process directly and write results
                    results = process_window(args.start, args.start, seq_length, basename, args.algorithm, 
                                args, full_sequence, chromosome, strand, result_queue)
                    
                    # Write results directly for single window
                    merged_filename = os.path.join(output_dir, f"{os.path.basename(basename)}_merged_results.txt")
                    with open(merged_filename, 'w') as f:
                        # Write header - basic coordinates first, structural details, then effective coords and sequences
                        f.write("Chromosome\tStrand\ti_start\ti_end\tj_start\tj_end\t"
                               "Score\tRawMatch\tPercMatch\tGaps\t"
                               "dG(kcal/mol)\tpercent_paired\tlongest_helix\t"
                               "orig_arm_length\teff_arm_length\t"
                               "eff_i_start\teff_i_end\teff_j_start\teff_j_end\t"
                               "i_seq\tj_seq\tstructure\n")
                        
                        while not result_queue.empty():
                            result = result_queue.get()
                            # Write result - basic coords, structural details, eff coords, sequences
                            f.write(f"{result['chromosome']}\t{result['strand']}\t"
                                   f"{result['i_start']}\t{result['i_end']}\t"
                                   f"{result['j_start']}\t{result['j_end']}\t"
                                   f"{result['score']}\t{result['raw_match']}\t"
                                   f"{result['match_perc']}\t{result['gap_numb']}\t"
                                   f"{result['energy']}\t{result['percent_paired']}\t{result['longest_helix']}\t"
                                   f"{result['orig_arm_length']}\t{result['eff_arm_length']}\t"
                                   f"{result['eff_i_start']}\t{result['eff_i_end']}\t"
                                   f"{result['eff_j_start']}\t{result['eff_j_end']}\t"
                                   f"{result['i_seq']}\t{result['j_seq']}\t"
                                   f"{result['structure']}\n")
                else:
                    # Normal processing for longer sequences
                    print(f"Scanning {cur_record.name} from {args.start} to {end_coordinate} with window size {args.w} and step size {args.step}")
                    
                    # Set up output file
                    merged_filename = os.path.join(output_dir, f"{os.path.basename(basename)}_merged_results.txt")
                    
                    # Start the writer process
                    writer_proc = multiprocessing.Process(target=result_writer, 
                                                        args=(merged_filename, result_queue, cpu_count))
                    writer_proc.start()
                    
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
                                            args, full_sequence, chromosome, strand, result_queue)))
                    # Close the pool and wait for all workers to finish
                    pool.close()
                    pool.join()
                    
                    # Signal writer that all workers are done
                    for _ in range(cpu_count):
                        result_queue.put("DONE")
                    
                    # Wait for writer to finish
                    writer_proc.join()

                    # Results are already written by the writer process
                    print(f"\nResults saved to: {merged_filename}")

                    # If we're only processing one specific sequence, stop after finding it
                    if hasattr(args, 'only_seq') and args.only_seq:
                        break

            # Now generate the BP file
            try:
                # Use the function from the previous script to generate BP file
                bp_filename = os.path.join(output_dir, f"{os.path.basename(basename)}.dsRNApredictions.bp")
                generate_bp_file(merged_filename, bp_filename)
            except NameError:
                print("BP file generation function not defined. Please add the generate_bp_file function to your script.")
            
                        
            # Print file names and paths
            print(f"\nResults written to {merged_filename}")
            if os.path.exists(bp_filename):
                print(f"Base Pair predictions written to {bp_filename}")
            
            # Check if results file is empty or has only headers
            try:
                results_df = pd.read_csv(merged_filename, sep="\t")
                if len(results_df) == 0:
                    print("\nNo dsRNA structures were found. Try adjusting your search parameters.")
            except Exception:
                pass
                
            # Check if any sequences were processed
            if sequence_count == 0:
                print("\nError: No valid sequences found in the input FASTA file.")
                sys.exit(1)
                
    except FileNotFoundError:
        print(f"Error: Could not open file '{args.filename}'. File not found.")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when trying to read '{args.filename}'.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file '{args.filename}': {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
            
# Run the main function
if __name__ == "__main__":
    main()
