#!/usr/bin/env python3
"""
dsRNAscan Results Browser - Simple Forna Viewer

Quick visualization of dsRNA structures from dsRNAscan output.
Supports optional RNA editing site annotation from BED or GFF3 files.
"""

import os
import sys
import json
import pandas as pd
import numpy as np
import socketserver
import webbrowser
import argparse
from pathlib import Path
import tempfile
import shutil
from http.server import SimpleHTTPRequestHandler
from collections import defaultdict
import bisect

class CORSHTTPRequestHandler(SimpleHTTPRequestHandler):
    """HTTP request handler with CORS support"""
    
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        super().end_headers()
    
    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()
    
    def handle_error(self, request, client_address):
        """Override to suppress broken pipe errors"""
        import errno
        exc_type, exc_value = sys.exc_info()[:2]
        if exc_type == BrokenPipeError or (exc_type == OSError and exc_value.errno == errno.EPIPE):
            # Ignore broken pipe errors
            pass
        else:
            super().handle_error(request, client_address)

class OptimizedEditingSites:
    """Optimized storage and lookup for editing sites using binary search"""
    
    def __init__(self):
        # Store as sorted arrays per chromosome/strand for binary search
        # Using numpy arrays for memory efficiency
        self.sites = defaultdict(lambda: defaultdict(lambda: {'positions': [], 'frequencies': []}))
        self.finalized = False
    
    def add_site(self, chrom, strand, position, frequency):
        """Add a site while maintaining sorted order"""
        if self.finalized:
            raise RuntimeError("Cannot add sites after finalization")
        
        sites = self.sites[chrom][strand]
        # Use bisect to maintain sorted order during insertion
        idx = bisect.bisect_left(sites['positions'], position)
        sites['positions'].insert(idx, position)
        sites['frequencies'].insert(idx, frequency)
    
    def finalize(self):
        """Convert lists to numpy arrays for memory efficiency"""
        print("Optimizing editing site storage...")
        for chrom in self.sites:
            for strand in self.sites[chrom]:
                data = self.sites[chrom][strand]
                if data['positions']:
                    # Use uint32 for positions (saves memory for large files)
                    data['positions'] = np.array(data['positions'], dtype=np.uint32)
                    # Use float16 for frequencies (sufficient precision, saves memory)
                    data['frequencies'] = np.array(data['frequencies'], dtype=np.float16)
        self.finalized = True
    
    def find_in_range(self, chrom, strand, start, end):
        """Binary search for sites in range - O(log n) instead of O(n)"""
        if chrom not in self.sites or strand not in self.sites[chrom]:
            return []
        
        sites = self.sites[chrom][strand]
        positions = sites['positions']
        
        if len(positions) == 0:
            return []
        
        # Binary search for range
        if self.finalized:
            # NumPy searchsorted for arrays
            start_idx = np.searchsorted(positions, start, side='left')
            end_idx = np.searchsorted(positions, end, side='right')
        else:
            # bisect for lists
            start_idx = bisect.bisect_left(positions, start)
            end_idx = bisect.bisect_right(positions, end)
        
        return [(int(positions[i]), float(sites['frequencies'][i])) 
                for i in range(start_idx, end_idx)]
    
    def get_summary(self):
        """Get summary statistics"""
        total_sites = 0
        chroms = set()
        for chrom in self.sites:
            chroms.add(chrom)
            for strand in self.sites[chrom]:
                total_sites += len(self.sites[chrom][strand]['positions'])
        return {'total_sites': total_sites, 'chromosomes': len(chroms)}

def parse_editing_file_optimized(editing_file, needed_chromosomes=None):
    """
    Parse BED or GFF3 file with editing sites using optimized storage.
    
    Args:
        editing_file: Path to BED or GFF3 file
        needed_chromosomes: Set of chromosomes to load (None = load all)
    
    Returns:
        OptimizedEditingSites object
    """
    editing_sites = OptimizedEditingSites()
    
    # First, count total lines for progress reporting
    print(f"Counting entries in {editing_file}...")
    total_lines = 0
    with open(editing_file, 'r') as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                total_lines += 1
    
    if total_lines > 10000:
        print(f"Large file detected: {total_lines:,} entries. Processing...")
    
    try:
        # Detect file format
        is_gff3 = False
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('##gff-version'):
                    is_gff3 = True
                    break
                elif not line.startswith('#') and line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) == 9 and parts[2] not in ['', '.']:
                        is_gff3 = True
                    break
        
        # Parse based on format
        processed = 0
        skipped = 0
        
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                chrom = parts[0]
                
                # Skip chromosomes we don't need (lazy loading)
                if needed_chromosomes and chrom not in needed_chromosomes:
                    skipped += 1
                    processed += 1
                    continue
                
                if is_gff3:
                    # GFF3 format parsing
                    if len(parts) < 9:
                        continue
                    
                    feature_type = parts[2]
                    
                    # Filter for editing-related features
                    if not any(term in feature_type.lower() for term in ['edit', 'modification', 'variant', 'snp', 'snv']):
                        attributes = parts[8]
                        if not any(term in attributes.lower() for term in ['edit', 'adar', 'apobec', 'a-to-i', 'c-to-u']):
                            skipped += 1
                            processed += 1
                            continue
                    
                    position = int(parts[3])  # GFF3 is 1-based
                    strand = parts[6]
                    
                    # Get frequency from score or attributes
                    frequency = 1.0
                    if parts[5] != '.':
                        try:
                            score = float(parts[5])
                            if score <= 1.0:
                                frequency = score
                            elif score <= 100:
                                frequency = score / 100.0
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                    
                    # Parse attributes for frequency
                    attributes = parts[8]
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            key = key.strip().lower()
                            if key in ['frequency', 'freq', 'confidence', 'score', 'editing_level']:
                                try:
                                    freq_val = float(value.strip('%'))
                                    if freq_val > 1.0:
                                        frequency = freq_val / 100.0
                                    else:
                                        frequency = freq_val
                                    break
                                except:
                                    pass
                
                else:
                    # BED format parsing
                    if len(parts) < 6:
                        continue
                    
                    position = int(parts[1]) + 1  # BED is 0-based, convert to 1-based
                    strand = parts[5]
                    
                    frequency = 1.0
                    if len(parts) > 4:
                        try:
                            score = float(parts[4])
                            if score <= 1.0:
                                frequency = score
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                
                # Add to optimized structure
                editing_sites.add_site(chrom, strand, position, frequency)
                
                processed += 1
                # Show progress for large files
                if total_lines > 10000 and processed % 10000 == 0:
                    pct = (processed / total_lines) * 100
                    print(f"  Processed {processed:,} / {total_lines:,} ({pct:.1f}%)", end='\r')
        
        # Clear progress line
        if total_lines > 10000:
            print(f"\n  Loaded {processed - skipped:,} sites, skipped {skipped:,}")
        
        # Finalize for memory efficiency
        editing_sites.finalize()
        
        # Report summary
        summary = editing_sites.get_summary()
        format_name = "GFF3" if is_gff3 else "BED"
        print(f"Loaded {summary['total_sites']:,} editing sites from {summary['chromosomes']} chromosomes ({format_name} format)")
        
        return editing_sites
    
    except Exception as e:
        print(f"Warning: Could not parse editing file: {e}")
        return OptimizedEditingSites()

def parse_editing_file(editing_file):
    """
    Parse BED or GFF3 file with editing sites.
    
    BED format: chr start end name score strand [frequency]
    GFF3 format: chr source type start end score strand phase attributes
    
    Returns dict: {chr: {strand: [(position, frequency)]}}
    """
    editing_sites = defaultdict(lambda: defaultdict(list))
    
    try:
        # Detect file format
        is_gff3 = False
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('##gff-version'):
                    is_gff3 = True
                    break
                elif not line.startswith('#') and line.strip():
                    # Check if it looks like GFF3 (9 tab-separated fields)
                    parts = line.strip().split('\t')
                    if len(parts) == 9 and parts[2] not in ['', '.']:
                        is_gff3 = True
                    break
        
        # Parse based on format
        with open(editing_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                
                if is_gff3:
                    # GFF3 format parsing
                    if len(parts) < 9:
                        continue
                    
                    chrom = parts[0]
                    feature_type = parts[2]
                    
                    # Filter for editing-related features
                    if not any(term in feature_type.lower() for term in ['edit', 'modification', 'variant', 'snp', 'snv']):
                        # Also check attributes for editing-related terms
                        attributes = parts[8]
                        if not any(term in attributes.lower() for term in ['edit', 'adar', 'apobec', 'a-to-i', 'c-to-u']):
                            continue
                    
                    # GFF3 is 1-based
                    position = int(parts[3])
                    strand = parts[6]
                    
                    # Try to get frequency from score or attributes
                    frequency = 1.0  # Default frequency
                    
                    # Check score field
                    if parts[5] != '.':
                        try:
                            score = float(parts[5])
                            if score <= 1.0:
                                frequency = score
                            elif score <= 100:
                                frequency = score / 100.0
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                    
                    # Parse attributes for frequency/confidence
                    attributes = parts[8]
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            key = key.strip().lower()
                            if key in ['frequency', 'freq', 'confidence', 'score', 'editing_level']:
                                try:
                                    freq_val = float(value.strip('%'))
                                    if freq_val > 1.0:
                                        frequency = freq_val / 100.0
                                    else:
                                        frequency = freq_val
                                    break
                                except:
                                    pass
                    
                else:
                    # BED format parsing
                    if len(parts) < 6:
                        continue
                    
                    chrom = parts[0]
                    # BED is 0-based, convert to 1-based for genomic coords
                    position = int(parts[1]) + 1  
                    strand = parts[5]
                    
                    # Try to get frequency from score or additional column
                    frequency = 1.0  # Default frequency
                    if len(parts) > 4:
                        try:
                            score = float(parts[4])
                            if score <= 1.0:
                                frequency = score
                            else:
                                frequency = score / 1000.0
                        except:
                            pass
                
                # Store as 1-based genomic position
                editing_sites[chrom][strand].append((position, frequency))
        
        # Sort positions for each chromosome/strand
        for chrom in editing_sites:
            for strand in editing_sites[chrom]:
                editing_sites[chrom][strand].sort(key=lambda x: x[0])
        
        # Report what was loaded
        total_sites = sum(len(sites) for chrom in editing_sites for sites in editing_sites[chrom].values())
        format_name = "GFF3" if is_gff3 else "BED"
        print(f"Loaded {total_sites} editing sites from {format_name} file")
        
        return editing_sites
    
    except Exception as e:
        print(f"Warning: Could not parse editing file: {e}")
        return {}

def map_editing_to_dsrna_optimized(dsrna_row, editing_sites):
    """
    Map editing sites to dsRNA structure positions using optimized binary search.
    
    For forward strand: positions map directly
    For reverse strand: need to reverse the mapping
    
    Returns list of [position_in_structure, frequency] pairs
    """
    if not editing_sites:
        return []
    
    chrom = str(dsrna_row['Chromosome'])
    strand = str(dsrna_row['Strand'])
    
    # Get the effective coordinates (trimmed regions)
    if 'eff_i_start' in dsrna_row and pd.notna(dsrna_row['eff_i_start']):
        i_start = int(dsrna_row['eff_i_start'])
        i_end = int(dsrna_row['eff_i_end'])
        j_start = int(dsrna_row['eff_j_start'])
        j_end = int(dsrna_row['eff_j_end'])
    else:
        # Fall back to original coordinates
        i_start = int(dsrna_row['i_start'])
        i_end = int(dsrna_row['i_end'])
        j_start = int(dsrna_row['j_start'])
        j_end = int(dsrna_row['j_end'])
    
    # Get sequences to determine lengths
    i_seq = str(dsrna_row['i_seq']) if pd.notna(dsrna_row['i_seq']) else ''
    j_seq = str(dsrna_row['j_seq']) if pd.notna(dsrna_row['j_seq']) else ''
    i_length = len(i_seq)
    j_length = len(j_seq)
    
    structure_edits = []
    
    # Use optimized binary search to find sites in i-arm range
    if strand == '+':
        # Forward strand
        i_sites = editing_sites.find_in_range(chrom, strand, i_start, i_end)
        for edit_pos, frequency in i_sites:
            structure_pos = edit_pos - i_start
            if 0 <= structure_pos < i_length:
                structure_edits.append([int(structure_pos), float(frequency)])
        
        # Find sites in j-arm range
        j_sites = editing_sites.find_in_range(chrom, strand, j_start, j_end)
        for edit_pos, frequency in j_sites:
            structure_pos = i_length + (edit_pos - j_start)
            if i_length <= structure_pos < (i_length + j_length):
                structure_edits.append([int(structure_pos), float(frequency)])
    
    else:  # strand == '-'
        # Reverse strand: coordinates go from high to low (3' to 5')
        # i_start is at the 3' end, i_end is at the 5' end (i_end < i_start)
        i_sites = editing_sites.find_in_range(chrom, strand, i_end, i_start)
        for edit_pos, frequency in i_sites:
            pos_from_3prime = edit_pos - i_end
            structure_pos = i_length - 1 - pos_from_3prime
            if 0 <= structure_pos < i_length:
                structure_edits.append([int(structure_pos), float(frequency)])
        
        # j-arm (j_start < j_end for minus)
        j_sites = editing_sites.find_in_range(chrom, strand, j_start, j_end)
        for edit_pos, frequency in j_sites:
            pos_from_3prime = edit_pos - j_start
            structure_pos = i_length + (j_length - 1 - pos_from_3prime)
            if i_length <= structure_pos < (i_length + j_length):
                structure_edits.append([int(structure_pos), float(frequency)])
    
    return structure_edits

def map_editing_to_dsrna(dsrna_row, editing_sites):
    """
    Map editing sites to dsRNA structure positions (legacy version for compatibility).
    
    For forward strand: positions map directly
    For reverse strand: need to reverse the mapping
    
    Returns list of [position_in_structure, frequency] pairs
    """
    # Check if we have the optimized structure
    if isinstance(editing_sites, OptimizedEditingSites):
        return map_editing_to_dsrna_optimized(dsrna_row, editing_sites)
    
    # Legacy code for backward compatibility
    if not editing_sites:
        return []
    
    chrom = str(dsrna_row['Chromosome'])
    strand = str(dsrna_row['Strand'])
    
    # Get editing sites for this chromosome/strand
    if chrom not in editing_sites or strand not in editing_sites[chrom]:
        return []
    
    chrom_edits = editing_sites[chrom][strand]
    
    # Get the effective coordinates (trimmed regions)
    if 'eff_i_start' in dsrna_row and pd.notna(dsrna_row['eff_i_start']):
        i_start = int(dsrna_row['eff_i_start'])
        i_end = int(dsrna_row['eff_i_end'])
        j_start = int(dsrna_row['eff_j_start'])
        j_end = int(dsrna_row['eff_j_end'])
    else:
        # Fall back to original coordinates
        i_start = int(dsrna_row['i_start'])
        i_end = int(dsrna_row['i_end'])
        j_start = int(dsrna_row['j_start'])
        j_end = int(dsrna_row['j_end'])
    
    structure_edits = []
    
    # Get sequences to determine lengths
    i_seq = str(dsrna_row['i_seq']) if pd.notna(dsrna_row['i_seq']) else ''
    j_seq = str(dsrna_row['j_seq']) if pd.notna(dsrna_row['j_seq']) else ''
    i_length = len(i_seq)
    j_length = len(j_seq)
    
    # Map editing sites to structure positions
    for edit_pos, frequency in chrom_edits:
        structure_pos = None
        
        if strand == '+':
            # Forward strand: direct mapping
            # Check if edit is in i-arm
            if i_start <= edit_pos <= i_end:
                # Position in i-arm (0-based for structure)
                structure_pos = edit_pos - i_start
            # Check if edit is in j-arm
            elif j_start <= edit_pos <= j_end:
                # Position in j-arm, offset by i-arm length
                structure_pos = i_length + (edit_pos - j_start)
        
        else:  # strand == '-'
            # Reverse strand: coordinates go from high to low (3' to 5')
            # Need to reverse the mapping
            
            # For minus strand, genomic coords are reversed:
            # i_start is actually at the 3' end of i-arm
            # j_end is actually at the 5' end of j-arm
            
            # Check if edit is in i-arm (remember: i_end < i_start for minus)
            if i_end <= edit_pos <= i_start:
                # Position from 3' end of i-arm, but we need position from 5' end
                # So we reverse it
                pos_from_3prime = edit_pos - i_end
                structure_pos = i_length - 1 - pos_from_3prime
            
            # Check if edit is in j-arm (j_start < j_end for minus) 
            elif j_start <= edit_pos <= j_end:
                # Position from 3' end of j-arm
                pos_from_3prime = edit_pos - j_start
                # In structure, j-arm comes after i-arm
                # But for minus strand, we need to reverse within j-arm
                structure_pos = i_length + (j_length - 1 - pos_from_3prime)
        
        # Add to list if within structure bounds
        if structure_pos is not None and 0 <= structure_pos < (i_length + j_length):
            structure_edits.append([int(structure_pos), frequency])
    
    return structure_edits

def process_dsrnascan_results(results_file, editing_sites=None):
    """
    Process dsRNAscan results file and convert to JSON format for web display
    """
    try:
        # Read the results file
        df = pd.read_csv(results_file, sep='\t')
        
        # Filter out empty rows
        df = df.dropna(subset=['Chromosome'])
        
        if df.empty:
            print(f"Warning: No valid results found in {results_file}")
            return []
        
        # Convert to list of dictionaries for JSON serialization
        results = []
        for idx, row in df.iterrows():
            result = {
                'id': f'dsRNA_{idx}',
                'chromosome': str(row['Chromosome']),
                'strand': str(row['Strand']),
                'score': float(row['Score']) if pd.notna(row['Score']) else 0,
                'i_start': int(row['i_start']),
                'i_end': int(row['i_end']),
                'j_start': int(row['j_start']),
                'j_end': int(row['j_end']),
                'i_seq': str(row['i_seq']) if pd.notna(row['i_seq']) else '',
                'j_seq': str(row['j_seq']) if pd.notna(row['j_seq']) else '',
                'structure': str(row['structure']) if pd.notna(row['structure']) else '',
                'dG': float(row['dG(kcal/mol)']) if pd.notna(row['dG(kcal/mol)']) else 0,
                'percent_paired': float(row['percent_paired']) if pd.notna(row['percent_paired']) else 0
            }
            
            # Add optional fields if they exist
            if 'longest_helix' in row and pd.notna(row['longest_helix']):
                result['longest_helix'] = int(row['longest_helix'])
            
            # Add effective coordinates if available
            if 'eff_i_start' in row and pd.notna(row['eff_i_start']):
                result['eff_i_start'] = int(row['eff_i_start'])
                result['eff_i_end'] = int(row['eff_i_end'])
                result['eff_j_start'] = int(row['eff_j_start'])
                result['eff_j_end'] = int(row['eff_j_end'])
            
            # Map editing sites if provided
            if editing_sites:
                result['editing_sites'] = map_editing_to_dsrna(row, editing_sites)
            
            results.append(result)
        
        return results
    
    except Exception as e:
        print(f"Error processing {results_file}: {e}")
        return []

def create_html_page(results_data, output_dir, has_editing=False):
    """
    Create a simple HTML page for viewing dsRNA structures with Forna
    """
    html_content = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>dsRNA Structure Viewer''' + (' with Editing Sites' if has_editing else '') + '''</title>
    
    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    
    <!-- D3 and Fornac for structure visualization -->
    <script src="https://d3js.org/d3.v3.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/ViennaRNA/fornac@master/dist/fornac.css">
    
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        
        .header {
            text-align: center;
            margin-bottom: 30px;
        }
        
        .controls {
            text-align: center;
            margin-bottom: 20px;
        }
        
        select {
            padding: 8px;
            font-size: 16px;
            margin: 0 10px;
        }
        
        /* Style for options with editing sites */
        select option[style*="background-color"] {
            background-color: #D5F4E6 !important;
            font-weight: 500;
        }
        
        /* Additional visual cue with emoji */
        .has-edits::before {
            content: "🟢 ";
        }
        
        .info {
            background: white;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        #structure-viewer {
            background: white;
            height: 600px;
            border: 1px solid #ddd;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        .metric {
            display: inline-block;
            margin-right: 30px;
            font-size: 16px;
        }
        
        .metric-label {
            font-weight: bold;
            color: #666;
        }
        
        .metric-value {
            color: #2c3e50;
            font-family: monospace;
            font-size: 18px;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>dsRNA Structure Viewer''' + (' with RNA Editing Sites' if has_editing else '') + '''</h1>
        <p>Visualize predicted double-stranded RNA structures''' + (' with editing annotations' if has_editing else '') + '''</p>
    </div>
    
    <div class="controls">
        <label>Select structure: </label>
        <select id="structure-select">
            <option value="">Choose a structure...</option>
        </select>
        <div id="edit-summary" style="margin-top: 10px; font-size: 14px; color: #27AE60;"></div>
    </div>
    
    <div class="info" id="info-panel" style="display: none;">
        <div class="metric">
            <span class="metric-label">Location:</span>
            <span class="metric-value" id="location"></span>
        </div>
        <div class="metric">
            <span class="metric-label">Free Energy:</span>
            <span class="metric-value" id="energy"></span>
        </div>
        <div class="metric">
            <span class="metric-label">Base Pairs:</span>
            <span class="metric-value" id="pairs"></span>
        </div>
        <div class="metric">
            <span class="metric-label">Score:</span>
            <span class="metric-value" id="score"></span>
        </div>
    </div>
    
    <div id="structure-viewer"></div>
    
    <script>
        // Load the results data
        const resultsData = ''' + json.dumps(results_data, ensure_ascii=True) + ''';
        
        // Function to annotate editing sites on the structure
        function annotateEditingSites(containerId, editingSites) {
            const svg = d3.select(`#${containerId} svg`);
            const allCircles = svg.selectAll('circle');
            
            // Get the most common radius (likely nucleotides)
            const radiusCounts = {};
            allCircles.each(function() {
                const r = parseFloat(d3.select(this).attr('r'));
                radiusCounts[r] = (radiusCounts[r] || 0) + 1;
            });
            
            let mostCommonRadius = 0;
            let maxCount = 0;
            for (const [radius, count] of Object.entries(radiusCounts)) {
                if (count > maxCount) {
                    maxCount = count;
                    mostCommonRadius = parseFloat(radius);
                }
            }
            
            // Select circles with the most common radius (nucleotides)
            const nucleotideCircles = [];
            allCircles.each(function() {
                const circle = d3.select(this);
                const r = parseFloat(circle.attr('r'));
                if (Math.abs(r - mostCommonRadius) < 0.1) {
                    nucleotideCircles.push(this);
                }
            });
            
            // Annotate editing sites
            editingSites.forEach(site => {
                let pos, frequency;
                
                if (typeof site === 'number') {
                    pos = site;
                    frequency = 1.0;
                } else if (Array.isArray(site)) {
                    pos = site[0];
                    frequency = site[1];
                } else {
                    return;
                }
                
                if (pos < nucleotideCircles.length) {
                    const circle = d3.select(nucleotideCircles[pos]);
                    
                    // Calculate color based on frequency
                    let strokeColor;
                    let strokeWidth = '3px';
                    
                    if (frequency >= 0.8) {
                        strokeColor = '#0B5345';  // Dark green
                    } else if (frequency >= 0.5) {
                        strokeColor = '#148F77';  // Medium-dark green
                    } else if (frequency >= 0.3) {
                        strokeColor = '#27AE60';  // Medium green
                    } else if (frequency >= 0.1) {
                        strokeColor = '#52BE80';  // Light-medium green
                    } else {
                        strokeColor = '#ABEBC6';  // Very light green
                        strokeWidth = '2px';
                    }
                    
                    circle
                        .style('stroke', strokeColor)
                        .style('stroke-width', strokeWidth)
                        .style('stroke-opacity', 1.0)
                        .classed('editing-site', true);
                    
                    // Add tooltip
                    let title = circle.select('title');
                    if (title.empty()) {
                        title = circle.append('title');
                    }
                    const freqText = ` (${(frequency * 100).toFixed(1)}%)`;
                    title.text(`Position ${pos + 1}: RNA editing site${freqText}`);
                }
            });
        }
        
        // Initialize when page loads
        document.addEventListener('DOMContentLoaded', function() {
            const select = document.getElementById('structure-select');
            
            // Count dsRNAs with editing sites
            let editedCount = 0;
            let totalEdits = 0;
            
            // Populate dropdown
            resultsData.forEach((result, index) => {
                const hasEdits = result.editing_sites && result.editing_sites.length > 0;
                const editCount = hasEdits ? result.editing_sites.length : 0;
                const editInfo = hasEdits ? ` [${editCount} edits]` : '';
                
                if (hasEdits) {
                    editedCount++;
                    totalEdits += editCount;
                }
                
                // Create option with styling for edited entries
                let optionHTML = '<option value="' + index + '"';
                if (hasEdits) {
                    // Add green background for entries with editing sites
                    optionHTML += ' style="background-color: #D5F4E6; font-weight: 500;"';
                }
                optionHTML += '>' + result.chromosome + ':' + result.i_start + '-' + result.j_end + 
                    ' (' + result.strand + ', dG=' + result.dG.toFixed(1) + editInfo + ')</option>';
                
                select.innerHTML += optionHTML;
            });
            
            // Show summary if there are editing sites
            if (editedCount > 0) {
                const summary = document.getElementById('edit-summary');
                summary.innerHTML = `<strong style="color: #27AE60;">✓</strong> ${editedCount} of ${resultsData.length} dsRNAs have editing sites (${totalEdits} total sites)`;
                summary.style.display = 'block';
            }
            
            // Handle selection
            select.addEventListener('change', function() {
                if (this.value === '') {
                    document.getElementById('info-panel').style.display = 'none';
                    document.getElementById('structure-viewer').innerHTML = '';
                    return;
                }
                
                const data = resultsData[parseInt(this.value)];
                displayInfo(data);
                visualizeStructure(data);
            });
        });
        
        function displayInfo(data) {
            document.getElementById('info-panel').style.display = 'block';
            document.getElementById('location').textContent = `${data.chromosome}:${data.i_start}-${data.j_end} (${data.strand})`;
            document.getElementById('energy').textContent = `${data.dG.toFixed(2)} kcal/mol`;
            document.getElementById('pairs').textContent = `${data.percent_paired.toFixed(1)}%`;
            document.getElementById('score').textContent = data.score;
            
            // Add simple strand indicator
            const infoPanel = document.getElementById('info-panel');
            const existingNote = infoPanel.querySelector('.strand-note');
            if (existingNote) {
                existingNote.remove();
            }
            
            // Show detailed structure info for RNA biologists
            const note = document.createElement('div');
            note.className = 'strand-note';
            note.style.cssText = 'margin-top: 10px; padding: 8px; background: #f0f0f0; border-radius: 4px; font-size: 14px; color: #555;';
            const iLength = data.i_seq ? data.i_seq.length : 0;
            const jLength = data.j_seq ? data.j_seq.length : 0;
            
            // Calculate additional metrics
            const totalLength = iLength + jLength;
            const longestHelix = data.longest_helix || 'N/A';
            const basePairs = Math.round(data.percent_paired * totalLength / 100);
            const editCount = data.editing_sites ? data.editing_sites.length : 0;
            
            // Build detailed info
            let infoHTML = '<strong>Structure Details:</strong><br>';
            infoHTML += `• i-arm: positions 1-${iLength} (${iLength} nt)<br>`;
            infoHTML += `• j-arm: positions ${iLength+1}-${totalLength} (${jLength} nt)<br>`;
            infoHTML += `• Total length: ${totalLength} nt<br>`;
            infoHTML += `• Base pairs: ${basePairs} (${data.percent_paired.toFixed(1)}%)<br>`;
            infoHTML += `• Longest helix: ${longestHelix} bp<br>`;
            
            if (editCount > 0) {
                infoHTML += `• <strong style="color: #27AE60;">RNA editing sites: ${editCount}</strong><br>`;
            }
            
            // Add genomic coordinates if available
            if (data.i_start && data.j_end) {
                const genomicSpan = data.j_end - data.i_start + 1;
                infoHTML += `• Genomic span: ${genomicSpan.toLocaleString()} bp`;
            }
            
            note.innerHTML = infoHTML;
            infoPanel.appendChild(note);
        }
        
        function visualizeStructure(data) {
            // Clear previous visualization
            $('#structure-viewer').empty();
            
            // Split the RNAduplex structure on & to get individual arm structures
            const structureParts = data.structure.split('&');
            let struct_i = structureParts[0];
            let struct_j = structureParts[1];
            
            // Get sequences (already effective/trimmed from dsRNAscan)
            let eff_i_seq = data.i_seq;
            let eff_j_seq = data.j_seq;
            
            // Combine sequences and structures for full dsRNA visualization
            let sequence = eff_i_seq + eff_j_seq;
            let structure = struct_i + struct_j;
            
            // Create container for Forna
            const container = d3.select('#structure-viewer')
                .append('div')
                .attr('id', 'forna-container')
                .style('width', '100%')
                .style('height', '600px');
            
            // Initialize Forna
            const rnaViz = new fornac.FornaContainer("#forna-container", {
                applyForce: true,
                allowPanningAndZooming: true,
                initialSize: [800, 580],
                friction: 0.35,
                middleCharge: -30,
                otherCharge: -30
            });
            
            // Add the RNA structure
            const options = {
                structure: structure,
                sequence: sequence,
                name: data.id,
                labelInterval: 10
            };
            
            rnaViz.addRNA(structure, options);
            
            // Apply editing sites after structure loads
            if (data.editing_sites && data.editing_sites.length > 0) {
                setTimeout(() => {
                    annotateEditingSites('forna-container', data.editing_sites);
                }, 1500);
            }
            
            // Center and zoom
            setTimeout(() => {
                rnaViz.centerView();
                rnaViz.zoomToFit();
            }, 1000);
        }
    </script>
</body>
</html>'''
    
    # Write HTML file
    html_path = os.path.join(output_dir, 'index.html')
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    return html_path

def main():
    parser = argparse.ArgumentParser(description='Browse dsRNAscan results with Forna visualization')
    parser.add_argument('output_dir', nargs='?', default='.', 
                       help='dsRNAscan output directory (default: current directory)')
    parser.add_argument('--editing-file', type=str,
                       help='BED or GFF3 file with RNA editing sites')
    parser.add_argument('--large-editing-file', action='store_true',
                       help='Force optimized parsing for large editing files')
    parser.add_argument('--port', type=int, default=8080,
                       help='Port for the web server (default: 8080)')
    parser.add_argument('--no-browser', action='store_true',
                       help='Do not automatically open the web browser')
    
    args = parser.parse_args()
    
    # Find results files in the output directory
    output_path = Path(args.output_dir)
    if not output_path.exists():
        print(f"Error: Directory {args.output_dir} does not exist")
        sys.exit(1)
    
    # Look for merged_results.txt files (the actual dsRNAscan output)
    results_files = list(output_path.glob('*_merged_results.txt'))
        
    if not results_files:
        print(f"Error: No dsRNAscan results files (*_merged_results.txt) found in {args.output_dir}")
        print("Make sure you're in a dsRNAscan output directory (e.g., dsrnascan_YYYYMMDD_HHMMSS/)")
        sys.exit(1)
    
    # First scan dsRNAscan results to find which chromosomes we need
    needed_chromosomes = set()
    for rf in results_files:
        try:
            df = pd.read_csv(rf, sep='\t', usecols=['Chromosome'])
            needed_chromosomes.update(df['Chromosome'].unique())
        except:
            pass
    
    print(f"Found dsRNAs on {len(needed_chromosomes)} chromosomes")
    
    # Parse editing sites if provided (with lazy loading for needed chromosomes only)
    editing_sites = None
    if args.editing_file:
        print(f"Loading editing sites from {args.editing_file}...")
        
        # Check file size or use forced optimization
        file_size = os.path.getsize(args.editing_file)
        use_optimized = args.large_editing_file or file_size > 10_000_000  # 10MB
        
        if use_optimized:
            print("Using optimized parser for large file...")
            editing_sites = parse_editing_file_optimized(args.editing_file, needed_chromosomes)
        else:
            # For smaller files, still use optimized parser but load all chromosomes
            editing_sites = parse_editing_file_optimized(args.editing_file, None)
    
    # Process all results files
    all_results = []
    for rf in results_files:
        print(f"Processing {rf}...")
        results = process_dsrnascan_results(rf, editing_sites)
        all_results.extend(results)
    
    print(f"Loaded {len(all_results)} dsRNA predictions")
    
    if editing_sites:
        # Count how many dsRNAs have editing sites
        with_edits = sum(1 for r in all_results if r.get('editing_sites'))
        total_edits = sum(len(r.get('editing_sites', [])) for r in all_results)
        print(f"Mapped {total_edits} editing sites to {with_edits} dsRNA structures")
        
        # Show performance benefit if using optimized structure
        if isinstance(editing_sites, OptimizedEditingSites):
            summary = editing_sites.get_summary()
            print(f"  (Optimized lookup from {summary['total_sites']:,} total sites)")
    
    if not all_results:
        print("Error: No valid results found in the input files")
        sys.exit(1)
    
    # Create temporary directory for web files
    temp_dir = tempfile.mkdtemp(prefix='dsrna_browser_')
    
    try:
        # Create HTML page
        html_path = create_html_page(all_results, temp_dir, has_editing=bool(editing_sites))
        
        # Change to temp directory
        os.chdir(temp_dir)
        
        # Start web server
        with socketserver.TCPServer(("", args.port), CORSHTTPRequestHandler) as httpd:
            print(f"\nServer running at http://localhost:{args.port}/")
            print(f"Serving files from: {temp_dir}")
            print("\nPress Ctrl+C to stop the server")
            
            # Open browser if requested
            if not args.no_browser:
                webbrowser.open(f'http://localhost:{args.port}/')
            
            try:
                httpd.serve_forever()
            except KeyboardInterrupt:
                print("\nServer stopped.")
    
    finally:
        # Cleanup temp directory
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == '__main__':
    main()