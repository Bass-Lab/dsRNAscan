#!/usr/bin/env python3
"""
Prepare dsRNA subset files from main parquet file
This script creates pre-processed BED/parquet files for common dsRNA subsets
to enable faster loading during overlap analysis.

Usage:
    python prepare_data.py --input /path/to/dsrna_predictions.parquet
    python prepare_data.py --download  # Download from figshare/zenodo
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import pyarrow.parquet as pq
import urllib.request
import shutil

def download_dsrna_data(output_path="data/dsrna_predictions.parquet"):
    """Download the full dsRNA predictions file"""
    # This URL would be updated to point to the actual data repository
    # For now, we'll provide instructions for manual download
    print("=" * 80)
    print("DOWNLOADING dsRNA PREDICTIONS DATA")
    print("=" * 80)
    print("\nThe full dsRNA predictions file (~500MB) can be downloaded from:")
    print("\n  Option 1: Zenodo")
    print("  https://zenodo.org/records/XXXXXXX")  # To be updated with actual DOI
    print("\n  Option 2: Figshare")  
    print("  https://figshare.com/articles/dataset/XXXXXXX")  # To be updated
    print("\n  Option 3: Direct from paper supplementary data")
    print("  See: https://doi.org/10.1101/2025.01.24.634786")
    print("\nOnce downloaded, run:")
    print(f"  python prepare_data.py --input /path/to/downloaded/file.parquet")
    print("=" * 80)
    return False

def create_subsets(input_file, output_dir="data"):
    """Create subset files from main parquet"""
    
    print(f"Loading dsRNA data from {input_file}...")
    df = pd.read_parquet(input_file)
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Total dsRNAs loaded: {len(df):,}")
    
    # Define subset criteria
    subsets = {
        'conserved': {
            'description': 'Conserved dsRNAs (PhastCons > 0.5)',
            'filter': lambda df: (
                (df['i_phast17'] > 0.5) & 
                (df['j_phast17'] > 0.5) &
                (df['i_phast100'] > 0.5) & 
                (df['j_phast100'] > 0.5)
            )
        },
        'ml_structure': {
            'description': 'ML structure model predictions (> 0.247)',
            'filter': lambda df: df['pred_prob_editing_structure_only_alu100pct'] > 0.247
        },
        'ml_gtex': {
            'description': 'ML GTEx model predictions (> 0.251)',
            'filter': lambda df: df['pred_prob_editing_gtex_alu100pct'] > 0.251
        },
        'ml_structure_probing': {
            'description': 'ML structure probing model (> 0.0140)',
            'filter': lambda df: df['pred_3utr_All_power_weighted_advantage'] > 0.0140
        },
        'conserved_high_conf': {
            'description': 'Conserved + ML structure model',
            'filter': lambda df: (
                (df['i_phast17'] > 0.5) & 
                (df['j_phast17'] > 0.5) &
                (df['i_phast100'] > 0.5) & 
                (df['j_phast100'] > 0.5) &
                (df['pred_prob_editing_structure_only_alu100pct'] > 0.247)
            )
        },
        'alu': {
            'description': 'Alu-derived dsRNAs',
            'filter': lambda df: df['alu'] == 'Alu'
        },
        'nonalu': {
            'description': 'Non-Alu dsRNAs',
            'filter': lambda df: df['alu'] == 'Non-Alu'
        }
    }
    
    # Create subset files
    print("\nCreating subset files...")
    subset_info = []
    
    for subset_name, subset_config in subsets.items():
        print(f"\n{subset_name}: {subset_config['description']}")
        
        # Apply filter
        mask = subset_config['filter'](df)
        subset_df = df[mask].copy()
        
        print(f"  - {len(subset_df):,} dsRNAs selected")
        
        # Save as parquet (compact)
        parquet_file = output_dir / f"{subset_name}_dsrnas.parquet"
        subset_df.to_parquet(parquet_file, compression='snappy')
        print(f"  - Saved to {parquet_file}")
        
        # Skip BED file creation for large subsets (>100k rows)
        if len(subset_df) < 100000:
            bed_file = output_dir / f"{subset_name}_dsrnas.bed"
            create_bed_file(subset_df, bed_file)
            print(f"  - Created BED: {bed_file}")
        else:
            bed_file = None
            print(f"  - Skipped BED creation (too large: {len(subset_df):,} rows)")
        
        subset_info.append({
            'name': subset_name,
            'description': subset_config['description'],
            'count': len(subset_df),
            'parquet': str(parquet_file),
            'bed': str(bed_file) if bed_file else 'N/A (too large)'
        })
    
    # Create a sample subset for testing
    print("\nCreating sample data for testing...")
    sample_df = df.sample(n=min(1000, len(df)), random_state=42)
    sample_file = output_dir / "sample_dsrnas.parquet"
    sample_df.to_parquet(sample_file)
    print(f"  - Sample data (1000 dsRNAs): {sample_file}")
    
    # Save subset information
    info_file = output_dir / "subset_info.txt"
    with open(info_file, 'w') as f:
        f.write("dsRNA Subset Information\n")
        f.write("=" * 60 + "\n\n")
        for info in subset_info:
            f.write(f"{info['name']}:\n")
            f.write(f"  Description: {info['description']}\n")
            f.write(f"  Count: {info['count']:,}\n")
            f.write(f"  Files: {info['bed']}\n")
            f.write(f"         {info['parquet']}\n\n")
    
    print(f"\nSubset information saved to: {info_file}")
    
    # Skip creating BED for full dataset (too large)
    print("\nSkipping full dataset BED file (5M+ rows - use parquet instead)")
    
    return subset_info

def create_bed_file(df, output_file):
    """Create BED file from dsRNA dataframe"""
    with open(output_file, 'w') as f:
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Writing BED"):
            # Get strand, use '.' if not available
            strand = row['er_strand'] if 'er_strand' in row and pd.notna(row['er_strand']) else '.'
            
            # Write i-arm
            f.write(f"{row['er_chr']}\t{row['er_i_start']}\t{row['er_i_end']}\t")
            f.write(f"dsRNA_i_{idx}\t1000\t{strand}\n")
            
            # Write j-arm
            f.write(f"{row['er_chr']}\t{row['er_j_start']}\t{row['er_j_end']}\t")
            f.write(f"dsRNA_j_{idx}\t1000\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Prepare dsRNA subset files for overlap analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Prepare subsets from existing parquet file
  python prepare_data.py --input /path/to/dsrna_predictions.parquet
  
  # Download instructions
  python prepare_data.py --download
  
  # Use default path
  python prepare_data.py --input ../20250619.df_with_normalized_predictions.parquet
        """
    )
    
    parser.add_argument('--input', '-i', 
                       help='Path to input dsRNA predictions parquet file')
    parser.add_argument('--output-dir', '-o', default='data',
                       help='Output directory for subset files (default: data/)')
    parser.add_argument('--download', action='store_true',
                       help='Show download instructions for dsRNA data')
    
    args = parser.parse_args()
    
    if args.download:
        download_dsrna_data()
        return
    
    if not args.input:
        print("Error: Please provide input file with --input or use --download")
        print("\nTry looking for the file in these locations:")
        print("  - /Users/ryanandrews/Bioinformatics/20250619.df_with_normalized_predictions.parquet")
        print("  - ../dsrnascan_output/predictions.parquet")
        print("\nRun 'python prepare_data.py --help' for more information")
        sys.exit(1)
    
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Create subsets
    subset_info = create_subsets(args.input, args.output_dir)
    
    print("\n" + "=" * 60)
    print("DATA PREPARATION COMPLETE!")
    print("=" * 60)
    print(f"\nCreated {len(subset_info)} subset files in {args.output_dir}/")
    print("\nYou can now run overlap analysis with:")
    print("  python dsrna_overlap_analyzer.py your_features.bed")
    print("\nThe analyzer will automatically use the prepared subset files.")

if __name__ == "__main__":
    main()