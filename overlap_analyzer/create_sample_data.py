#!/usr/bin/env python3
"""
Create minimal sample data for testing the overlap analyzer
This creates a small parquet file with ~100 dsRNAs for quick testing
"""

import pandas as pd
import numpy as np
from pathlib import Path

def create_sample_dsrna_data():
    """Create a small sample dsRNA dataset"""
    
    np.random.seed(42)
    
    # Create sample data with required columns
    n_samples = 100
    
    # Generate chromosomes
    chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    
    data = {
        'er_chr': np.random.choice(chromosomes, n_samples),
        'er_i_start': np.random.randint(1000, 100000000, n_samples),
        'er_i_end': np.zeros(n_samples, dtype=int),
        'er_j_start': np.zeros(n_samples, dtype=int),
        'er_j_end': np.zeros(n_samples, dtype=int),
        'er_strand': np.random.choice(['+', '-'], n_samples),
        
        # Conservation scores (using both naming conventions)
        'i_phast17': np.random.uniform(0, 1, n_samples),
        'j_phast17': np.random.uniform(0, 1, n_samples),
        'i_phast100': np.random.uniform(0, 1, n_samples),
        'j_phast100': np.random.uniform(0, 1, n_samples),
        'PhastCons17_Ave_i': np.random.uniform(0, 1, n_samples),
        'PhastCons17_Ave_j': np.random.uniform(0, 1, n_samples),
        'PhastCons100_Ave_i': np.random.uniform(0, 1, n_samples),
        'PhastCons100_Ave_j': np.random.uniform(0, 1, n_samples),
        
        # ML predictions
        'pred_prob_editing_structure_only_alu100pct': np.random.uniform(0, 1, n_samples),
        'pred_prob_editing_gtex_alu100pct': np.random.uniform(0, 1, n_samples),
        'pred_3utr_All_power_weighted_advantage': np.random.uniform(-0.05, 0.15, n_samples),
        
        # Alu status
        'alu': np.random.choice(['Alu', 'Non-Alu'], n_samples),
        
        # Additional metadata
        'er_length_i': np.random.randint(30, 500, n_samples),
        'er_length_j': np.random.randint(30, 500, n_samples),
    }
    
    # Calculate end positions based on lengths
    data['er_i_end'] = data['er_i_start'] + data['er_length_i']
    data['er_j_start'] = data['er_i_end'] + np.random.randint(50, 5000, n_samples)
    data['er_j_end'] = data['er_j_start'] + data['er_length_j']
    
    df = pd.DataFrame(data)
    
    # Make column values consistent (both naming conventions should have same values)
    df['PhastCons17_Ave_i'] = df['i_phast17']
    df['PhastCons17_Ave_j'] = df['j_phast17']
    df['PhastCons100_Ave_i'] = df['i_phast100']
    df['PhastCons100_Ave_j'] = df['j_phast100']
    
    # Ensure some dsRNAs meet various filter criteria
    # Make 20% highly conserved
    conserved_idx = np.random.choice(n_samples, size=20, replace=False)
    for idx in conserved_idx:
        value_i17 = np.random.uniform(0.6, 1.0)
        value_j17 = np.random.uniform(0.6, 1.0)
        value_i100 = np.random.uniform(0.6, 1.0)
        value_j100 = np.random.uniform(0.6, 1.0)
        
        df.loc[idx, 'i_phast17'] = value_i17
        df.loc[idx, 'j_phast17'] = value_j17
        df.loc[idx, 'i_phast100'] = value_i100
        df.loc[idx, 'j_phast100'] = value_j100
        df.loc[idx, 'PhastCons17_Ave_i'] = value_i17
        df.loc[idx, 'PhastCons17_Ave_j'] = value_j17
        df.loc[idx, 'PhastCons100_Ave_i'] = value_i100
        df.loc[idx, 'PhastCons100_Ave_j'] = value_j100
    
    # Make 30% high ML confidence
    ml_idx = np.random.choice(n_samples, size=30, replace=False)
    for idx in ml_idx:
        df.loc[idx, 'pred_prob_editing_structure_only_alu100pct'] = np.random.uniform(0.3, 1.0)
    
    return df

def main():
    # Create sample data
    print("Creating sample dsRNA data...")
    df = create_sample_dsrna_data()
    
    # Create data directory
    data_dir = Path('data')
    data_dir.mkdir(exist_ok=True)
    
    # Save as parquet
    output_file = data_dir / 'sample_dsrnas_test.parquet'
    df.to_parquet(output_file, compression='snappy')
    
    print(f"Created sample data with {len(df)} dsRNAs")
    print(f"Saved to: {output_file}")
    
    # Show subset counts
    conserved = (
        (df['i_phast17'] > 0.5) & 
        (df['j_phast17'] > 0.5) &
        (df['i_phast100'] > 0.5) & 
        (df['j_phast100'] > 0.5)
    ).sum()
    
    ml_high = (df['pred_prob_editing_structure_only_alu100pct'] > 0.247).sum()
    alu = (df['alu'] == 'Alu').sum()
    
    print(f"\nSubset counts in sample:")
    print(f"  - Conserved: {conserved}")
    print(f"  - ML high confidence: {ml_high}")
    print(f"  - Alu-derived: {alu}")
    print(f"  - Non-Alu: {len(df) - alu}")
    
    print("\nYou can now test the overlap analyzer with:")
    print("  python dsrna_overlap_analyzer.py examples/microexons_30bp.bed \\")
    print("    --dsrna-file data/sample_dsrnas_test.parquet")

if __name__ == "__main__":
    main()