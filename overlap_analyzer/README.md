# dsRNA Overlap Analyzer

Statistical enrichment analysis tool for identifying genomic features that significantly overlap with dsRNA regions predicted by dsRNAscan.

## Overview

This tool analyzes whether your genomic features (e.g., RNA-binding protein sites, regulatory elements, disease variants) are enriched or depleted in dsRNA regions compared to random expectation. It uses permutation testing to calculate robust statistical significance.

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run with example data
python dsrna_overlap_analyzer.py examples/microexons_30bp.bed

# Analyze your own features
python dsrna_overlap_analyzer.py your_features.bed --subset conserved_high_conf
```

## Installation

```bash
cd overlap_analyzer
pip install -r requirements.txt
```

Requirements:
- Python 3.7+
- pandas, numpy, pybedtools, tqdm
- pyarrow (for parquet file support)

## Data Setup

### Option 1: Quick Start with Local Data
If you have the dsRNA predictions parquet file:
```bash
# Prepare subset files for faster analysis
python prepare_data.py --input /path/to/dsrna_predictions.parquet

# This creates optimized subset files in data/ directory
```

### Option 2: Download Data
```bash
# Get download instructions
./download_data.sh
# or
python prepare_data.py --download
```

### Option 3: Use dsRNAscan Output Directly
```bash
# The analyzer can use dsRNAscan output directly
python dsrna_overlap_analyzer.py features.bed \
    --dsrna-file ../dsrnascan_output/predictions.parquet
```

After data preparation, the `data/` directory will contain:
- Pre-filtered subset files (BED and parquet formats)
- `subset_info.txt` with subset statistics
- Sample data for testing

## Usage

### Basic Usage

```bash
python dsrna_overlap_analyzer.py <input_bed> [options]
```

### Input Format

The tool accepts BED format files with at least 3 columns (chr, start, end):
```
chr1    1000    2000    feature1    100    +
chr2    5000    6000    feature2    200    -
```

### dsRNA Data Options

The analyzer looks for dsRNA data in this order:
1. User-specified `--dsrna-file` argument
2. Pre-processed BED files in `data/` directory
3. Default parquet file with all predictions
4. Environment variable `$DSRNA_DATA_PATH`

### Available dsRNA Subsets

- `conserved_high_conf` - High-confidence conserved dsRNAs (default)
  - PhastCons conservation > 0.5
  - ML structure model confidence > 0.247
  
- `conserved` - All conserved dsRNAs
  - PhastCons conservation > 0.5
  
- `ml_structure` - ML structure model predictions
  - Structure-only model confidence > 0.247
  
- `ml_gtex` - GTEx model predictions  
  - GTEx model confidence > 0.251
  
- `all` - All predicted dsRNAs without filtering

### Command Line Options

```
--subset SUBSET          dsRNA subset to use (default: conserved_high_conf)
--dsrna-file FILE       Custom dsRNA predictions file (parquet or BED)
--permutations N        Number of permutations for significance testing (default: 100)
--output-format FORMAT  Output format: csv, json, or both (default: csv)
--output-prefix PREFIX  Prefix for output files (default: overlap_results)
--min-overlap N         Minimum overlap in base pairs (default: 1)
--reciprocal FRAC      Minimum reciprocal overlap fraction (0-1)
--strand-specific      Consider strand in overlap analysis
--verbose              Enable detailed logging
--debug                Enable debug mode with intermediate files
```

## Output

### Console Output
- Summary statistics with feature and dsRNA counts
- Overlap results with Z-scores and p-values
- Interpretation of enrichment/depletion
- Recommendations for additional analyses

### File Outputs

**CSV Format** (`overlap_results.csv`):
```csv
metric,value
total_features,1234
total_dsrnas,50000
observed_overlaps,567
expected_overlaps,234.5
z_score,15.6
p_value,1.2e-10
fold_enrichment,2.42
```

**JSON Format** (`overlap_results.json`):
```json
{
  "metadata": {...},
  "overlap_stats": {...},
  "significance": {...},
  "permutation_stats": {...}
}
```

## Examples

### Analyzing RNA-Binding Protein Sites
```bash
# RBP binding sites from eCLIP
python dsrna_overlap_analyzer.py rbp_peaks.bed --subset conserved --permutations 1000
```

### Testing Disease Variants
```bash
# GWAS variants in dsRNA regions
python dsrna_overlap_analyzer.py gwas_snps.bed --subset ml_gtex --min-overlap 1
```

### Regulatory Elements
```bash
# Enhancers overlapping with dsRNAs
python dsrna_overlap_analyzer.py enhancers.bed --reciprocal 0.5 --strand-specific
```

### Using Custom dsRNA Predictions
```bash
# From dsRNAscan output
python dsrna_overlap_analyzer.py features.bed \
    --dsrna-file ../dsrnascan_output/predictions.parquet \
    --permutations 500
```

## Statistical Methods

### Permutation Testing
1. Calculate observed overlaps between features and dsRNAs
2. Generate N random permutations by shuffling genomic positions
3. Calculate overlap for each permutation
4. Compute Z-score: (observed - mean_expected) / std_expected
5. Calculate empirical p-value from permutation distribution

### Significance Interpretation
- **Z-score > 3**: Strong enrichment (p < 0.001)
- **Z-score > 2**: Significant enrichment (p < 0.05)
- **Z-score < -2**: Significant depletion (p < 0.05)
- **-2 ≤ Z-score ≤ 2**: No significant enrichment/depletion

## Integration with dsRNAscan

This tool is designed to work seamlessly with dsRNAscan predictions:

```bash
# Run dsRNAscan
dsrnascan genome.fa --output-dir dsrna_predictions/

# Convert to parquet (if needed)
python convert_to_parquet.py dsrna_predictions/*_merged_results.txt

# Analyze overlaps
cd overlap_analyzer
python dsrna_overlap_analyzer.py my_features.bed \
    --dsrna-file ../dsrna_predictions/predictions.parquet
```

## Performance Tips

1. **Large Files**: Use pre-processed BED files (faster than parquet)
2. **Memory**: Process in chunks with `--chunk-size` for very large datasets
3. **Speed**: Reduce permutations for initial exploration (--permutations 100)
4. **Accuracy**: Increase permutations for publication (--permutations 1000+)

## Troubleshooting

### "No dsRNA data file found"
- Check file paths and permissions
- Set environment variable: `export DSRNA_DATA_PATH=/path/to/dsrna.parquet`
- Use `--dsrna-file` to specify location explicitly

### "Memory Error"
- Use a smaller subset: `--subset conserved_high_conf`
- Reduce permutations temporarily
- Process features in batches

### "No overlaps found"
- Check chromosome naming (chr1 vs 1)
- Verify coordinate systems match (0-based vs 1-based)
- Try `--min-overlap 1` for single-bp overlaps

## Citation

If you use this tool in your research, please cite:

Comprehensive mapping of human dsRNAome reveals conservation, neuronal enrichment, and intermolecular interactions
https://doi.org/10.1101/2025.01.24.634786

## Support

- Issues: [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)
- Main tool: [dsRNAscan](https://github.com/Bass-Lab/dsRNAscan)