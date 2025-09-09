# dsRNAscan

[![CI Tests](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml)
[![Python](https://img.shields.io/badge/Python-3.8%2B%20(Linux)%20|%203.9%2B%20(macOS)-blue.svg)](https://www.python.org/downloads/)
[![Platforms](https://img.shields.io/badge/Platforms-Linux%20|%20macOS-green.svg)](https://github.com/Bass-Lab/dsRNAscan)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with special support for **G-U wobble base pairing**.

You can browse human genome results at dsrna.chpc.utah.edu 

### Install from PyPI 
```bash
pip install dsrnascan
# Note: You may need to compile einverted with G-U patch for non-macOS systems
```

### Basic Usage
```bash
# Scan a genome/sequence for dsRNA structures
dsrnascan input.fasta # This uses defaults of -w 10000 -s 150 --score 50 -c 4 (cpus)

# Process specific chromosome, using 8 cpus (-c)
dsrnascan genome.fasta --only_seq chr21 -c 8

# Use custom parameters for detecting smaller structures (like a minimum of 15bp)
dsrnascan sequence.fasta -w 5000 --min_bp 15
```

## 📋 Requirements

### Platform Compatibility
- **Linux**: ✅ Python 3.8+ 
- **macOS**: ✅ Python 3.9+ (3.8 not supported)
- **Windows**: ❌ Not supported (use WSL or Docker)

### Dependencies (automatically installed):
  - numpy ≥1.19
  - pandas ≥1.1
  - biopython ≥1.78
  - ViennaRNA ≥2.4

### Important: einverted Binary

dsRNAscan requires the `einverted` tool from EMBOSS with our **G-U wobble patch** for accurate RNA structure detection. 

**Option 1: Automatic** (included binary)
- The package includes a pre-compiled einverted for some systems (MacOS ARM)
- It will be used automatically on compatible systems

**Option 2: System Installation** (Linux/Other)
```bash
# Ubuntu/Debian
sudo apt-get install emboss

# macOS with Homebrew
brew install emboss

# Conda (recommended for bioinformatics workflows)
conda install -c bioconda emboss
```

**Note:** System-installed EMBOSS won't have the G-U patch. For full RNA functionality with G-U wobble pairs, compile from source:

```bash
# Compile with G-U patch (optional but recommended)
cd dsRNAscan
DSRNASCAN_COMPILE_FULL=true pip install .
```

## Detailed Usage

### Command-Line Options

```bash
dsrnascan --help
```

### Complete Parameter Reference

**Core Parameters:**
- `-w`: Window size for scanning (default: 10000)
- `-s/--step`: Step size between windows (default: 150)
- `-t`: Folding temperature in Celsius (default: 37)

**Structure Requirements:**
- `--min_bp`: Minimum number of base pairs required (default: 25) - **Recommended**
- `--score`: Minimum score threshold for inverted repeat (default: 75) - *Deprecated, use --min_bp*
- `--paired_cutoff`: Minimum percentage of paired bases (default: 70)
- `--min`: Minimum length of inverted repeat (default: 30)
- `--max`: Maximum length of inverted repeat (default: 10000)
- `--max_span`: Maximum span of inverted repeat (default: window size)

**Region Selection:**
- `--only_seq`: Process only this specific sequence/chromosome (based on fasta header)
- `--start`: Starting coordinate for scan (default: 0, 1-based)
- `--end`: Ending coordinate for scan (default: 0 = end of sequence)

**Strand Options:**
- `--forward-only`: Process forward strand only
- `--reverse-only`: Process reverse strand only
- Default: both strands are processed

**Scoring Parameters:**
- `--match`: Match score (default: 3)
- `-x/--mismatch`: Mismatch score (default: -4)
- `--gaps`: Gap penalty (default: 12)

**Algorithm Options:**
- `--algorithm`: Inverted repeat algorithm (einverted or iupacpal)
- `--eliminate-nested`: Remove nested dsRNAs (default: True)
- `--chunk-size`: Windows per chunk for DataFrame processing (default: 10000)

**Output Options:**
- `--output-dir`: Output directory (default: dsrnascan_YYYYMMDD_HHMMSS)
- `--output_label`: Label for output files (default: sequence header)
- `--clean`: Clean up temporary files after processing

**Performance:**
- `-c/--cpus`: Number of CPUs to use (default: 4)

**Other Options:**
- `--version`: Show program version
- `-h/--help`: Show help message
- `--batch`: DEPRECATED - only with --legacy flag
- `--legacy`: DEPRECATED - use legacy non-DataFrame approach (slower)

### Output Files

dsRNAscan generates several output files in a timestamped directory:

1. **`*_merged_results.txt`**: Tab-delimited file with all predicted dsRNAs
   
   **Column Groups:**
   - **Genomic Coordinates** (Columns 1-6): `Chromosome`, `Strand`, `i_start`, `i_end`, `j_start`, `j_end`
   - **einverted Results** (Columns 7-10): `Score`, `RawMatch`, `PercMatch`, `Gaps`
     - These come from the inverted repeat detection by einverted
   - **RNAduplex Results** (Columns 11-19): `dG(kcal/mol)`, `percent_paired`, `longest_helix`, `eff_i_start`, `eff_i_end`, `eff_j_start`, `eff_j_end`, `i_seq`, `j_seq`, `structure`
     - These come from RNA secondary structure prediction by RNAduplex
     - The effective coordinates show the trimmed regions that form the optimal dsRNA structure
     - Sequences are reported in 5' to 3' RNA orientation (reverse complement for minus strand)
   
2. **`*.dsRNApredictions.bp`**: IGV-compatible visualization file
   - Load in IGV to visualize dsRNA locations on genome

### Example Workflows

```bash
# 1. Basic genome-wide scan with 16 CPUs
dsrnascan genome.fa -c 16 --output-dir results/

# 2. Scan specific genomic region (e.g., 200kb region on chr21)
dsrnascan hg38_chromosomes.fa.gz \
    --only_seq chr21 \
    --start 33455482 \
    --end 33655482 \
    -w 10000 -s 5000 \
    --score 75

# 3. Scan multiple chromosomes
dsrnascan genome.fa --only_seq chr1,chr2,chr3 -c 8

# 4. Sensitive scan for shorter dsRNAs
dsrnascan sequence.fa \
    -w 5000 -s 100 \
    --score 30 \
    --min 20 \
    --paired_cutoff 60

# 5. Process RNA-seq assembled transcripts
dsrnascan transcripts.fa \
    -w 1000 -s 50 \
    --paired_cutoff 60 \
    --min 25

# 6. Scan both strands (forward and reverse)
dsrnascan sequence.fa --both

# 7. Scan only reverse strand
dsrnascan sequence.fa --reverse

# 8. Quick test run on small region
dsrnascan test.fa -w 100 -s 50 --score 15 --min 10
```

#### Region-Specific Scanning

dsRNAscan supports scanning specific genomic regions, which is useful for:
- Focusing on regions of interest (e.g., gene loci, QTL regions)
- Testing parameters on small regions before genome-wide runs
- Reducing computational time for targeted analysis

```bash
# Scan a 1MB region on chromosome 21
dsrnascan hg38.fa.gz \
    --only_seq chr21 \
    --start 30000000 \
    --end 31000000 \
    -w 10000 -s 1000

# Scan around a specific gene (e.g., 50kb upstream and downstream)
# If gene is at chr1:1000000-1050000
dsrnascan genome.fa \
    --only_seq chr1 \
    --start 950000 \
    --end 1100000
```

## Installation Troubleshooting

### Note: einverted with G-U Wobble Pairing Support

**IMPORTANT**: dsRNAscan requires a patched version of einverted that recognizes G-U wobble base pairs as matches. Standard EMBOSS einverted treats G-U as mismatches, which misses many RNA structures.

#### Option 1: Use Pre-compiled Binary (macOS ARM64 only)
The PyPI package includes a pre-compiled einverted for macOS ARM64 (Apple Silicon).

#### Option 2: Compile Patched einverted (Recommended for other platforms)
```bash
# The package includes the patch and compilation script
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
./compile_patched_einverted.sh
```

This script will:
1. Download EMBOSS 6.6.0 source code
2. Apply the G-U wobble pairing patch (`einverted.patch`)
3. Compile einverted with RNA-aware scoring
4. Install it to `dsrnascan/tools/einverted`

#### Option 3: Manual Compilation
```bash
# Download and extract EMBOSS
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar -xzf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0/emboss

# Apply the G-U patch (included in dsRNAscan package)
patch -p0 < /path/to/dsrnascan/einverted.patch

# Compile just einverted
gcc -O2 -o einverted einverted.c \
    -I../ajax/core -I../ajax/ajaxdb -I../ajax/acd \
    -L../ajax/core/.libs -L../ajax/ajaxdb/.libs -L../ajax/acd/.libs \
    -lajax -lajaxdb -lacd -lm -lz

# Copy to dsRNAscan tools directory
cp einverted /path/to/dsrnascan/tools/
```

#### Option 4: Use Standard EMBOSS (Not Recommended)
```bash
conda install -c bioconda emboss
```
**Warning**: Standard einverted will miss RNA structures with G-U wobble pairs, significantly reducing sensitivity for dsRNA detection.

### "einverted binary not found" Error
If you get this error, einverted is not in your PATH. Solutions:
1. Compile the patched version as shown above
2. Set environment variable: `export EINVERTED_PATH=/path/to/einverted`
3. Or use the package without functional testing: `dsrnascan --help` works without einverted

### "ModuleNotFoundError: No module named 'ViennaRNA'"
Install ViennaRNA Python bindings:
```bash
# Via conda (recommended)
conda install -c bioconda viennarna

# Via pip
pip install ViennaRNA
```

### Installation on HPC/Cluster

**Important**: Cluster EMBOSS modules have standard einverted which lacks G-U wobble support. You need to compile the patched version:

```bash
# Load Python module
module load python/3.8  # or your cluster's Python module

# Install dsRNAscan
pip install --user dsrnascan

# Clone repo to get compilation script and patch
git clone https://github.com/Bass-Lab/dsRNAscan.git ~/dsRNAscan_source

# Compile patched einverted in your home directory
cd ~/dsRNAscan_source
./compile_patched_einverted.sh

# Copy the compiled einverted to a location in your PATH or set environment variable
mkdir -p ~/bin
cp dsrnascan/tools/einverted ~/bin/
export PATH=$HOME/bin:$PATH

# Or set EINVERTED_PATH environment variable
export EINVERTED_PATH=$HOME/dsRNAscan_source/dsrnascan/tools/einverted

# Add to your ~/.bashrc or job submission script
echo 'export EINVERTED_PATH=$HOME/dsRNAscan_source/dsrnascan/tools/einverted' >> ~/.bashrc
```

For job submission scripts:
```bash
#!/bin/bash
#SBATCH --job-name=dsrnascan
#SBATCH --cpus-per-task=16

module load python/3.8
export EINVERTED_PATH=$HOME/bin/einverted  # Use your compiled version

dsrnascan genome.fa -c 16 --output-dir results/
```


## Using dsRNAscan as a Python Module

While primarily designed as a standalone tool, dsRNAscan can be imported and used in Python scripts:

```python
# Method 1: Simple usage
from dsrnascan import main
import sys

# Simulate command line arguments
sys.argv = ['dsrnascan', 'input.fasta', '-w', '1000', '--score', '30']
main()

# Method 2: Using subprocess for better control
import subprocess
result = subprocess.run(['dsrnascan', 'input.fasta', '--score', '30'], 
                       capture_output=True, text=True)

# Method 3: Parse results programmatically
import pandas as pd
import glob

# Run dsRNAscan
subprocess.run(['dsrnascan', 'input.fasta'])

# Find and read results
output_dir = sorted(glob.glob('dsrnascan_*'))[-1]
results = pd.read_csv(f"{output_dir}/*_merged_results.txt", sep='\t')
```

For more examples, see `using_dsrnascan_as_module.py` in the repository.

## Citation

If you use dsRNAscan in your research, please cite:
Comprehensive mapping of human dsRNAome reveals conservation, neuronal enrichment, and intermolecular interactions

https://doi.org/10.1101/2025.01.24.634786

## Additional Tools 

### dsrna-browse - Interactive Results Viewer with RNA Editing Support

**dsrna-browse** is an interactive web-based viewer for dsRNAscan results, featuring:
- Fornac RNA secondary structure visualization
- Interactive dropdown selection of dsRNA predictions
- Detailed structure metrics (free energy, base pairs, helix length)
- RNA editing site annotation from BED or GFF3 files

#### Basic Usage
```bash
# Browse results in current directory
dsrna-browse

# Browse results in specific output directory
dsrna-browse dsrnascan_20250120_143022/

# Use custom port
dsrna-browse --port 8888

# Don't auto-open browser
dsrna-browse --no-browser
```

#### With RNA Editing Sites
```bash
# Annotate with editing sites from BED file
dsrna-browse dsrnascan_output/ --editing-file editing_sites.bed

# Annotate with editing sites from GFF3 file
dsrna-browse dsrnascan_output/ --editing-file editing_sites.gff3
```

The viewer supports both BED and GFF3 formats for editing sites:
- **BED format**: chr, start, end, name, score (0-1000), strand
- **GFF3 format**: Automatically detects editing-related features

Editing sites are visualized with green gradient coloring:
- Dark green: High-frequency sites (≥80%)
- Medium green: Medium-frequency sites (30-80%)
- Light green: Low-frequency sites (<30%)

The viewer will:
1. Process all `*_merged_results.txt` files in the directory
2. Map editing sites to dsRNA structure positions (strand-aware)
3. Start a local web server (default port 8080)
4. Open your browser to display an interactive interface
5. Show RNA structures with editing annotations

Press Ctrl+C to stop the server when done.

### overlap_analyzer 


Statistical enrichment analysis for genomic features overlapping with dsRNA predictions. See [overlap_analyzer/README.md](overlap_analyzer/README.md) for details.

Note: overlap_analyzer is not included in the PyPI package to reduce size. Clone the repository to access it.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)
- **Documentation**: [GitHub Wiki](https://github.com/Bass-Lab/dsRNAscan/wiki)

## Acknowledgments

- EMBOSS team for the einverted tool
- ViennaRNA team for RNA folding algorithms

---
**Note**: This tool is for research purposes. Ensure you understand the parameters for your specific use case.