# dsRNAscan

[![CI Tests](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with special support for **G-U wobble base pairing**.

### Install from PyPI 
```bash
pip install dsrnascan
```

### Basic Usage
```bash
# Scan a genome/sequence for dsRNA structures
dsrnascan input.fasta # This uses defaults of -w 10000 -s 150 --score 50

# Process specific chromosome
dsrnascan genome.fasta --chr chr21 -c 8

# Use custom parameters for sensitive detection
dsrnascan sequence.fasta -w 5000 --min 20 --score 30
```

## 📋 Requirements

- **Python 3.8+**
- **Dependencies** (automatically installed):
  - numpy ≥1.19
  - pandas ≥1.1
  - biopython ≥1.78
  - ViennaRNA ≥2.4

### Important: einverted Binary

dsRNAscan requires the `einverted` tool from EMBOSS with our **G-U wobble patch** for accurate RNA structure detection. 

**Option 1: Automatic** (macOS with included binary)
- The package includes a pre-compiled einverted for macOS ARM64
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

Key parameters:
- `-w/--window`: Window size for scanning (default: 10000)
- `-s/--step`: Step size between windows (default: 150)
- `--score`: Minimum score threshold for inverted repeats (default: 50)
- `--min/--max`: Min/max length of inverted repeats (default: 30/10000)
- `--paired_cutoff`: Minimum percentage of paired bases (default: 70%)
- `-c/--cpus`: Number of CPUs to use (default: 4)
- `--chr`: Specific chromosome to process
- `--reverse`: Scan reverse strand

### Output Files

dsRNAscan generates several output files in a timestamped directory:

1. **`*_merged_results.txt`**: Tab-delimited file with all predicted dsRNAs
   - Columns include: coordinates, scores, sequences, structures, folding energy
   
2. **`*.dsRNApredictions.bp`**: IGV-compatible visualization file
   - Load in IGV to visualize dsRNA locations on genome

### Example Workflow

```bash
# 1. Basic genome scan
dsrnascan human_genome.fa -c 16 --output-dir results/

# 2. Scan specific region with sensitive parameters
dsrnascan chr21.fa -w 5000 -s 100 --score 30 --min 20

# 3. Process RNA-seq assembled transcripts
dsrnascan transcripts.fa -w 1000 --paired_cutoff 60

# 4. Scan both strands
dsrnascan sequence.fa --reverse
```

## Installation Troubleshooting

### "einverted binary not found"
The package needs einverted from EMBOSS. Solutions:
1. Install EMBOSS: `conda install -c bioconda emboss`
2. Or compile during install: `DSRNASCAN_COMPILE_FULL=true pip install .`
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
```bash
module load python/3.8  # or your Python module
module load emboss      # if available
pip install --user git+https://github.com/Bass-Lab/dsRNAscan.git
```


## Using dsRNAscan as a Python Module

While primarily designed as a command-line tool, dsRNAscan can be imported and used in Python scripts:

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

## Algorithm Details

dsRNAscan uses a multi-step approach:

1. **Window Extraction**: Divides genome into overlapping windows
2. **Inverted Repeat Detection**: Uses modified einverted with G-U wobble support
3. **Structure Prediction**: Validates structures with RNAduplex (ViennaRNA)
4. **Filtering**: Applies score and pairing percentage cutoffs
5. **Parallel Processing**: Distributes windows across multiple CPUs

The key innovation is the **G-U wobble patch** for einverted, allowing detection of RNA-specific base pairs crucial for identifying functional dsRNA structures.

## Citation

If you use dsRNAscan in your research, please cite:
Comprehensive mapping of human dsRNAome reveals conservation, neuronal enrichment, and intermolecular interactions

https://doi.org/10.1101/2025.01.24.634786



## Additional Tools

- **overlap_analyzer/** - Statistical enrichment analysis for genomic features overlapping with dsRNA predictions. See [overlap_analyzer/README.md](overlap_analyzer/README.md) for details.

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