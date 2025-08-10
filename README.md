# dsRNAscan

[![CI Tests](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci-simple.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with special support for **G-U wobble base pairing** critical for RNA analysis.

## 🚀 Quick Start

### Install from PyPI (Coming Soon)
```bash
pip install dsrnascan
```

### Install from GitHub  
```bash
# Direct from GitHub with automatic einverted compilation
pip install --no-binary :all: git+https://github.com/Bass-Lab/dsRNAscan.git

# Or if you prefer using pre-built wheels (faster but einverted needs separate setup)
pip install git+https://github.com/Bass-Lab/dsRNAscan.git
```

### Install from Local Files
```bash
# Option 1: Clone and install with einverted compilation
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
pip install --no-binary :all: .

# Option 2: Quick install (uses pre-compiled binaries if available)
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
pip install .

# Option 3: Development mode (editable install)
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
pip install -e .

# Option 4: Manual einverted compilation then install
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
./compile_patched_einverted.sh  # Compile einverted with G-U patch
pip install .
```

### Basic Usage
```bash
# Scan a genome/sequence for dsRNA structures
dsrnascan input.fasta -w 10000 -s 150 --score 50

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

## 🧬 Key Features

- **G-U Wobble Base Pairing**: Modified einverted algorithm specifically for RNA
- **Parallel Processing**: Multi-CPU support for genome-scale analysis
- **Flexible Windowing**: Customizable window and step sizes
- **RNA Structure Prediction**: Integration with ViennaRNA for structure validation
- **Multiple Output Formats**: Tab-delimited results and IGV visualization files

## 📖 Detailed Usage

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

## 🔧 Installation Troubleshooting

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

## 🧪 Testing

Run test with sample data:
```bash
# Create test file
echo ">test_sequence" > test.fasta
echo "GGGGGGGGGGAAAAAAAAAAAAAACCCCCCCCCC" >> test.fasta

# Run dsRNAscan
dsrnascan test.fasta -w 100 -s 50 --score 15
```

## 📚 Algorithm Details

dsRNAscan uses a multi-step approach:

1. **Window Extraction**: Divides genome into overlapping windows
2. **Inverted Repeat Detection**: Uses modified einverted with G-U wobble support
3. **Structure Prediction**: Validates structures with RNAduplex (ViennaRNA)
4. **Filtering**: Applies score and pairing percentage cutoffs
5. **Parallel Processing**: Distributes windows across multiple CPUs

The key innovation is the **G-U wobble patch** for einverted, allowing detection of RNA-specific base pairs crucial for identifying functional dsRNA structures.

## 📄 Citation

If you use dsRNAscan in your research, please cite:
```
Bass Lab. dsRNAscan: A tool for genome-wide prediction of double-stranded RNA structures.
https://github.com/Bass-Lab/dsRNAscan
```

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## 📝 License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

- **Issues**: [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)
- **Documentation**: [GitHub Wiki](https://github.com/Bass-Lab/dsRNAscan/wiki)

## Acknowledgments

- EMBOSS team for the einverted tool
- ViennaRNA team for RNA folding algorithms
- All contributors to the project

---
**Note**: This tool is for research purposes. Ensure you understand the parameters for your specific use case.