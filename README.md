# dsRNAscan

[![CI](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci.yml/badge.svg)](https://github.com/Bass-Lab/dsRNAscan/actions/workflows/ci.yml)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool designed to identify and characterize **double-stranded RNAs (dsRNAs)** encoded by genomic sequences. This package also bundles a modified version of **einverted** from the EMBOSS suite for detecting inverted repeats that may form dsRNA structures.

---

## Installation & Running

### 1. Conda-Based Build (for Modified einverted + EMBOSS)

We provide a Conda-based recipe that compiles EMBOSS (including our modified einverted.c) and installs dsRNAscan:

1. Install :
   ```bash
   conda install -c conda-forge conda-build
   ```

2. Clone & build locally:
   ```bash
   git clone https://github.com/Bass-Lab/dsRNAscan.git
   cd dsRNAscan
   conda build conda-recipe
   ```

3. Install the newly built package in your environment:
   ```bash
   conda install --use-local dsrnascan
   ```

4. Verify:
   ```bash
   dsrnascan --help
   einverted --help
   ```

---

### 2. Using `python` Directly

If you prefer to run `dsRNAscan.py` with Python:

```bash
# Clone the repo
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan

# (Optional) Install requirements with pip
pip install numpy pandas biopython viennarna

# Run dsRNAscan via Python:
python dsRNAscan.py --help
```

### 3. Making dsRNAscan.py Executable

If you'd rather run dsRNAscan.py as a standalone script:

```bash
# Clone the repo
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan

# Make the script executable
chmod +x dsRNAscan.py

# (Optional) Install dependencies
pip install numpy pandas biopython viennarna

# Execute dsRNAscan directly:
./dsRNAscan.py --help
```

## Requirements

- Python 3.x
- Dependencies:
  - numpy
  - pandas
  - biopython
  - ViennaRNA (for advanced RNA folding support)
  - argparse, glob, multiprocessing (these are typically included with Python)

If you build dsRNAscan with our Conda recipe, these dependencies are handled automatically. If installing manually with pip, you can run:

```bash
pip install numpy pandas biopython viennarna
```

---

## Usage

### Basic Command-Line Usage:

```bash
./dsRNAscan.py <input_fasta_file> --chr <chromosome_name> [other flags...]
```

**Note:** Input files can be compressed (.gz) or uncompressed FASTA format.

### Key Parameters:

- `--input`: Path to the input FASTA file.
- `--chr`: Chromosome name. If "header", uses the sequence header as the chromosome name.
- `--start`: Starting coordinate for scanning (default: 0).
- `--end`: Ending coordinate for scanning (default: 0 → full sequence).
- `--step-size`: Step size for each window (default: 150).
- `--window-size`: Window size for scanning (default: 10,000).
- `--min-length`: Minimum length for dsRNA prediction (default: 30).
- `--max_span`: Maximum span of inverted repeats (default: 10,000).
- `--gaps`: Gap penalty (default: 12).
- `--match`: Match score (default: 3).
- `--mismatch`: Mismatch penalty (default: -4).
- `--paired_cutoff`: Minimum % of base-pairing to keep (default: 70).
- `--algorithm`: "einverted" or "iupacpal" (default: "einverted").
- `--reverse`: Flag to scan the reverse strand.
- `--cpus`: Number of CPUs for parallel processing (default: 4).

### Example:

```bash
./dsRNAscan.py --input genome.fasta --chr header --window-size 10000 --step-size 150
```

---

## Output

1. **Results File** (`<basename>_merged_results.txt`)
   
   Tab-delimited file with the following columns:
   - **Chromosome**: Chromosome/sequence name
   - **Strand**: + (forward) or - (reverse)
   - **Score**: Einverted alignment score
   - **RawMatch**: Number of matching base pairs
   - **PercMatch**: Percentage of matching base pairs
   - **Gaps**: Number of gaps in alignment
   - **i_start, i_end**: Start and end positions of the first arm
   - **j_start, j_end**: Start and end positions of the second arm
   - **eff_i_start, eff_i_end**: Effective coordinates of first arm after RNAduplex
   - **eff_j_start, eff_j_end**: Effective coordinates of second arm after RNAduplex
   - **i_seq, j_seq**: Sequences of both arms
   - **structure**: RNA secondary structure in dot-bracket notation
   - **dG(kcal/mol)**: Free energy of the structure
   - **percent_paired**: Percentage of bases that are paired
   - **longest_helix**: Length of the longest continuous helix
   - **i_eff_seq, j_eff_seq**: Effective sequences after RNAduplex adjustment
   - **FORNA_Link**: URL for FORNA visualization
   - **orig_arm_length**: Original arm length from einverted (v0.3.0+)
   - **eff_arm_length**: Effective arm length after RNAduplex refinement (v0.3.0+)

2. **BP File** (`<basename>.dsRNApredictions.bp`)
   
   Tab-delimited file for IGV visualization with:
   - Color definitions for different pairing percentages
   - Base pair coordinates in IGV-compatible format

---

## Genome-Wide Scanning

### Quick Start for Large Genomes

For scanning compressed genome files:
```bash
# No need to decompress - dsRNAscan now supports .gz files directly!
# Fast initial scan
dsrnascan genome.fa.gz --chr chr1 -w 20000 -s 10000 --min 100 -c 8

# Detailed scan of specific region
dsrnascan genome.fa --chr chr1 --start 1000000 --end 2000000 -w 10000 -s 2000
```

### Parallel Processing
```bash
# Process multiple chromosomes in parallel
parallel -j 4 'dsrnascan chr{}.fa --chr chr{} -w 10000 -s 5000' ::: {1..22} X Y
```

### Performance Tips
- Use larger windows (`-w 20000`) and steps (`-s 10000`) for initial screening
- Increase `--score` threshold to reduce false positives in repetitive regions
- Process chromosomes separately to manage memory usage
- See `GENOME_SCANNING_GUIDE.md` for detailed strategies

---

## Logging

dsRNAscan.py can produce console output logs (and you can add more debug statements if needed).

---

## Manual EMBOSS/einverted Build

If you don't use Conda, you can build EMBOSS with our patched einverted.c via the provided build.sh:

1. Clone this repository
2. Run build.sh in conda-recipe/ (or wherever you placed it), which downloads EMBOSS, patches config.sub, and compiles the local einverted.
3. Use the newly compiled einverted in your dsRNAscan runs.

---

## Contributing

We welcome pull requests! Please ensure your code is well-documented and tested.

---

## License

This project is licensed under the GPL License.

---

## Contact

For questions or more information, contact:
Ryan Andrews ryan.andrews@biochem.utah.edu
