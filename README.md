# dsRNAscan

**dsRNAscan** is a bioinformatics tool designed to identify and characterize **double-stranded RNAs (dsRNAs)** encoded by genomic sequences. This package also bundles a modified version of **einverted** from the EMBOSS suite for detecting inverted repeats that may form dsRNA structures.

---

## Installation & Running

### 1. Using `python` Directly

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

### 2. Making dsRNAscan.py Executable

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

### 3. Conda-Based Build (for Modified einverted + EMBOSS)

We provide a Conda-based recipe that compiles EMBOSS (including our modified einverted.c) and installs dsRNAscan:

1. Install conda-build:
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
./dsRNAscan.py --input <input_fasta_file> --chr <chromosome_name> [other flags...]
```

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

1. **Results File** (`<basename>_ein_results.txt`)
   
   Fields include:
   - Chromosome, Strand, Score, RawMatch, PercMatch, Gaps
   - i_start, i_end, j_start, j_end
   - i_seq, j_seq
   - structure (dot-bracket)
   - dG(kcal/mol)
   - percent_paired
   - link to FORNA visualiztion 
   - longest_helix length

2. **BP File** (`<basename>.dsRNApredictions.bp`)
   
   A tab-delimited file describing base pairs in a format suitable for IGV or other genome browsers.

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

This project is licensed under the MIT License.

---

## Contact

For questions or more information, contact:
Ryan Andrews ryan.andrews@biochem.utah.edu
