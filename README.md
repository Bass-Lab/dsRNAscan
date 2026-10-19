# dsRNAscan

[![CI Tests](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with support for **G-U wobble base pairing** and **ML-based scoring**.

Browse human genome results at [dsrna.chpc.utah.edu](http://dsrna.chpc.utah.edu)

## Install

```bash
pip install dsrnascan
```

Linux and macOS, Python 3.8+ (Windows: use WSL). If ViennaRNA fails via pip, install it with conda: `conda install -c bioconda viennarna`

## Usage

```bash
# Basic scan (defaults: -w 10000 -s 500 --score 75 -c 4)
dsrnascan input.fasta

# Scan a region
dsrnascan genome.fasta --only_seq chr21 --start 33455482 --end 33655482

# Faster fold (~99% identical results)
dsrnascan genome.fasta --banded -c 8

# Settings used for the human dsRNAome resource
dsrnascan genome.fasta --v1
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-w` | 10000 | Window size (bp) |
| `-s/--step` | 500 | Step size between windows |
| `-c/--cpus` | 4 | Number of CPUs |
| `--score` | 75 | Minimum einverted score |
| `--paired_cutoff` | 70 | Minimum % paired bases |
| `--only_seq`, `--start/--end` | Full sequence | Chromosome and region to scan |
| `--banded` | False | Faster fold (~99% identical results) |
| `--v1` | False | Settings used for the human dsRNAome resource |
| `--format` | bedpe | Output format: `bedpe`, `gff3`, or `both` |

Run `dsrnascan --help` for all options.

## Output

- **`*_merged_results.txt`**: one row per dsRNA with coordinates, einverted score, RNAduplex fold, structure and ML scores
- **`*.bedpe`**: paired arm coordinates (default)
- **`*.gff3`**: arms as linked features for IGV (`--format gff3` or `both`)
- **`*.bp`**: base-pair arcs for IGV

## Browse Results

```bash
# Interactive viewer with Forna RNA structure visualization
dsrna-browse results_directory/

# With RNA editing site annotations (BED or GFF3)
dsrna-browse results_directory/ --editing-file editing_sites.bed
```

## Citation

If you use dsRNAscan, please cite:

> Comprehensive mapping of human dsRNAome reveals conservation, neuronal enrichment, and intermolecular interactions
> https://doi.org/10.1101/2025.01.24.634786

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE). **Issues:** [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)

## Acknowledgments

- EMBOSS einverted (search) and ViennaRNA duplexfold (folding), both reimplemented in dsRNAscan's compiled engine
