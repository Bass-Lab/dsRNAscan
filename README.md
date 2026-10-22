# dsRNAscan

[![CI Tests](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml/badge.svg)](https://github.com/rjandr/dsRNAscan/actions/workflows/ci-simple.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**dsRNAscan** is a bioinformatics tool for genome-wide identification of **double-stranded RNA (dsRNA) structures**. It uses a sliding window approach to detect inverted repeats that can form dsRNA secondary structures, with support for **G-U wobble base pairing** and **ML-based scoring**.

Browse human genome results at [dsrna.chpc.utah.edu](http://dsrna.chpc.utah.edu)

## Install

```bash
pip install dsrnascan
```

Linux and macOS, Python 3.8+ (Windows: use WSL).

## Usage

```bash
# Scan a genome or transcriptome
dsrnascan input.fasta

# Scan a region
dsrnascan genome.fasta --only_seq chr21 --start 33455482 --end 33655482
```

By default dsRNAscan uses the settings of the human dsRNAome resource (`--v1`, on by default): windows of
10000 and 3000 nt, step 150, and a banded fold that gives ~99% identical results to the exact one. Any
setting you change is listed in `dsrnascan.log`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-w` | 10000 and 3000 | Window size (one window of this size if set) |
| `-s/--step` | 150 | Step size between windows |
| `-c/--cpus` | 4 (whole allocation in a SLURM job) | Number of CPUs |
| `--score` | 75 | Minimum einverted score |
| `--percent_paired_cutoff` | 70 | Minimum % paired bases after folding |
| `--only_seq`, `--start/--end` | Full sequence | Sequence and region to scan |
| `--exact` | off | Exact fold (slower) |
| `--format` | bedpe | Output format: `bedpe`, `gff3`, or `both` |

Run `dsrnascan --help` for all options.

## Output

- **`*.dsrnas.parquet`**: one row per dsRNA with coordinates, search score, duplex fold, structure and ML scores
  (`pd.read_parquet` in Python, `arrow::read_parquet` in R; `--txt` writes the same table as tab-separated text)
- **`*.bedpe`**: paired arm coordinates (default)
- **`*.gff3`**: arms as linked features for IGV (`--format gff3` or `both`)
- **`*.bp`**: base-pair arcs for IGV

ML scores (`likely_edited`, `likely_forms`) come from structure-only models trained on human data, where
inverted Alus dominate. Treat them as a rough guide, especially outside Alus and in other species. Both scores run 0-1 (cutoffs: stability 0.553, probing 0.440).

## Browse, filter and export

```bash
# Interactive browser: filters, sorting, structure view, export of the filtered set
dsrna-browse results_folder/

# With RNA editing sites (BED or GFF3), counted per arm and shown on the structure
dsrna-browse results_folder/ --editing-file editing_sites.bed

# The same filters from the command line
dsrna-browse results_folder/ --min percent_paired 85 --max "dG(kcal/mol)" -150 --export strong.gff3
dsrna-browse results_folder/ --region chr19:1-3000000 --likely-edited --export region.tsv
```

Export formats: `.tsv`, `.parquet`, `.bedpe`, `.bed`, `.gff3`, `.bp`, `.fa`, `.dbn`. `--columns` lists the filterable
columns; `--where` takes an SQL condition on them.

## Citation

If you use dsRNAscan, please cite:

> dsRNAscan maps human dsRNAome, revealing conservation, intermolecular dsRNA, and correlates of ADAR dependency.
> *Molecular Cell* (2026). https://doi.org/10.1016/j.molcel.2026.08.001

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE). **Issues:** [GitHub Issues](https://github.com/Bass-Lab/dsRNAscan/issues)

## Acknowledgments

dsRNAscan's compiled engine reimplements:

- **EMBOSS einverted** (inverted-repeat search): Rice P, Longden I, Bleasby A (2000) *Trends Genet* 16:276-277
- **ViennaRNA duplexfold** (duplex folding): Lorenz R et al. (2011) *Algorithms Mol Biol* 6:26, with
  Turner 2004 energy parameters (Mathews DH et al. (2004) *PNAS* 101:7287)
