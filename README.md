# dsRNAscan

**dsRNAscan** is a bioinformatics tool designed to identify and characterize **double-stranded RNAs (dsRNAs)** encoded by genomic sequences. It scans the genome for potential dsRNA-forming regions by predicting intra-molecular base pairing, evaluating various structural characteristics, and identifying key features relevant to dsRNA biology, such as RNA editing by ADAR enzymes.

## Requirements

- **Python 3.x**
- **Dependencies**:
  - `numpy`
  - `pandas`
  - `BioPython`
  - `ViennaRNA` (RNA structure prediction library)
  - `argparse`, `glob`, `multiprocessing`

To install these dependencies, run:

```bash
pip install numpy pandas biopython RNAhon RNA tensorflow
```

## Installation

Clone the repository:

```bash
git clone https://github.com/Bass-Lab/dsRNAscan.git
cd dsRNAscan
```

## Usage

Basic Command-Line Usage

```bash
python dsRNAscan.py --input <input_fasta_file> 
```

Key Parameters

	•	--input: Path to the input FASTA file containing genomic sequences.
	•	--chr: Chromosome name. If “header”, it uses the sequence header as the chromosome name.
	•	--start: Starting coordinate for scanning (default: 0).
	•	--end: Ending coordinate for scanning (default: 0, which sets the end to the length of the sequence).
	•	--step-size: Step size for moving the window across sequences (default: 150 bp).
	•	--window-size: Sliding window size for scanning genomic sequences (default: 10,000 bp).
	•	--min-length: Minimum length for dsRNA prediction (default: 30 bp).
	•	--max_span: Maximum span of inverted repeats (default: 10,000 bp).
	•	--gaps: Gap penalty for the inverted repeat finder (default: 12).
	•	--match: Match score for base pairing (default: 3).
	•	--mismatch: Mismatch penalty (default: -4).
	•	--paired_cutoff: Minimum percentage of paired bases to retain results (default: 70%).
	•	--algorithm: Algorithm for inverted repeat detection (default: “einverted”). Options include "einverted" or "iupacpal".
	•	--reverse: An optional flag to run the analysis on the reverse strand.
	•	--cpus: Number of CPUs to use for parallel processing (default: 4).

Example

```bash
python dsRNAscan.py --input genome.fasta --window-size 10000 --step-size 500
```

This command scans the given genome FASTA file with a window size of 10,000 bp, moving 500 bp at a time (this is also the default setting).

Output

	•	Results File (<basename>_ein_results.txt): Contains dsRNA predictions with the following columns:
	•	Chromosome: Chromosome of the predicted dsRNA.
	•	Strand: Strand information (+ or -).
	•	Score: Inverted repeat score.
	•	RawMatch: Number of matches vs. total bases.
	•	PercMatch: Percentage of base pairing.
	•	Gaps: Number of gaps in the sequence alignment.
	•	i_start, i_end: Start and end coordinates of the first dsRNA arm.
	•	j_start, j_end: Start and end coordinates of the second dsRNA arm.
	•	i_seq, j_seq: Sequences of the predicted dsRNA arms.
	•	structure: Predicted secondary structure in dot-bracket notation.
	•	dG(kcal/mol): Calculated free energy of the dsRNA structure.
	•	percent_paired: Percentage of paired bases in the predicted dsRNA.
 
	•	Base Pair Prediction File (<basename>.dsRNApredictions.bp): Provides base pair predictions in a tab-separated format.


Logging

The script generates a log file named process.log, detailing each processing step, potential errors, and runtime information. This is useful for tracking progress and troubleshooting.

Contributing

Contributions are welcome! Please fork this repository, make your changes, and submit a pull request. Ensure that your code is well-documented and tested.

License

This project is licensed under the MIT License.

Contact

For questions or further information, please contact Ryan Andrews at Ryan.andrews@biochem.utah.edu

