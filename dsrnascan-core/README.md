# dsrnascan-core

Compiled scanning engine for [dsRNAscan](https://github.com/rjandr/dsRNAscan). It sweeps each
sequence once and reproduces dsRNAscan 0.5.x's output exactly (einverted's per-window reporting,
RNAduplex folds, filters and nested elimination), several times faster and in bounded memory.
Design and validation: [SPEC.md](SPEC.md), `validation/`.

## Install

It is built into dsRNAscan (`pip install dsrnascan`) as `dsrnascan._core` and used for every scan
at 37 C. For development, from the repository root: `pip install maturin && maturin develop --release`.

The einverted pipeline remains available as a legacy path: `dsrnascan --legacy-einverted`, and
automatically for `-t` other than 37. `--banded` folds each dsRNA only near its alignment (faster;
~1 in 500 dsRNAs gets slightly different fold values).

## Standalone CLI (development and benchmarking)

    cargo build --release
    ./target/release/dsrnascan-core genome.fa > merged_results.tsv      # --help for options

Credits: duplex folding is a port of ViennaRNA's `duplexfold` (Lorenz et al. 2011) with Turner 2004
parameters from NNDB (Mathews et al. 2004; Turner & Mathews 2010). Search recurrence and reporting
follow EMBOSS einverted (Rice et al. 2000).
