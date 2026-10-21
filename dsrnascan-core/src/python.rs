//! Python binding (`--features python`, built into the dsrnascan package by maturin): `dsrnascan._core.scan`.
//! dsRNAscan's Python side reads the FASTA and applies --only_seq/--start/--end, calls `scan`, and
//! keeps writing its files and ML scores from the rows returned here.

use crate::{complement, dsrnas, dsrnas_both, dsrnas_both_stream, encode, init_threads, par_map, round2, Dsrna, Params};
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// One dsRNA in dsRNAscan's merged_results column order.
type Row = (
    String, String, usize, usize, usize, usize, i32, String, f64, u32, f64, usize, f64, usize,
    usize, usize, usize, usize, String, String, String,
);

/// Scan `records` = [(name, offset, sequence)] on `strands` ("both", "plus" or "minus").
/// Coordinates are 1-based on the full reference: `offset` is where `sequence` starts in it.
/// `lengths`: each record's full chromosome length, so `minus_revcomp` region scans put their
/// windows where a whole-chromosome run would (optional; without it windows start at the ends).
/// `temperature`: folding temperature in Celsius.
/// `on_rows`, `chunk_bp`: if given, results are not returned but passed to `on_rows({column: list})`
/// one stretch of `chunk_bp` nt at a time (sequences in turn, plus strand then minus), so memory
/// holds one stretch of results; same rows in the same order as the returned columns.
/// Returns {column: list} with dsRNAscan's merged_results columns (pandas.DataFrame-ready).
#[pyfunction]
#[pyo3(signature = (records, strands="both", threads=0, window=10000, max_span=10000, step=500,
                    threshold=75, match_score=3, mismatch=-4, gap=12, min_match_perc=70.0,
                    min_percent_paired=70.0, min_bp=25, eliminate_nested=true, banded=false,
                    windows=Vec::new(), minus_revcomp=false, lengths=Vec::new(), temperature=37.0,
                    on_rows=None, chunk_bp=5_000_000))]
#[allow(clippy::too_many_arguments)]
fn scan<'py>(
    py: Python<'py>,
    records: Vec<(String, usize, Vec<u8>)>,
    strands: &str,
    threads: usize,
    window: usize,
    max_span: usize,
    step: usize,
    threshold: i32,
    match_score: i32,
    mismatch: i32,
    gap: i32,
    min_match_perc: f64,
    min_percent_paired: f64,
    min_bp: usize,
    eliminate_nested: bool,
    banded: bool,
    windows: Vec<usize>,
    minus_revcomp: bool,
    lengths: Vec<usize>,
    temperature: f64,
    on_rows: Option<PyObject>,
    chunk_bp: usize,
) -> PyResult<Bound<'py, PyDict>> {
    if !(0.0..=100.0).contains(&temperature) {
        return Err(pyo3::exceptions::PyValueError::new_err("temperature must be between 0 and 100 C"));
    }
    let strands: Vec<char> = match strands {
        "both" => vec!['+', '-'],
        "plus" => vec!['+'],
        "minus" => vec!['-'],
        s => return Err(pyo3::exceptions::PyValueError::new_err(format!("strands must be both/plus/minus, not {s}"))),
    };
    let p = Params {
        match_score,
        mismatch,
        gap,
        threshold,
        max_span,
        min_match_perc,
        min_percent_paired,
        min_bp,
        overlap_frac: eliminate_nested.then_some(1.0),
        band_margin: banded.then_some(30),
        window_step: Some(step),
        window,
        windows,
        minus_revcomp,
        temperature,
        threads: if threads > 0 { threads } else { Params::default().threads },
        ..Params::default()
    };
    init_threads(p.threads);
    if let Some(cb) = on_rows {
        // stream: one sequence at a time (each runs on all threads), results passed on per stretch
        let mut err: Option<PyErr> = None;
        py.allow_threads(|| {
            for (r, (name, o, seq)) in records.iter().enumerate() {
                let p = Params { grid: if p.minus_revcomp { lengths.get(r).map(|&l| (*o, l)) } else { None }, ..p.clone() };
                let mut emit = |s: char, v: Vec<Dsrna>| {
                    if err.is_some() || v.is_empty() {
                        return;
                    }
                    let rows = rows_of(name, *o, s, v);
                    Python::with_gil(|py| {
                        if let Err(e) = columns(py, &rows).and_then(|d| cb.call1(py, (d,))) {
                            err = Some(e);
                        }
                    });
                };
                if strands.len() == 2 {
                    dsrnas_both_stream(&encode(seq), &p, chunk_bp, &mut emit);
                } else {
                    let enc = if strands[0] == '-' { complement(&encode(seq)) } else { encode(seq) };
                    emit(strands[0], dsrnas(&enc, strands[0] == '-', &p));
                }
                if err.is_some() {
                    break;
                }
            }
        });
        return match err {
            Some(e) => Err(e),
            None => columns(py, &[]),
        };
    }
    // both strands of a sequence run as one job (fused into one sweep in the default mode)
    let lanes: Vec<char> = if strands.len() == 2 { vec!['*'] } else { strands };
    let jobs: Vec<(usize, char)> = (0..records.len()).flat_map(|r| lanes.iter().map(move |&s| (r, s))).collect();
    // the scan runs on the core's own threads; Python can do other work meanwhile
    let rows = py.allow_threads(|| {
        par_map(&jobs, |&(r, s)| {
            let (name, o, seq) = &records[r];
            // with minus_revcomp, windows sit on the whole chromosome's grid (`Params::grid`)
            let p = Params { grid: if p.minus_revcomp { lengths.get(r).map(|&l| (*o, l)) } else { None }, ..p.clone() };
            let enc = if s == '-' { complement(&encode(seq)) } else { encode(seq) };
            let sets: Vec<(char, Vec<Dsrna>)> = if s == '*' {
                let [a, b] = dsrnas_both(&enc, &p);
                vec![('+', a), ('-', b)]
            } else {
                vec![(s, dsrnas(&enc, s == '-', &p))]
            };
            sets.into_iter().flat_map(|(s, v)| rows_of(name, *o, s, v)).collect::<Vec<Row>>()
        })
    });
    let rows: Vec<Row> = rows.into_iter().flatten().collect();
    columns(py, &rows)
}

/// A sequence's dsRNAs on one strand as rows; `o` shifts coordinates to the full reference.
fn rows_of(name: &str, o: usize, s: char, v: Vec<Dsrna>) -> Vec<Row> {
    v.into_iter()
        .map(|d| {
            let h = &d.hit;
            (
                name.to_string(), s.to_string(),
                h.i_start + o, h.i_end + o, h.j_start + o, h.j_end + o,
                h.score, format!("{}/{}", h.matches, h.matches + h.mismatches),
                round2(h.match_perc()), h.gaps, d.energy as f64 / 100.0,
                d.base_pairs, d.percent_paired, d.longest_helix,
                d.eff[0] + o, d.eff[1] + o, d.eff[2] + o, d.eff[3] + o,
                d.i_seq, d.j_seq, d.structure,
            )
        })
        .collect()
}

/// Rows as {column: list} in dsRNAscan's column order (pandas.DataFrame-ready).
fn columns<'py>(py: Python<'py>, rows: &[Row]) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new_bound(py);
    macro_rules! col {
        ($name:expr, $i:tt) => {
            out.set_item($name, rows.iter().map(|r| r.$i.clone()).collect::<Vec<_>>())?
        };
    }
    col!("Chromosome", 0); col!("Strand", 1); col!("i_start", 2); col!("i_end", 3);
    col!("j_start", 4); col!("j_end", 5); col!("Score", 6); col!("RawMatch", 7);
    col!("PercMatch", 8); col!("Gaps", 9); col!("dG(kcal/mol)", 10); col!("base_pairs", 11);
    col!("percent_paired", 12); col!("longest_helix", 13); col!("eff_i_start", 14);
    col!("eff_i_end", 15); col!("eff_j_start", 16); col!("eff_j_end", 17); col!("i_seq", 18);
    col!("j_seq", 19); col!("structure", 20);
    Ok(out)
}

#[pymodule]
#[pyo3(name = "_core")]
fn dsrnascan_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(scan, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
