//! Development tool: read "strand1<TAB>strand2[<TAB>gaps<TAB>margin]" lines on stdin, write
//! "energy<TAB>structure<TAB>i<TAB>j<TAB>pairs<TAB>longest_helix" (energy in kcal/mol), one per line.
//! With gaps and margin, folds only the band around the search path (`path_band`).
//! Used by validation/duplex_parity.py to compare against RNA.duplexfold.

use dsrnascan_core::duplex::{duplexfold, duplexfold_band, path_band};
use std::io::{BufRead, Write};

fn main() {
    let stdin = std::io::stdin();
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        let (a, b) = (f[0].as_bytes(), f[1].as_bytes());
        let d = if f.len() >= 4 {
            let (gaps, margin) = (f[2].parse().unwrap(), f[3].parse().unwrap());
            let (lo, hi) = path_band(a.len(), b.len(), gaps, margin);
            duplexfold_band(a, b, lo, hi)
        } else {
            duplexfold(a, b)
        };
        match d {
            Some(d) => writeln!(
                out,
                "{:.2}\t{}\t{}\t{}\t{}\t{}",
                d.energy as f64 / 100.0, d.structure, d.i, d.j, d.pairs, d.longest_helix
            ),
            None => writeln!(out, "NA\tNA\tNA\tNA\t0\t0"),
        }
        .unwrap();
    }
}
