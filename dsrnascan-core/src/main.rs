//! Command-line scanner: prints dsRNAs as TSV (dsRNAscan's merged_results columns), or raw
//! search hits with --no-select.

use dsrnascan_core::{candidates, complement, dsrnas, dsrnas_both, einverted, encode, init_threads, par_map, round2, Params};
use std::fmt::Write as _;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;

const USAGE: &str = "usage: dsrnascan-core <fasta> [options]
  search:   --threshold N (75)  --match N (3)  --mismatch N (-4)  --gap N (12)  --no-gu
            -w/--max-span N (10000)  --min-match-perc F (70)
  windows:  -s/--step N (500)  --window N (max span)  --windows W1,W2,..
            --minus-revcomp (scan the minus strand on its reverse complement)
  filters:  --percent-paired F (70)  --min-bp N (25)
  fold:     --banded (faster, ~99% identical)  --band-margin N (30)  -t/--temperature C (37)
  presets:  --v1 (settings used for the human dsRNAome resource)
  region:   --only-seq NAME  --start N  --end N (1-based; output in full-sequence coordinates)
  other:    --strand plus|minus|both (both)  --threads N (all cores)  --block-size N (auto)
            --no-select (raw search hits only)
  experimental: --windowless  --greed F  --no-rescue  --arm-overlap F  --overlap-frac F  --compact";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut p = Params::default();
    let mut strand = String::from("both");
    let mut raw = false;
    let (mut only_seq, mut start, mut end) = (None::<String>, 0usize, 0usize);
    let mut path = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut val = || it.next().unwrap_or_else(|| die(&format!("{a} needs a value")));
        match a.as_str() {
            "--threshold" => p.threshold = num(val()),
            "--match" => p.match_score = num(val()),
            "--mismatch" => p.mismatch = num(val()),
            "--gap" => p.gap = num(val()),
            "-w" | "--max-span" => p.max_span = num(val()),
            "--window" => p.window = num(val()),
            "--min-match-perc" => p.min_match_perc = num(val()),
            "--percent-paired" => p.min_percent_paired = num(val()),
            "--min-bp" => p.min_bp = num(val()),
            "--overlap-frac" => p.overlap_frac = Some(num(val())),
            "--compact" => p.overlap_frac = Some(0.5),
            "--greed" => p.greed = num(val()),
            "--arm-overlap" => p.arm_frac = num(val()),
            "--no-rescue" => p.rescue = false,
            "--no-gu" => p.gu = false,
            "--minus-revcomp" => p.minus_revcomp = true,
            "--windows" => p.windows = val().split(',').map(|w| num(w.trim())).collect(),
            "--v1" => {
                // settings used for the human dsRNAome resource
                p.windows = vec![10_000, 3_000];
                p.window_step = Some(150);
                p.min_bp = 0;
                p.min_match_perc = 0.0;
                p.minus_revcomp = true;
                p.band_margin = Some(p.band_margin.unwrap_or(30)); // banded fold: 2.3x faster, 99.5% identical
            }
            "-s" | "--step" => p.window_step = Some(num(val())),
            "--windowless" => p.window_step = None,
            "--banded" => p.band_margin = Some(p.band_margin.unwrap_or(30)),
            "--band-margin" => p.band_margin = Some(num(val())),
            "-t" | "--temperature" => {
                p.temperature = num(val());
                if !(0.0..=100.0).contains(&p.temperature) {
                    die("--temperature must be between 0 and 100 C");
                }
            }
            "--no-select" => raw = true,
            "--only-seq" | "--only_seq" => only_seq = Some(val().clone()),
            "--start" => start = num(val()),
            "--end" => end = num(val()),
            "--strand" => strand = val().clone(),
            "--threads" => p.threads = num(val()),
            "--block-size" => p.block_size = num(val()),
            "-h" | "--help" => die(USAGE),
            s if s.starts_with('-') => die(&format!("unknown option {s}\n{USAGE}")),
            s => path = Some(s.to_string()),
        }
    }
    let path = path.unwrap_or_else(|| die(USAGE));
    let strands: &[char] = match strand.as_str() {
        "plus" => &['+'],
        "minus" => &['-'],
        "both" => &['+', '-'],
        s => die(&format!("bad --strand {s}")),
    };

    init_threads(p.threads);
    let t0 = Instant::now();
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    if raw {
        writeln!(out, "seq\tstrand\ti_start\ti_end\tj_start\tj_end\tscore\tmatches\ttotal\tmatch_perc\tgaps\ti_seq\tj_seq").unwrap();
    } else {
        writeln!(out, "Chromosome\tStrand\ti_start\ti_end\tj_start\tj_end\tScore\tRawMatch\tPercMatch\tGaps\t\
                       dG(kcal/mol)\tbase_pairs\tpercent_paired\tlongest_helix\teff_i_start\teff_i_end\t\
                       eff_j_start\teff_j_end\ti_seq\tj_seq\tstructure").unwrap();
    }
    // Every (sequence, strand) is an independent job on the shared pool; each job's own blocks,
    // windows and folds run on the same pool, so short and long sequences balance out. Output is
    // written in input order afterwards, so it does not depend on the number of threads.
    // (reported name, offset of the scanned sequence in that reference, scanned sequence,
    //  reference length when known: a samtools region record does not say how long its chromosome is)
    let mut records: Vec<(String, usize, Vec<u8>, Option<usize>)> = Vec::new();
    for (header, seq) in read_fasta(&path) {
        let region = region_header(&header, seq.len());
        let (name, off) = region.clone().unwrap_or((header.clone(), 0));
        if only_seq.as_ref().is_some_and(|o| *o != header && *o != name) {
            continue;
        }
        // 0.5.x: start_pos = max(0, start - 1); end 0 = to the end
        let (s0, e0) = (start.saturating_sub(1).min(seq.len()), if end > 0 { end.min(seq.len()) } else { seq.len() });
        if p.minus_revcomp && region.is_some() {
            eprintln!("note: {header}: chromosome length unknown, so windows start at this record's ends, \
                       not on the whole-chromosome grid a full run would use");
        }
        records.push((name, off + s0, seq[s0..e0.max(s0)].to_vec(), region.is_none().then_some(seq.len())));
    }
    let n_bases: usize = records.iter().map(|(_, _, s, _)| s.len()).sum();
    // both strands of a sequence share one fused sweep ('*' job) unless raw candidates are asked for
    let fused = !raw && strands.len() == 2;
    let lanes: &[char] = if fused { &['*'] } else { strands };
    let jobs: Vec<(usize, char)> = (0..records.len()).flat_map(|r| lanes.iter().map(move |&s| (r, s))).collect();
    let chunks: Vec<(usize, String)> = par_map(&jobs, |&(r, s)| {
        let (name, o, seq, len) = &records[r];
        let o = *o; // every coordinate is reported on the full reference
        // reverse-complement scans reproduce whole-chromosome runs, so their windows sit on that
        // run's grid even when only a region is scanned
        let p = Params { grid: if p.minus_revcomp { len.map(|l| (o, l)) } else { None }, ..p.clone() };
        let enc = if s == '-' { complement(&encode(seq)) } else { encode(seq) };
        let mut text = String::new();
        let mut n = 0;
        if raw {
            let letters: Vec<u8> = enc.iter().map(|&b| b"ACGTN"[b as usize]).collect();
            let hits = match p.window_step {
                Some(step) => einverted::window_hits(&enc, &p, step), // einverted's raw window hits
                None => candidates(&enc, &p),
            };
            for h in hits {
                let mut a = letters[h.i_start - 1..h.i_end].to_vec();
                let mut b = letters[h.j_start - 1..h.j_end].to_vec();
                if s == '-' {
                    a.reverse();
                    b.reverse();
                }
                writeln!(text, "{name}\t{s}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}",
                         h.i_start + o, h.i_end + o, h.j_start + o, h.j_end + o, h.score, h.matches,
                         h.matches + h.mismatches, h.match_perc(), h.gaps,
                         String::from_utf8_lossy(&a), String::from_utf8_lossy(&b)).unwrap();
                n += 1;
            }
        } else {
            let sets: Vec<(char, Vec<dsrnascan_core::Dsrna>)> = if s == '*' {
                let [a, b] = dsrnas_both(&enc, &p);
                vec![('+', a), ('-', b)]
            } else {
                vec![(s, dsrnas(&enc, s == '-', &p))]
            };
            for (s, d) in sets.into_iter().flat_map(|(s, v)| v.into_iter().map(move |d| (s, d))) {
                let h = &d.hit;
                let total = h.matches + h.mismatches;
                writeln!(text, "{name}\t{s}\t{}\t{}\t{}\t{}\t{}\t{}/{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                         h.i_start + o, h.i_end + o, h.j_start + o, h.j_end + o, h.score, h.matches, total,
                         py_float(round2(h.match_perc())), h.gaps, py_float(d.energy as f64 / 100.0),
                         d.base_pairs, py_float(d.percent_paired), d.longest_helix,
                         d.eff[0] + o, d.eff[1] + o, d.eff[2] + o, d.eff[3] + o, d.i_seq, d.j_seq, d.structure).unwrap();
                n += 1;
            }
        }
        (n, text)
    });
    let mut n_out = 0;
    for (n, text) in chunks {
        out.write_all(text.as_bytes()).unwrap();
        n_out += n;
    }
    out.flush().unwrap();
    eprintln!(
        "{n_out} {} | {n_bases} bp x {} strand(s) | {:.2} s | {} threads",
        if raw { "candidates" } else { "dsRNAs" },
        strands.len(),
        t0.elapsed().as_secs_f64(),
        p.threads
    );
}

/// Python's float repr for the values dsRNAscan writes: "80.0", "91.15", "-205.8".
fn py_float(x: f64) -> String {
    if x.fract() == 0.0 {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}

/// A samtools-style region name ("chr19:10000001-10300000", 1-based inclusive) whose length matches
/// the sequence: the chromosome name and the 0-based offset of the sequence in it.
fn region_header(name: &str, len: usize) -> Option<(String, usize)> {
    let (chrom, range) = name.rsplit_once(':')?;
    let (a, b) = range.split_once('-')?;
    let (a, b): (usize, usize) = (a.replace(',', "").parse().ok()?, b.replace(',', "").parse().ok()?);
    (a >= 1 && b >= a && b - a + 1 == len && !chrom.is_empty()).then(|| (chrom.to_string(), a - 1))
}

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
    let mut recs: Vec<(String, Vec<u8>)> = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if let Some(h) = line.strip_prefix('>') {
            recs.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(r) = recs.last_mut() {
            r.1.extend(line.trim_end().bytes());
        }
    }
    recs
}

fn num<T: std::str::FromStr>(s: &str) -> T {
    s.parse().unwrap_or_else(|_| die(&format!("not a number: {s}")))
}

fn die(msg: &str) -> ! {
    eprintln!("{msg}");
    std::process::exit(2)
}
