//! dsrnascan-core: find intramolecular dsRNAs in one sweep of the sequence (see SPEC.md).
//!
//! Per sequence and strand:
//!   1. `scan`   - one pass over pair space, block-parallel; best cell per path origin, with the
//!                 path's gap count (which fixes its match/mismatch counts); the match-percentage
//!                 cutoff constrains which cell is a path's best
//!   2. `select` - box rule, strongest first; each candidate that is not already covered by a kept
//!                 hit is folded (`duplex`, a port of RNA.duplexfold) and kept only if the fold
//!                 passes the percent-paired and base-pair filters

pub mod duplex;
pub mod einverted;
pub mod energy;
pub mod scan;
#[cfg(feature = "python")]
mod python;
pub mod select;
mod turner2004;

use rayon::prelude::*;

/// Search, filter and selection parameters. Defaults reproduce dsRNAscan 0.5.x's output
/// (-w 10000 -s 500 --score 75 --min_bp 25 --paired_cutoff 70; see `einverted.rs`).
#[derive(Clone, Debug)]
pub struct Params {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap: i32,
    pub threshold: i32,
    /// Same meaning as einverted's -maxrepeat (dsRNAscan's -w): outer 3' position minus outer 5'
    /// position is at most `max_span - 2`.
    pub max_span: usize,
    /// Minimum aligned match percentage (dsRNAscan --paired_cutoff), a search constraint.
    pub min_match_perc: f64,
    /// Minimum percent paired of the folded duplex (dsRNAscan --percent_paired_cutoff).
    pub min_percent_paired: f64,
    /// Minimum base pairs in the folded duplex (dsRNAscan --min_bp).
    pub min_bp: usize,
    /// Drop a hit when at least this fraction of BOTH its arms lies inside the matching arms of a
    /// stronger kept hit. 1.0 = drop only fully nested hits (0.5.x). `None` = no selection at all.
    pub overlap_frac: Option<f64>,
    /// Greed: drop a hit when one of its arms overlaps an arm of a stronger kept hit by at least
    /// `arm_frac` and it scores below `greed` x that hit's score (see select.rs). 0 = off.
    pub greed: f64,
    pub arm_frac: f64,
    /// Under greed, keep a weaker partner whose stronger rival's other arm lies outside its span,
    /// as 0.5.x's overlapping windows did.
    pub rescue: bool,
    /// Reproduce 0.5.x (default 500): einverted's reporting replayed on windows of `max_span`
    /// every this many nt (`einverted.rs`), then 0.5.x's fold, filters and nested elimination.
    /// `None` = windowless selection (`overlap_frac`, `greed`, `rescue`).
    pub window_step: Option<usize>,
    /// Window width for the 0.5.x replay (dsRNAscan -w); 0 = `max_span`. `max_span` is then
    /// einverted's -maxrepeat inside each window (0.5.x: min(-w, --max_span)).
    pub window: usize,
    /// Fold only the band around the search path, widened by this many nt (`None` = exact fold).
    pub band_margin: Option<usize>,
    /// Rows of pair space owned by each parallel block (0 = automatic).
    pub block_size: usize,
    /// Size of the global thread pool (set once by the caller, see `init_threads`).
    pub threads: usize,
    /// Search scores G-U as a match (dsRNAscan's einverted); false = stock EMBOSS einverted.
    pub gu: bool,
    /// Several window sizes, merged before folding (`--windows`, `--v1`); empty = just `window`.
    /// Each run uses einverted -maxrepeat = min(its window, max_span).
    pub windows: Vec<usize>,
    /// Scan the minus strand on the reverse complement, as dsRNAscan v0.1.0 did (0.5.x complements
    /// without reversing, which puts its windows at different places on the minus strand).
    pub minus_revcomp: bool,
    /// (offset of the scanned sequence in its chromosome, chromosome length): put windows where a
    /// whole-chromosome run puts them, every `window_step` from the chromosome's first nt (plus,
    /// complemented minus) or from its last nt (reverse-complement minus), so a region scan gives
    /// that run's hits. `None` = windows start at the scanned sequence's own ends (0.5.x).
    pub grid: Option<(usize, usize)>,
    /// Folding temperature in Celsius (energies rescaled as ViennaRNA does; 37 = Turner 2004 as published).
    pub temperature: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            match_score: 3,
            mismatch: -4,
            gap: 12,
            threshold: 75,
            max_span: 10_000,
            min_match_perc: 70.0,
            min_percent_paired: 70.0,
            min_bp: 25,
            overlap_frac: Some(1.0),
            greed: 1.0,
            arm_frac: 0.5,
            rescue: true,
            window_step: Some(500),
            window: 0,
            band_margin: None,
            block_size: 0,
            threads: std::thread::available_parallelism().map_or(1, |n| n.get()),
            gu: true,
            windows: Vec::new(),
            minus_revcomp: false,
            grid: None,
            temperature: 37.0,
        }
    }
}

/// One search hit (an inverted repeat). Coordinates are 1-based and inclusive on the scanned
/// sequence. Following dsRNAscan's naming, `i` is the lower-coordinate arm and `j` the upper.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hit {
    pub i_start: usize,
    pub i_end: usize,
    pub j_start: usize,
    pub j_end: usize,
    pub score: i32,
    pub matches: u32,
    pub mismatches: u32,
    pub gaps: u32,
}

impl Hit {
    /// einverted's match percentage: matches / (matches + mismatches); gaps not counted.
    pub fn match_perc(&self) -> f64 {
        let total = self.matches + self.mismatches;
        if total == 0 {
            0.0
        } else {
            100.0 * self.matches as f64 / total as f64
        }
    }
}

/// A reported dsRNA: the search hit plus its duplex fold, as dsRNAscan 0.5.x reports it.
#[derive(Clone, Debug)]
pub struct Dsrna {
    pub hit: Hit,
    /// duplex free energy, dcal/mol
    pub energy: i32,
    pub base_pairs: usize,
    pub percent_paired: f64,
    pub longest_helix: usize,
    /// eff_i_start, eff_i_end, eff_j_start, eff_j_end: the duplex-trimmed extent
    pub eff: [usize; 4],
    /// duplex-trimmed arms (RNA letters), each read 5'->3' along the scanned strand
    pub i_seq: String,
    pub j_seq: String,
    pub structure: String,
}

/// numpy.round(x, 2), which dsRNAscan applies before comparing against cutoffs.
pub fn round2(x: f64) -> f64 {
    (x * 100.0).round_ties_even() / 100.0
}

/// A=0 C=1 G=2 T/U=3, anything else 4 (pairs with nothing, like einverted).
pub fn encode(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => 4,
        })
        .collect()
}

/// Complement without reversing, so coordinates stay on the forward axis (as dsRNAscan does).
pub fn complement(enc: &[u8]) -> Vec<u8> {
    enc.iter().map(|&b| if b < 4 { 3 - b } else { 4 }).collect()
}

/// Search hits of one encoded sequence that meet the match constraint, unselected.
pub fn candidates(enc: &[u8], p: &Params) -> Vec<Hit> {
    let n = enc.len();
    let halo = p.max_span + 2;
    // A path can only reach rows < origin + halo, and only be reached from rows > row - halo,
    // so computing [own.start - halo, own.end + halo) makes every owned origin exact.
    par_map(&blocks(n, p), |&(s, e)| {
        scan::scan_block(enc, p, s.saturating_sub(halo)..(e + halo).min(n), s..e)
    })
    .into_iter()
    .flatten()
    .map(|c| c.to_hit(p))
    .collect()
}

/// The full pipeline for one strand: search, then selection with lazy folding and filters.
/// `minus`: `enc` is the complement of the forward sequence and arms are read 5'->3' along it.
pub fn dsrnas(enc: &[u8], minus: bool, p: &Params) -> Vec<Dsrna> {
    let mut out: Vec<Dsrna> = match p.window_step {
        Some(step) => {
            // 0.5.x: the windows' einverted hits, then match pre-filter + fold + filters, then
            // nested elimination
            return finish_windows(enc, strand_window_hits(enc, minus, p, step), minus, p);
        }
        None => {
            let rule = p.overlap_frac.map(|f| select::Rule {
                both_frac: f,
                arm_frac: p.arm_frac,
                greed: p.greed,
                rescue: p.rescue,
            });
            select::select_lazy(candidates(enc, p), rule, p.max_span, |h| fold(enc, h, minus, p))
                .into_iter()
                .map(|(_, d)| d)
                .collect()
        }
    };
    out.sort_by_key(|d| (d.hit.i_start, d.hit.j_end, d.hit.i_end, d.hit.j_start));
    out
}

/// Both strands of one sequence: [plus, minus]. With 0.5.x windows (the default) the two strands
/// share one fused sweep; otherwise each strand is scanned on its own.
pub fn dsrnas_both(enc: &[u8], p: &Params) -> [Vec<Dsrna>; 2] {
    let comp = complement(enc);
    let Some(step) = p.window_step else {
        return [dsrnas(enc, false, p), dsrnas(&comp, true, p)];
    };
    if p.windows.len() > 1 || p.minus_revcomp || p.grid.is_some() {
        let (hp, hm) = (strand_window_hits(enc, false, p, step), strand_window_hits(&comp, true, p, step));
        return [finish_windows(enc, hp, false, p), finish_windows(&comp, hm, true, p)];
    }
    let [hp, hm] = einverted::window_hits_both(&enc, &single_run(p, window_list(p)[0]), step);
    [finish_windows(enc, hp, false, p), finish_windows(&comp, hm, true, p)]
}

/// The window sizes to run: `windows`, or the single `window` (0 = `max_span`).
fn window_list(p: &Params) -> Vec<usize> {
    if p.windows.is_empty() {
        vec![if p.window > 0 { p.window } else { p.max_span }]
    } else {
        p.windows.clone()
    }
}

/// Parameters of one window-size run: einverted -maxrepeat = min(window, max_span), as 0.5.x.
fn single_run(p: &Params, w: usize) -> Params {
    Params { window: w, max_span: w.min(p.max_span), windows: Vec::new(), ..p.clone() }
}

/// einverted's window hits for one strand over every window size, merged. `enc` is the strand as
/// scanned by 0.5.x (complemented for the minus strand); with `minus_revcomp` the minus strand is
/// instead scanned on the reverse complement and its hits mirrored back onto the forward axis.
/// With `grid`, the windows start on the whole chromosome's grid instead of at the sequence's end.
fn strand_window_hits(enc: &[u8], minus: bool, p: &Params, step: usize) -> Vec<Hit> {
    let n = enc.len();
    let rc: Option<Vec<u8>> = (minus && p.minus_revcomp).then(|| enc.iter().rev().copied().collect());
    let seq = rc.as_deref().unwrap_or(enc);
    // first window start: distance from the grid's origin (chromosome start, or its end on the
    // reverse complement) to the scanned sequence, rounded up to the next grid point
    let phase = p.grid.map_or(0, |(off, len)| {
        let from = if rc.is_some() { len - (off + n) } else { off };
        (step - from % step) % step
    });
    if phase >= seq.len() {
        return Vec::new();
    }
    // every window size from one shared sweep (see einverted::window_hits_multi)
    let runs: Vec<(usize, Params)> = window_list(p).into_iter().map(|w| (w, single_run(p, w))).collect();
    let mut hits = einverted::window_hits_multi(&seq[phase..], &runs, step);
    for h in &mut hits {
        (h.i_start, h.i_end, h.j_start, h.j_end) = (h.i_start + phase, h.i_end + phase, h.j_start + phase, h.j_end + phase);
    }
    if rc.is_some() {
        // position x on the reverse complement is n - x + 1 forward; the arms swap roles
        for h in &mut hits {
            (h.i_start, h.i_end, h.j_start, h.j_end) = (n + 1 - h.j_end, n + 1 - h.j_start, n + 1 - h.i_end, n + 1 - h.i_start);
        }
        hits.sort_by_key(|h| (h.i_start, h.j_end, h.i_end, h.j_start, h.score));
    }
    hits
}

/// Both strands like `dsrnas_both`, but handing each strand's dsRNAs to `emit` `chunk` nt at a time
/// (plus strand first, by i_start), so only one stretch of folded results is held at once. Same
/// dsRNAs in the same order as `dsrnas_both`.
pub fn dsrnas_both_stream(enc: &[u8], p: &Params, chunk: usize, emit: &mut dyn FnMut(char, Vec<Dsrna>)) {
    let comp = complement(enc);
    let Some(step) = p.window_step else {
        emit('+', dsrnas(enc, false, p));
        emit('-', dsrnas(&comp, true, p));
        return;
    };
    let (hp, hm) = if p.windows.len() > 1 || p.minus_revcomp || p.grid.is_some() {
        (strand_window_hits(enc, false, p, step), strand_window_hits(&comp, true, p, step))
    } else {
        let [a, b] = einverted::window_hits_both(enc, &single_run(p, window_list(p)[0]), step);
        (a, b)
    };
    finish_windows_chunked(enc, hp, false, p, chunk, &mut |v| emit('+', v));
    finish_windows_chunked(&comp, hm, true, p, chunk, &mut |v| emit('-', v));
}

/// `finish_windows` one stretch of `chunk` nt (by i_start) at a time. Nested elimination still sees
/// every container: a container starts no later than the dsRNA it contains and ends at or after its
/// i arm, so the containers from earlier stretches are among the kept hits that reach this stretch.
fn finish_windows_chunked(enc: &[u8], mut hits: Vec<Hit>, minus: bool, p: &Params, chunk: usize, emit: &mut dyn FnMut(Vec<Dsrna>)) {
    hits.sort_by_key(|h| h.i_start); // stable, as nested_05x's sort
    let chunk = chunk.max(1);
    let mut tail: Vec<Hit> = Vec::new(); // folded hits of earlier stretches that may still contain later ones
    let mut k = 0;
    while k < hits.len() {
        let lo = hits[k].i_start / chunk * chunk;
        let hi = lo.saturating_add(chunk);
        let end = k + hits[k..].partition_point(|h| h.i_start < hi);
        let folded: Vec<Dsrna> = par_map(&hits[k..end], |h| fold(enc, h, minus, p)).into_iter().flatten().collect();
        let mut out = if p.overlap_frac.is_some() {
            tail.retain(|t| t.i_end >= lo);
            let all: Vec<&Hit> = tail.iter().chain(folded.iter().map(|d| &d.hit)).collect();
            let nested = select::nested_flags(&all)[tail.len()..].to_vec();
            tail.extend(folded.iter().map(|d| d.hit.clone()));
            folded.into_iter().zip(nested).filter(|(_, x)| !x).map(|(d, _)| d).collect()
        } else {
            folded
        };
        out.sort_by_key(|d| (d.hit.i_start, d.hit.j_end, d.hit.i_end, d.hit.j_start));
        emit(out);
        k = end;
    }
}

/// 0.5.x after the windows: match pre-filter + fold + filters, then nested elimination.
fn finish_windows(enc: &[u8], hits: Vec<Hit>, minus: bool, p: &Params) -> Vec<Dsrna> {
    let folded: Vec<Dsrna> = par_map(&hits, |h| fold(enc, h, minus, p)).into_iter().flatten().collect();
    let mut out = if p.overlap_frac.is_some() { select::nested_05x(folded) } else { folded };
    out.sort_by_key(|d| (d.hit.i_start, d.hit.j_end, d.hit.i_end, d.hit.j_start));
    out
}

/// One arm, read 5'->3' along the scanned strand (reversed on the minus strand).
fn arm(enc: &[u8], start: usize, end: usize, minus: bool) -> Vec<u8> {
    let mut a: Vec<u8> = enc[start - 1..end].iter().map(|&b| b"ACGTN"[b as usize]).collect();
    if minus {
        a.reverse();
    }
    a
}

/// Fold a hit's arms as dsRNAscan 0.5.x does with RNA.duplexfold; `None` if it fails the filters.
fn fold(enc: &[u8], h: &Hit, minus: bool, p: &Params) -> Option<Dsrna> {
    if round2(h.match_perc()) < p.min_match_perc {
        return None; // only reachable when the cutoff is applied after selection
    }
    let (a, b) = (arm(enc, h.i_start, h.i_end, minus), arm(enc, h.j_start, h.j_end, minus));
    let en = energy::energy(p.temperature);
    let d = match p.band_margin {
        None => duplex::duplexfold(en, &a, &b)?,
        Some(m) => {
            let (lo, hi) = duplex::path_band(a.len(), b.len(), h.gaps as usize, m);
            duplex::duplexfold_band(en, &a, &b, lo, hi)?
        }
    };
    let (s1, s2) = d.structure.split_once('&').unwrap();
    let (len1, len2) = (s1.len(), s2.len());
    let percent_paired = round2(d.pairs as f64 * 2.0 / (len1 + len2) as f64 * 100.0);
    if percent_paired < p.min_percent_paired || d.pairs < p.min_bp {
        return None;
    }
    // trim positions within the arms, 1-based (dsRNAscan's i_trim_*/j_trim_*)
    let (its, ite, jts, jte) = (d.i + 1 - len1, d.i, d.j, d.j + len2 - 1);
    let eff = if minus {
        // arms were reversed, so map trim positions back from the arm ends
        [h.i_end + 1 - ite, h.i_end + 1 - its, h.j_end + 1 - jte, h.j_end + 1 - jts]
    } else {
        [h.i_start + its - 1, h.i_start + ite - 1, h.j_start + jts - 1, h.j_start + jte - 1]
    };
    let rna = |s: &[u8]| s.iter().map(|&c| if c == b'T' { 'U' } else { c as char }).collect::<String>();
    Some(Dsrna {
        hit: h.clone(),
        energy: d.energy,
        base_pairs: d.pairs,
        percent_paired,
        longest_helix: d.longest_helix,
        eff,
        i_seq: rna(&a[its - 1..ite]),
        j_seq: rna(&b[jts - 1..jte]),
        structure: d.structure,
    })
}

/// Order-preserving parallel map on the global rayon pool. Results do not depend on the number of
/// threads. Every item is its own task, so a run of slow items (huge gappy folds) spreads over the
/// pool instead of pinning one thread.
pub fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    items.par_iter().with_max_len(1).map(f).collect()
}

/// Size the global thread pool once, before any work (0 = all cores).
pub fn init_threads(n: usize) {
    let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
}

/// Sweep blocks for a sequence of length n: `p.block_size`, or automatically about two blocks per
/// thread, at least four halos long (a block also sweeps one halo of rows before it) and at most
/// 4 Mb, so a short sequence still uses every core and a chromosome pays well under 1% for halos.
pub fn blocks(n: usize, p: &Params) -> Vec<(usize, usize)> {
    let lo = 4 * (p.max_span + 2);
    let bs = match p.block_size {
        0 => n.div_ceil(2 * rayon::current_num_threads()).clamp(lo, lo.max(4_000_000)),
        b => b,
    };
    (0..n).step_by(bs.max(1)).map(|s| (s, (s + bs).min(n))).collect()
}
