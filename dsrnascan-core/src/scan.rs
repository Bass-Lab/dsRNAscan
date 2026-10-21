//! The sweep over pair space.
//!
//! Row `q` is a 3'-arm position; column `k` pairs it with 5'-arm position `p = q - 1 - k`.
//! The recurrence is einverted's (Smith-Waterman, linear gap, G-U counts as a match), with
//! its tie rules: 5' bulge beats 3' bulge beats extending the helix. Unlike einverted, we keep
//! the best cell of EVERY path origin, not one per position, so two hits sharing an arm
//! both survive to selection.
//!
//! Each cell carries its path's origin and gap count. Because the forward tie rules are the
//! same as einverted's traceback rules, following the chosen predecessors from a best cell is
//! exactly einverted's traceback, so origin + gap count + score determine the whole hit and
//! no traceback pass is needed.

use crate::{Hit, Params};
use std::collections::HashMap;
use std::ops::Range;

/// Best cell of one path; `origin` encodes the innermost pair as `q_in * stride + k_in`.
#[derive(Clone, Copy, Debug)]
pub struct Candidate {
    pub score: i32,
    pub q_out: usize,
    pub k_out: usize,
    pub origin: u64,
    pub gaps: u32,
    /// the score includes one match einverted does not report (an adjacent-base pair; see COL0_META)
    pub hidden: bool,
}

impl Candidate {
    /// einverted's preference: higher score, then the later row, then the first column in a row.
    fn beats(&self, b: &Candidate) -> bool {
        (self.score, self.q_out, std::cmp::Reverse(self.k_out))
            > (b.score, b.q_out, std::cmp::Reverse(b.k_out))
    }

    /// Recover the hit from the box corners, score and gap count.
    pub fn to_hit(&self, p: &Params) -> Hit {
        let stride = columns(p) as u64 + 1;
        let (q_in, k_in) = ((self.origin / stride) as i64, (self.origin % stride) as i64);
        let (q_out, k_out) = (self.q_out as i64, self.k_out as i64);
        let (p_in, p_out) = (q_in - 1 - k_in, q_out - 1 - k_out);
        let (len5, len3, g) = (p_in - p_out + 1, q_out - q_in + 1, self.gaps as i64);
        // 5' bulges use only 5' bases, 3' bulges only 3' bases, helix steps one of each.
        let g5 = (g + len5 - len3) / 2;
        let steps = len5 - g5;
        let (m, mm) = (p.match_score as i64, p.mismatch as i64);
        let counted = self.score as i64 - if self.hidden { m } else { 0 };
        let matches = (counted + p.gap as i64 * g - steps * mm) / (m - mm);
        debug_assert_eq!(matches * m + (steps - matches) * mm - p.gap as i64 * g, counted);
        Hit {
            i_start: (p_out + 1) as usize,
            i_end: (p_in + 1) as usize,
            j_start: (q_in + 1) as usize,
            j_end: (q_out + 1) as usize,
            score: self.score,
            matches: matches as u32,
            mismatches: (steps - matches) as u32,
            gaps: self.gaps,
        }
    }
}

pub(crate) fn pair_table(p: &Params) -> [[i32; 5]; 5] {
    let mut t = [[p.mismatch; 5]; 5];
    // A-U, U-A, C-G, G-C, then G-U and U-G (dsRNAscan's patched einverted; stock EMBOSS omits them)
    let pairs: &[(usize, usize)] = &[(0, 3), (3, 0), (1, 2), (2, 1), (2, 3), (3, 2)];
    for &(a, b) in &pairs[..if p.gu { 6 } else { 4 }] {
        t[a][b] = p.match_score;
    }
    t
}

/// Number of columns per row (einverted: -maxrepeat minus its two sentinel cells).
pub(crate) fn columns(p: &Params) -> usize {
    p.max_span.saturating_sub(2)
}

// Per-cell path metadata packed in one u64:
//   rows since the origin (21 bits) | origin column (21) | hidden-match flag (1) | gap count (21).
// Rows-since-origin (not the absolute row) keeps it small: +1 whenever a step comes from the row above.
const GAP_BITS: u32 = 21;
const HIDDEN: u64 = 1 << GAP_BITS;
const K_BITS: u32 = 21;
const K_SHIFT: u32 = GAP_BITS + 1;
const ROW_ONE: u64 = 1 << (K_SHIFT + K_BITS);
const GAP_MASK: u64 = (1 << GAP_BITS) - 1;
const K_MASK: u64 = (1 << K_BITS) - 1;

// einverted's traceback never visits column 0 (adjacent bases, `while (max > 0 && j >= 1)`), so a
// path whose innermost pair is there is reported from the next cell outward, its score still
// including the column-0 match. Out of column 0 a path goes either by a helix step into column 2,
// for which a column-0 cell carries metadata that becomes "fresh origin at column 2, one hidden
// match" once the step adds ROW_ONE, or by a bulge into column 1, handled where column 1 is
// computed.
const COL0_META: u64 = ((2u64 << K_SHIFT) | HIDDEN).wrapping_sub(ROW_ONE);

/// Largest max_span the packing supports (rows-since-origin and origin column both < 2^21).
pub const MAX_SPAN_LIMIT: usize = 1 << K_BITS;

/// What the sweep does with a cell at or above threshold. Cells arrive row by row (q ascending),
/// columns ascending within a row.
pub(crate) trait Sink {
    fn cell(&mut self, q: usize, k: usize, score: i32, meta: u64);
}

/// Origin of a cell's path, from its packed metadata: (q_in, k_in, gaps, hidden match).
pub(crate) fn unpack(q: usize, meta: u64) -> (usize, usize, u32, bool) {
    (q - (meta / ROW_ONE) as usize, ((meta >> K_SHIFT) & K_MASK) as usize, (meta & GAP_MASK) as u32, meta & HIDDEN != 0)
}

/// Sweep `rows`; return the best cell of every path whose origin row lies in `own`
/// and whose best score reaches the threshold.
pub fn scan_block(seq: &[u8], p: &Params, rows: Range<usize>, own: Range<usize>) -> Vec<Candidate> {
    let stride = columns(p) as u64 + 1;
    let mut sink = PathSink { p, stride, held: None, best: HashMap::new() };
    sweep(seq, p, rows, &mut sink);
    if let Some(h) = sink.held.take() {
        flush(h, &mut sink.best);
    }
    let mut out: Vec<Candidate> = sink
        .best
        .into_values()
        .filter(|c| own.contains(&((c.origin / stride) as usize)))
        .collect();
    out.sort_by_key(|c| (c.q_out, c.k_out, c.origin));
    out
}

/// One step of a row's running maximum: the first cell, scanning partners outward, that beats
/// every cell before it (and the threshold). This is einverted's per-row `maxval` bookkeeping,
/// kept for every prefix so it can be cut at any window's first partner.
#[derive(Clone, Copy, Debug)]
pub struct Step {
    pub k: u32,
    pub score: i32,
    pub meta: u64,
}

/// Running-maximum steps of a contiguous range of rows: 16 bytes per step plus 4 per row.
pub struct RowSteps {
    first: usize,    // first row covered
    start: Vec<u32>, // row first + r has steps[start[r]..start[r + 1]]; one entry per row + an end marker
    steps: Vec<Step>,
}

impl RowSteps {
    /// Row q's best cell among columns k <= limit (first occurrence on ties), if any.
    pub fn best(&self, q: usize, limit: usize) -> Option<Step> {
        let r = q.checked_sub(self.first)?;
        let (a, b) = (*self.start.get(r)? as usize, *self.start.get(r + 1)? as usize);
        let row = &self.steps[a..b];
        let n = row.partition_point(|s| s.k as usize <= limit);
        n.checked_sub(1).map(|i| row[i])
    }
}

/// Sweep `rows` and record the running-maximum steps of the rows in `own` (`own` within `rows`),
/// kept only where a window starting on the `grid` (window starts every `grid` nt) can look them up.
/// `caps`: extra column limits to keep answers for (smaller runs sharing this sweep, whose own
/// sweeps would stop at column cap; see `einverted::window_hits_multi`).
pub fn row_steps(seq: &[u8], p: &Params, rows: Range<usize>, own: Range<usize>, grid: usize, window: usize, caps: &[usize]) -> RowSteps {
    let out = RowSteps { first: own.start, start: Vec::with_capacity(own.len() + 1), steps: Vec::new() };
    let mut sink = StepSink { thr: p.threshold, own: own.clone(), max: 0, grid: grid.max(1), w: window, caps: caps.to_vec(), out };
    sweep(seq, p, rows, &mut sink);
    sink.compress_last();
    sink.open_rows_to(own.end); // trailing rows with no cell at or above threshold: no steps
    let end = sink.out.steps.len() as u32;
    sink.out.start.push(end);
    sink.out
}

/// `row_steps` for both strands in one sweep of the plus-strand sequence: [plus, minus].
pub fn row_steps_both(seq: &[u8], p: &Params, rows: Range<usize>, own: Range<usize>, grid: usize, window: usize) -> [RowSteps; 2] {
    let mk = || StepSink {
        thr: p.threshold,
        own: own.clone(),
        max: 0,
        grid: grid.max(1),
        w: window,
        caps: Vec::new(),
        out: RowSteps { first: own.start, start: Vec::with_capacity(own.len() + 1), steps: Vec::new() },
    };
    let mut sinks = [mk(), mk()];
    sweep_lanes(seq, p, [pair_table(p), minus_pair_table(p)], rows, &mut sinks);
    sinks.map(|mut sink| {
        sink.compress_last();
        sink.open_rows_to(own.end);
        let end = sink.out.steps.len() as u32;
        sink.out.start.push(end);
        sink.out
    })
}

struct StepSink {
    thr: i32,
    own: Range<usize>,
    max: i32,
    grid: usize,
    w: usize,
    caps: Vec<usize>,
    out: RowSteps,
}

impl StepSink {
    /// Open every row before `q` that is not open yet; each starts at the current step count, so
    /// rows that never get a cell stay empty.
    fn open_rows_to(&mut self, q: usize) {
        self.compress_last();
        let n = self.out.steps.len() as u32;
        while self.out.first + self.out.start.len() < q {
            self.out.start.push(n);
        }
    }
}

impl StepSink {
    /// Keep only the steps of the last opened row that some window can ask for. A window starting
    /// at ws asks row q for its best cell with k <= q - ws - 1, and window starts lie on the grid,
    /// so a row needs at most w / grid entries however long its running maximum is.
    fn compress_last(&mut self) {
        let Some(&a) = self.out.start.last() else { return };
        let a = a as usize;
        if self.out.steps.len() - a <= 1 {
            return;
        }
        let q = self.out.first + self.out.start.len() - 1;
        let row = self.out.steps.split_off(a);
        // the column limits windows can ask row q about, plus the caps of runs sharing the sweep
        let mut limits: Vec<usize> = self.caps.clone();
        let mut ws = q / self.grid * self.grid; // latest window start at or before q
        loop {
            if q > ws {
                limits.push(q - ws - 1);
            }
            if ws < self.grid || ws - self.grid + self.w <= q {
                break; // no earlier window contains row q
            }
            ws -= self.grid;
        }
        limits.sort_unstable();
        let mut kept = usize::MAX;
        for lim in limits {
            let n = row.partition_point(|s| s.k as usize <= lim);
            if n > 0 && n - 1 != kept {
                kept = n - 1;
                self.out.steps.push(row[kept]);
            }
        }
    }
}

impl Sink for StepSink {
    fn cell(&mut self, q: usize, k: usize, score: i32, meta: u64) {
        if !self.own.contains(&q) {
            return;
        }
        if self.out.first + self.out.start.len() <= q {
            self.open_rows_to(q + 1); // row q's first cell at or above threshold
            self.max = self.thr - 1;
        }
        if score > self.max {
            self.max = score;
            self.out.steps.push(Step { k: k as u32, score, meta });
        }
    }
}

struct PathSink<'a> {
    p: &'a Params,
    stride: u64,
    held: Option<Candidate>, // one-entry cache: cells of one path arrive together
    best: HashMap<u64, Candidate>,
}

impl Sink for PathSink<'_> {
    fn cell(&mut self, q: usize, k: usize, score: i32, meta: u64) {
        record(self.p, self.stride, q, k, score, meta, &mut self.held, &mut self.best);
    }
}

/// The DP over `rows` of pair space, calling `sink` for every cell at or above threshold.
fn sweep<S: Sink>(seq: &[u8], p: &Params, rows: Range<usize>, sink: &mut S) {
    sweep_lanes(seq, p, [pair_table(p)], rows, std::array::from_mut(sink));
}

/// Pair scores of the minus strand read on the plus-strand sequence: the minus strand is scanned
/// as the complement, and complementing keeps every Watson-Crick pair but turns G-U into C-A.
pub(crate) fn minus_pair_table(p: &Params) -> [[i32; 5]; 5] {
    let plus = pair_table(p);
    let comp = |b: usize| if b < 4 { 3 - b } else { 4 };
    std::array::from_fn(|x| std::array::from_fn(|y| plus[comp(x)][comp(y)]))
}

/// `L` DPs over the same sequence in one pass, one per score table (lane): with tables
/// [plus, minus] the two strands share every sequence read and loop step, and their identical
/// per-cell arithmetic is laid out for the compiler to vectorise. Each lane computes exactly what
/// a single-table sweep would.
fn sweep_lanes<S: Sink, const L: usize>(
    seq: &[u8],
    p: &Params,
    tables: [[[i32; 5]; 5]; L],
    rows: Range<usize>,
    sinks: &mut [S; L],
) {
    let kk = columns(p);
    if kk < 2 {
        return;
    }
    assert!(p.max_span < MAX_SPAN_LIMIT, "max_span must be < {MAX_SPAN_LIMIT}");
    let n = seq.len();
    let rows = rows.start.min(n)..rows.end.min(n);
    // Each row reads its 5' partners as one forward slice of a reversed copy. The copy covers only
    // these rows and the kk bases before them (one per task in flight, so never the whole sequence):
    // rev[end - q + k] == seq[q - 1 - k].
    let (base, end) = (rows.start.saturating_sub(kk), rows.end);
    let rev: Vec<u8> = seq[base..end].iter().rev().copied().collect();
    let (gap, thr) = (p.gap, p.threshold);
    let (mut prev, mut cur) = (vec![[0i32; L]; kk], vec![[0i32; L]; kk]);
    let (mut pmeta, mut cmeta) = (vec![[0u64; L]; kk], vec![[0u64; L]; kk]);

    for q in rows {
        let kmax = q.min(kk); // valid k: 0..kmax (p >= 0)
        if kmax == 0 {
            prev.fill([0; L]);
            continue;
        }
        let row: [&[i32; 5]; L] = std::array::from_fn(|l| &tables[l][seq[q] as usize]);
        let partners = &rev[end - q..end - q + kmax];
        cur[0] = std::array::from_fn(|l| row[l][partners[0] as usize]);
        cmeta[0] = [COL0_META; L];
        if kmax < 2 {
            std::mem::swap(&mut prev, &mut cur);
            std::mem::swap(&mut pmeta, &mut cmeta);
            continue;
        }
        // Column 1 has no helix step (its diagonal would pair a base with itself). A bulge out of
        // column 0 makes this cell the reported origin with the column-0 match hidden; left and
        // up give the same metadata, so einverted's tie order does not matter here.
        let (mut c, mut m) = ([0i32; L], [0u64; L]);
        for l in 0..L {
            let s = row[l][partners[1] as usize];
            c[l] = cur[0][l].max(prev[0][l]) - gap;
            m[l] = (1u64 << K_SHIFT) | HIDDEN | 1;
            if s > c[l] {
                c[l] = s;
                m[l] = 1u64 << K_SHIFT;
            }
        }
        cur[1] = c;
        cmeta[1] = m;
        for l in 0..L {
            if c[l] >= thr {
                offer(&mut sinks[l], q, 1, c[l], m[l]);
            }
        }
        // The cell to the left (k - 1) is (c, m). Every choice below is a plain select, no branches.
        for k in 2..kmax {
            // SAFETY: 2 <= k < kmax <= kk == len of every row array, and partners.len() == kmax.
            let (b, up, mup, pd, mdp) = unsafe {
                (
                    *partners.get_unchecked(k) as usize,
                    *prev.get_unchecked(k - 1),
                    *pmeta.get_unchecked(k - 1),
                    *prev.get_unchecked(k - 2),
                    *pmeta.get_unchecked(k - 2),
                )
            };
            let mut hit = false;
            for l in 0..L {
                let s = row[l][b];
                let take_up = up[l] > c[l]; // tie keeps the 5' bulge, as einverted does
                c[l] = if take_up { up[l] } else { c[l] } - gap;
                m[l] = if take_up { mup[l] + ROW_ONE } else { m[l] } + 1;
                let ext = pd[l] > 0;
                let d = s + if ext { pd[l] } else { 0 };
                let md = if ext { mdp[l].wrapping_add(ROW_ONE) } else { (k as u64) << K_SHIFT }; // wraps only for COL0_META
                let take_d = d > c[l]; // tie keeps the bulge, as einverted does
                c[l] = if take_d { d } else { c[l] };
                m[l] = if take_d { md } else { m[l] };
                hit |= c[l] >= thr;
            }
            unsafe {
                *cur.get_unchecked_mut(k) = c;
                *cmeta.get_unchecked_mut(k) = m;
            }
            if hit {
                for l in 0..L {
                    if c[l] >= thr {
                        offer(&mut sinks[l], q, k, c[l], m[l]);
                    }
                }
            }
        }
        std::mem::swap(&mut prev, &mut cur);
        std::mem::swap(&mut pmeta, &mut cmeta);
    }
}

/// Kept out of line so the hot loop stays branch-free; most cells never get here.
#[cold]
#[inline(never)]
fn offer<S: Sink>(sink: &mut S, q: usize, k: usize, score: i32, meta: u64) {
    sink.cell(q, k, score, meta);
}

/// Offer a cell at or above threshold as its path's best.
#[allow(clippy::too_many_arguments)]
fn record(
    p: &Params,
    stride: u64,
    q: usize,
    k: usize,
    score: i32,
    meta: u64,
    held: &mut Option<Candidate>,
    best: &mut HashMap<u64, Candidate>,
) {
    let (q_in, k_in, gaps, hidden) = unpack(q, meta);
    let cand = Candidate { score, q_out: q, k_out: k, origin: q_in as u64 * stride + k_in as u64, gaps, hidden };
    // The match cutoff constrains which cell is a path's best, rather than filtering afterwards:
    // otherwise a path whose best-scoring extent runs into low-identity flanks would fail the
    // cutoff and take its high-identity core down with it.
    // Rounded to 2 decimals before comparing, as dsRNAscan 0.5.x does.
    if p.min_match_perc > 0.0 && crate::round2(cand.to_hit(p).match_perc()) < p.min_match_perc {
        return;
    }
    match held {
        Some(h) if h.origin == cand.origin => {
            if cand.beats(h) {
                *h = cand
            }
        }
        _ => {
            if let Some(h) = held.replace(cand) {
                flush(h, best);
            }
        }
    }
}

fn flush(c: Candidate, best: &mut HashMap<u64, Candidate>) {
    best.entry(c.origin)
        .and_modify(|b| {
            if c.beats(b) {
                *b = c
            }
        })
        .or_insert(c);
}
