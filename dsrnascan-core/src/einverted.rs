//! dsRNAscan 0.5.x's windows, reproduced from one sweep.
//!
//! 0.5.x runs einverted on overlapping windows. Inside a window, einverted's DP values equal the
//! sweep's: a cell depends only on cells with a larger 5' position and a smaller 3' position, and
//! those lie in the same window. So the windows differ only in einverted's *reporting*, which
//! reads one number per row: the row's best cell among partners inside the window. The sweep
//! records each row's running maximum (`scan::RowSteps`); here einverted's report loop is replayed
//! per window on those rows, ported from EMBOSS einverted (EMBOSS 6.6.0 einverted.c with
//! dsRNAscan's G-U patch; Rice et al. 2000, Trends Genet 16:276): its main loop and
//! einverted_report's window-edge check.

use crate::scan::{self, columns, unpack, Candidate, RowSteps};
use crate::{blocks, par_map, Hit, Params};

/// The hits einverted reports for the window of `len` rows starting at 0-based position `ws`.
pub fn window_report(steps: &RowSteps, ws: usize, len: usize, p: &Params, s: &mut Scratch) -> Vec<Candidate> {
    let m = p.max_span; // einverted's maxsave (-maxrepeat)
    let stride = columns(p) as u64 + 1;
    // einverted clears back/localMax per record; cell is only read where back points
    s.back.fill(0);
    s.local_max.fill(0);
    let (back, local_max, cell) = (&mut s.back, &mut s.local_max, &mut s.cell);
    let mut last_reported: i64 = -1;
    let mut out = Vec::new();

    for i in 0..len + m {
        let irel = i % m;
        if back[irel] != 0 {
            // a 5' position has left the span: report the best registered row up to its row
            let jback = back[irel];
            let lo = i as i64 - m as i64;
            let mut jmax = jback;
            let mut ok = false;
            while !ok && back[irel] != 0 {
                let mut imax = 0;
                jmax = jback;
                let mut j = jback as i64;
                while j > lo {
                    let ts = local_max[j as usize % m];
                    if ts != 0 && ts > imax {
                        jmax = j as usize;
                        imax = ts;
                    }
                    j -= 1;
                }
                if imax == 0 {
                    ok = true;
                    continue;
                }
                // einverted_report: the traceback fails if the path reaches below row i - maxsave
                let (k, meta) = cell[jmax % m];
                let (q_in, k_in, gaps, hidden) = unpack(jmax, meta);
                ok = q_in as i64 >= lo;
                if ok {
                    out.push(Candidate {
                        score: imax,
                        q_out: ws + jmax,
                        k_out: k,
                        origin: (ws + q_in) as u64 * stride + k_in as u64,
                        gaps,
                        hidden,
                    });
                } else {
                    local_max[jmax % m] = 0;
                }
            }
            last_reported = jmax as i64;
            let mut j = jmax as i64;
            while j >= lo {
                let s = j.rem_euclid(m as i64) as usize;
                local_max[s] = 0;
                back[s] = 0;
                j -= 1;
            }
        }
        if i >= len {
            continue; // report only
        }
        // row i: its best cell among partners inside the window (columns k <= i - 1)
        // capped at this run's columns (smaller runs can share a longer-span sweep)
        let best = if i >= 2 { steps.best(ws + i, (i - 1).min(columns(p).saturating_sub(1))) } else { None };
        match best {
            Some(st) => {
                let jstart = i - st.k as usize - 1; // the cell's 5' partner
                let slot = jstart % m;
                let jb = back[slot];
                if jstart as i64 > last_reported && (jb == 0 || local_max[jb % m] < st.score) {
                    back[slot] = i;
                    local_max[irel] = st.score;
                    cell[irel] = (st.k as usize, st.meta);
                }
            }
            None => local_max[irel] = 0,
        }
    }
    out
}

/// Every hit any 0.5.x window would report (windows of `p.window` every `step`, full-length
/// windows only as in 0.5.x, one window if the sequence is shorter), deduplicated.
///
/// Memory is bounded by chunk size, not sequence length. Window starts are taken in chunks; a
/// chunk sweeps only the rows its windows cover (plus the halo that makes them exact), replays
/// its windows one after another, and drops its rows. Chunks run in parallel, one per thread.
pub fn window_hits(enc: &[u8], p: &Params, step: usize) -> Vec<Hit> {
    window_hits_lanes(enc, p, step, false).pop().unwrap()
}

/// `window_hits` for both strands from one fused sweep of the plus-strand sequence: [plus, minus].
/// Identical to two separate calls (the minus one on the complement), at much less sweep cost.
pub fn window_hits_both(enc: &[u8], p: &Params, step: usize) -> [Vec<Hit>; 2] {
    let mut v = window_hits_lanes(enc, p, step, true);
    let minus = v.pop().unwrap();
    [v.pop().unwrap(), minus]
}

fn window_hits_lanes(enc: &[u8], p: &Params, step: usize, both: bool) -> Vec<Vec<Hit>> {
    let n = enc.len();
    let win = if p.window > 0 { p.window } else { p.max_span };
    let (w, halo) = (win.min(n), p.max_span + 2);
    let step = step.max(1);
    let starts: Vec<usize> = if n < win { vec![0] } else { (0..=n - win).step_by(step).collect() };
    // at most 250 kb of window starts per chunk (rows held: chunk + one window)
    let span = blocks(n, p).first().map_or(n, |&(s, e)| e - s).min(250_000);
    let chunks: Vec<&[usize]> = starts.chunks((span / step).max(1)).collect();
    let per_chunk: Vec<Vec<Vec<Candidate>>> = par_map(&chunks, |starts| {
        let (first, last) = (starts[0], *starts.last().unwrap());
        let own = first..(last + w).min(n);
        // a row's DP values are exact once the sweep starts max_span rows earlier
        let rows = first.saturating_sub(halo)..own.end;
        let lanes: Vec<RowSteps> = if both {
            scan::row_steps_both(enc, p, rows, own, step, win).into()
        } else {
            vec![scan::row_steps(enc, p, rows, own, step, win, &[])]
        };
        let mut scratch = Scratch::new(p.max_span);
        lanes
            .iter()
            .map(|steps| starts.iter().flat_map(|&ws| window_report(steps, ws, w, p, &mut scratch)).collect())
            .collect()
    });
    (0..if both { 2 } else { 1 })
        .map(|lane| {
            let mut hits: Vec<Hit> = per_chunk.iter().flat_map(|c| c[lane].iter().map(|c| c.to_hit(p))).collect();
            hits.sort_by_key(|h| (h.i_start, h.j_end, h.i_end, h.j_start, h.score));
            hits.dedup();
            hits
        })
        .collect()
}

/// einverted's per-record tables, reused across the windows of a chunk.
pub struct Scratch {
    back: Vec<usize>,          // by 5' position mod m: the row registered for it (0 = none)
    local_max: Vec<i32>,       // by row mod m: that row's registered score
    cell: Vec<(usize, u64)>,   // by row mod m: the registered cell (column, metadata)
}

impl Scratch {
    pub fn new(m: usize) -> Self {
        Scratch { back: vec![0; m], local_max: vec![0; m], cell: vec![(0, 0); m] }
    }
}

/// Several window sizes from ONE sweep of one strand, merged (`--windows`, `--v1`). A cell's value
/// and path depend only on cells with the same or a smaller span, so a run with a smaller span
/// limit computes exactly the shared sweep's cells below its limit. The sweep uses the largest
/// span; each run replays its own windows, looking rows up no further than its own column cap.
/// `runs` = (window size, that run's Params). Same hits as running each size on its own.
pub fn window_hits_multi(enc: &[u8], runs: &[(usize, Params)], step: usize) -> Vec<Hit> {
    let n = enc.len();
    let step = step.max(1);
    let shared = runs.iter().max_by_key(|(_, q)| q.max_span).map(|(_, q)| q.clone()).unwrap();
    let wmax = runs.iter().map(|&(w, _)| w).max().unwrap();
    let halo = shared.max_span + 2;
    let caps: Vec<usize> = runs.iter().map(|(_, q)| columns(q).saturating_sub(1)).collect();
    let all: Vec<Vec<usize>> =
        runs.iter().map(|&(w, _)| if n < w { vec![0] } else { (0..=n - w).step_by(step).collect() }).collect();
    // chunk the position axis; each chunk replays every run's windows that start inside it
    let span = blocks(n, &shared).first().map_or(n, |&(s, e)| e - s).min(250_000).max(step);
    let bounds: Vec<usize> = (0..n.max(1)).step_by(span).collect();
    let per_chunk: Vec<Vec<Hit>> = par_map(&bounds, |&c0| {
        let mine: Vec<(usize, &[usize])> = all
            .iter()
            .enumerate()
            .map(|(r, st)| (r, &st[st.partition_point(|&s| s < c0)..st.partition_point(|&s| s < c0 + span)]))
            .filter(|(_, st)| !st.is_empty())
            .collect();
        let Some(end) = mine.iter().map(|&(r, st)| (st[st.len() - 1] + runs[r].0).min(n)).max() else {
            return Vec::new();
        };
        let first = mine.iter().map(|&(_, st)| st[0]).min().unwrap();
        let steps = scan::row_steps(enc, &shared, first.saturating_sub(halo)..end, first..end, step, wmax, &caps);
        let mut out = Vec::new();
        for (r, st) in mine {
            let (w, q) = (runs[r].0, &runs[r].1);
            let mut scratch = Scratch::new(q.max_span);
            for &ws in st {
                out.extend(window_report(&steps, ws, w.min(n), q, &mut scratch).into_iter().map(|c| c.to_hit(q)));
            }
        }
        out
    });
    let mut hits: Vec<Hit> = per_chunk.into_iter().flatten().collect();
    hits.sort_by_key(|h| (h.i_start, h.j_end, h.i_end, h.j_start, h.score));
    hits.dedup();
    hits
}
