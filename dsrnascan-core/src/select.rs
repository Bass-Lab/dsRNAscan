//! Selection. Each hit is a box in the dot plot (5' arm x 3' arm), kept strongest-first.
//! Two rules decide whether a stronger kept hit `y` suppresses a hit `x`:
//!  - box rule: at least `both_frac` of BOTH of x's arms lies inside the matching arms of y
//!    (1.0 = x is nested in y). Always applies.
//!  - greed: one of x's arms overlaps the same arm of y (5' vs 5', 3' vs 3') by at least
//!    `arm_frac` of x's arm, and x scores below `greed` x y's score. `greed = 1.0` (default) is
//!    einverted-like, one partner per arm; `0.0` turns it off, so arm-sharing hits survive
//!    however weak they are.

use crate::{par_map, Dsrna, Hit};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug)]
pub struct Rule {
    pub both_frac: f64,
    pub arm_frac: f64,
    pub greed: f64,
    /// keep a weaker partner whose stronger rival's other arm lies outside its span (0.5.x windows)
    pub rescue: bool,
}

fn overlap(a0: usize, a1: usize, b0: usize, b1: usize) -> usize {
    (a1.min(b1) + 1).saturating_sub(a0.max(b0))
}

fn suppresses(y: &Hit, x: &Hit, r: &Rule) -> bool {
    let (l5, l3) = ((x.i_end - x.i_start + 1) as f64, (x.j_end - x.j_start + 1) as f64);
    let o55 = overlap(x.i_start, x.i_end, y.i_start, y.i_end) as f64;
    let o33 = overlap(x.j_start, x.j_end, y.j_start, y.j_end) as f64;
    if o55 > 0.0 && o33 > 0.0 && o55 >= r.both_frac * l5 && o33 >= r.both_frac * l3 {
        return true;
    }
    if (x.score as f64) >= r.greed * y.score as f64 {
        return false;
    }
    // Same-role arms only (5' vs 5', 3' vs 3'), as einverted tracks them. An arm used as the 5'
    // arm of one hit and the 3' arm of another (alternating Alus: B pairs with A upstream and with
    // C downstream) is not suppressed: counting it lost 261/725 of 0.5.x's dsRNAs on chr19.
    let shares5 = o55 > 0.0 && o55 >= r.arm_frac * l5;
    let shares3 = o33 > 0.0 && o33 >= r.arm_frac * l3;
    if !r.rescue {
        return shares5 || shares3;
    }
    // Rescue, as 0.5.x's overlapping windows did: a weaker partner of a shared arm survives when
    // the stronger hit's other arm lies outside the weaker hit's span, because some window then
    // held the weaker pair without the stronger one. 305 of the 372 0.5.x dsRNAs that plain
    // greed loses on chr19 are this case.
    let inside = |a0: usize, a1: usize| x.i_start <= a0 && a1 <= x.j_end;
    (shares5 && inside(y.j_start, y.j_end)) || (shares3 && inside(y.i_start, y.i_end))
}

/// Kept hits indexed by the start of each arm. Arms are shorter than max_span, so any kept arm
/// overlapping [s, e] starts within [s - max_span, e].
struct Kept {
    width: usize,
    grid: HashMap<usize, Vec<usize>>,
    hits: Vec<Hit>,
}

impl Kept {
    fn near(&self, s: usize, e: usize) -> impl Iterator<Item = &Hit> + '_ {
        (s.saturating_sub(self.width) / self.width..=e / self.width)
            .filter_map(|b| self.grid.get(&b))
            .flatten()
            .map(|&k| &self.hits[k])
    }
    fn covered(&self, h: &Hit, r: &Rule) -> bool {
        self.near(h.i_start, h.i_end).chain(self.near(h.j_start, h.j_end)).any(|y| suppresses(y, h, r))
    }
    fn add(&mut self, h: Hit) {
        let k = self.hits.len();
        self.grid.entry(h.i_start / self.width).or_default().push(k);
        if h.j_start / self.width != h.i_start / self.width {
            self.grid.entry(h.j_start / self.width).or_default().push(k);
        }
        self.hits.push(h);
    }
}

/// Strongest first: higher score, then position (deterministic).
fn strength(a: &Hit, b: &Hit) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.i_start.cmp(&b.i_start))
        .then(a.j_end.cmp(&b.j_end))
        .then(a.i_end.cmp(&b.i_end))
        .then(a.j_start.cmp(&b.j_start))
}

fn order(hits: &mut [Hit]) {
    hits.sort_by(strength);
}

/// dsRNAscan 0.5.x's nested elimination (`eliminate_nested_dsrnas`), exactly: a dsRNA is dropped
/// when another dsRNA's arms contain both of its arms (i in i, j in j), whatever the scores; of
/// coordinate-identical dsRNAs the first is kept. Input is one strand of one sequence.
pub fn nested_05x(mut ds: Vec<Dsrna>) -> Vec<Dsrna> {
    if ds.len() <= 1 {
        return ds;
    }
    ds.sort_by_key(|d| d.hit.i_start);
    let starts: Vec<usize> = ds.iter().map(|d| d.hit.i_start).collect();
    let max_arm = ds.iter().map(|d| d.hit.i_end - d.hit.i_start).max().unwrap();
    let nested: Vec<bool> = (0..ds.len())
        .map(|k| {
            let h = &ds[k].hit;
            // a container starts within max_arm before h and no later than h
            let lo = starts.partition_point(|&s| s < h.i_start.saturating_sub(max_arm));
            let hi = starts.partition_point(|&s| s <= h.i_start);
            (lo..hi).any(|c| {
                let y = &ds[c].hit;
                let contains = y.i_end >= h.i_end && y.j_start <= h.j_start && y.j_end >= h.j_end;
                let identical = (y.i_start, y.i_end, y.j_start, y.j_end) == (h.i_start, h.i_end, h.j_start, h.j_end);
                c != k && contains && !(identical && c > k)
            })
        })
        .collect();
    ds.into_iter().zip(nested).filter(|(_, x)| !x).map(|(d, _)| d).collect()
}

/// Selection without filters.
pub fn select(hits: Vec<Hit>, rule: Rule, max_span: usize) -> Vec<Hit> {
    select_lazy(hits, Some(rule), max_span, |_| Some(())).into_iter().map(|(h, _)| h).collect()
}

/// The box rule with lazy evaluation. Candidates are visited strongest-first. One already covered
/// by a kept hit is skipped without calling `eval`; otherwise `eval` decides, and `None` (fails
/// the filters) is never kept and never suppresses anything. This equals filtering every candidate
/// first and then selecting, but evaluates only candidates that could be kept.
///
/// Parallelism by optimistic replay: run the selection assuming every unevaluated candidate
/// passes, evaluate all of those it would keep in one parallel pass, and replay with the results.
/// Only candidates uncovered by a failed evaluation need another round. The last replay decides
/// every candidate from real results, so the output equals the sequential definition and does not
/// depend on the number of threads. `rule = None` keeps every candidate that passes `eval`.
pub fn select_lazy<T: Send>(
    mut hits: Vec<Hit>,
    rule: Option<Rule>,
    max_span: usize,
    eval: impl Fn(&Hit) -> Option<T> + Sync + Send,
) -> Vec<(Hit, T)> {
    order(&mut hits);
    let Some(rule) = rule else {
        let vals = par_map(&hits, &eval);
        return hits.into_iter().zip(vals).filter_map(|(h, v)| v.map(|v| (h, v))).collect();
    };
    let mut known: HashMap<usize, Option<T>> = HashMap::new();
    loop {
        let mut kept = Kept { width: max_span.max(1), grid: HashMap::new(), hits: Vec::new() };
        let (mut keep, mut needed) = (Vec::new(), Vec::new());
        for (k, h) in hits.iter().enumerate() {
            if kept.covered(h, &rule) {
                continue;
            }
            match known.get(&k) {
                Some(None) => continue, // failed the filters: never kept, never suppresses
                Some(Some(_)) => keep.push(k),
                None => needed.push(k), // not evaluated yet: assume it passes
            }
            kept.add(h.clone());
        }
        if needed.is_empty() {
            return keep
                .into_iter()
                .map(|k| (hits[k].clone(), known.remove(&k).flatten().unwrap()))
                .collect();
        }
        let vals = par_map(&needed, |&k| eval(&hits[k]));
        known.extend(needed.into_iter().zip(vals));
    }
}
