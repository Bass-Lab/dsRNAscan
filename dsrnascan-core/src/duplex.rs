//! RNA duplex minimum free energy: a port of ViennaRNA's `duplexfold` (src/ViennaRNA/duplex.c,
//! ViennaRNA 2.7.0; Lorenz et al. 2011, Algorithms Mol Biol 6:26) with Turner 2004 parameters
//! (`turner2004.rs`). Same recursion, integer arithmetic and tie rules, so energy, structure and
//! trim positions match `RNA.duplexfold`. Strand 1 pairs 5'->3' against strand 2 3'->5'; no
//! intramolecular pairs; loops of up to 30 unpaired nucleotides.

use crate::turner2004 as t;

const MAXLOOP: usize = 30;
const INF: i32 = 10_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Duplex {
    /// dcal/mol (kcal/mol x 100)
    pub energy: i32,
    /// `RNA.duplexfold`'s dot-bracket, including one flanking unpaired nt per side when present
    pub structure: String,
    /// `RNA.duplexfold`'s `.i`: 1-based index of the last strand-1 nt the structure covers
    pub i: usize,
    /// `RNA.duplexfold`'s `.j`: 1-based index of the first strand-2 nt the structure covers
    pub j: usize,
    pub pairs: usize,
    /// longest run of directly stacked pairs (the ML models' `bp_count`)
    pub longest_helix: usize,
}

fn base(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

/// Pair index AU CG GC GU UA UG (strand-1 base, strand-2 base), as in turner2004.rs.
fn pair_type(a: u8, b: u8) -> Option<usize> {
    match (a, b) {
        (0, 3) => Some(0),
        (1, 2) => Some(1),
        (2, 1) => Some(2),
        (2, 3) => Some(3),
        (3, 0) => Some(4),
        (3, 2) => Some(5),
        _ => None,
    }
}

/// The same pair read from the other strand (AU <-> UA, CG <-> GC, GU <-> UG).
const REV: [usize; 6] = [4, 2, 1, 5, 0, 3];

fn is_au_gu(p: usize) -> bool {
    matches!(p, 0 | 3 | 4 | 5)
}

/// ViennaRNA scores an unknown base (N) as the least favourable of A/C/G/U for the term it sits in.
#[inline]
fn worst<const K: usize>(bases: [u8; K], f: impl Fn([usize; K]) -> i32) -> i32 {
    if bases.iter().all(|&b| b < 4) {
        return f(bases.map(|b| b as usize));
    }
    worst_with_n(bases, f)
}

#[cold]
fn worst_with_n<const K: usize>(bases: [u8; K], f: impl Fn([usize; K]) -> i32) -> i32 {
    let ns: Vec<usize> = (0..K).filter(|&k| bases[k] > 3).collect();
    let mut idx = bases.map(|b| (b as usize).min(3));
    let mut worst = i32::MIN;
    for code in 0..4usize.pow(ns.len() as u32) {
        let mut c = code;
        for &k in &ns {
            idx[k] = c % 4;
            c /= 4;
        }
        worst = worst.max(f(idx));
    }
    worst
}

/// Helix end: pair (top, bottom) drawn 5'->3' on top; x = nt 3' of top, y = nt 5' of bottom.
fn end_term(top: u8, bot: u8, x: Option<u8>, y: Option<u8>) -> i32 {
    let p = pair_type(top, bot).expect("end of a duplex must be a pair");
    let au = if is_au_gu(p) { t::TERM_AU } else { 0 };
    au + match (x, y) {
        (Some(x), Some(y)) => worst([x, y], |[x, y]| t::TMM[p][x][y]),
        (Some(x), None) => worst([x], |[x]| t::D3[p][x]),
        (None, Some(y)) => worst([y], |[y]| t::D5[p][y]),
        (None, None) => 0,
    }
}

/// End at strand 1's 5' side (the duplex's first pair), drawn with strand 2 on top.
fn first_end(s1: &[u8], s2: &[u8], i: usize, j: usize) -> i32 {
    end_term(s2[j], s1[i], s2.get(j + 1).copied(), i.checked_sub(1).map(|k| s1[k]))
}

/// End at strand 1's 3' side (the duplex's last pair).
fn last_end(s1: &[u8], s2: &[u8], i: usize, j: usize) -> i32 {
    end_term(s1[i], s2[j], s1.get(i + 1).copied(), j.checked_sub(1).map(|k| s2[k]))
}

/// Energy of the loop closed by outer pair (k, l) and inner pair (i, j), k < i, l > j.
fn loop_energy(s1: &[u8], s2: &[u8], k: usize, l: usize, i: usize, j: usize) -> i32 {
    let (n1, n2) = (i - k - 1, l - j - 1);
    let outer = pair_type(s1[k], s2[l]).unwrap();
    let inner = pair_type(s1[i], s2[j]).unwrap();
    if n1 == 0 && n2 == 0 {
        return t::STACK[outer][inner];
    }
    if n1 == 0 || n2 == 0 {
        let n = n1 + n2;
        if n == 1 {
            return t::BULGE[1] + t::STACK[outer][inner]; // the helix stacks across a 1-nt bulge
        }
        let au = |p| if is_au_gu(p) { t::TERM_AU } else { 0 };
        return t::BULGE[n] + au(outer) + au(inner);
    }
    let (x1, y1) = (s1[k + 1], s2[l - 1]); // loop nt next to the outer pair (strand 1, strand 2)
    let (x2, y2) = (s2[j + 1], s1[i - 1]); // next to the inner pair, read from the far side
    match (n1, n2) {
        (1, 1) => worst([x1, y1], |[x, y]| t::INT11[outer][inner][x][y]),
        (1, 2) => worst([s2[l - 2], x1, y1], |[f, x, y]| t::INT21[outer][inner][f][x][y]),
        (2, 1) => worst([s1[i - 2], x2, y2], |[f, x, y]| t::INT21[REV[inner]][REV[outer]][f][x][y]),
        (2, 2) => worst([x1, y1, y2, x2], |[a, b, c, d]| t::INT22[outer][inner][4 * a + b][4 * c + d]),
        _ => {
            let asym = (t::ASYM_PER_NT * (n1 as i32 - n2 as i32).abs()).min(t::ASYM_CAP);
            let mm = if n1.min(n2) == 1 {
                2 * t::MM_1XN
            } else if (n1.min(n2), n1.max(n2)) == (2, 3) {
                worst([x1, y1, x2, y2], |[a, b, c, d]| t::MM23[outer][a][b] + t::MM23[REV[inner]][c][d])
            } else {
                worst([x1, y1, x2, y2], |[a, b, c, d]| t::MM_GENERIC[a][b] + t::MM_GENERIC[c][d])
            };
            let au = |p| if is_au_gu(p) { t::INT_AU } else { 0 };
            t::INT_INIT[n1 + n2] + asym + mm + au(outer) + au(inner)
        }
    }
}

/// DP table restricted to an anti-diagonal band: cell (i, j) exists iff lo <= i + j <= hi.
/// A helix runs along an anti-diagonal; bulges and asymmetric loops drift off it.
struct Band {
    lo: usize,
    hi: usize,
    n2: usize,
    row_start: Vec<usize>, // row i's cells begin at c[row_start[i]]: storage is exactly the band
    c: Vec<i32>,
}

impl Band {
    fn new(lo: usize, hi: usize, n1: usize, n2: usize) -> Band {
        let mut b = Band { lo, hi, n2, row_start: Vec::with_capacity(n1), c: Vec::new() };
        let mut total = 0;
        for i in 0..n1 {
            b.row_start.push(total);
            total += b.cols(i).len();
        }
        b.c = vec![INF; total];
        b
    }
    /// Columns of row i inside the band.
    #[inline]
    fn cols(&self, i: usize) -> std::ops::Range<usize> {
        let a = self.lo.saturating_sub(i);
        let b = (self.hi + 1).saturating_sub(i).min(self.n2);
        a.min(b)..b
    }
    /// Position of (i, j) in the table, if inside the band.
    #[inline]
    fn slot(&self, i: usize, j: usize) -> Option<usize> {
        let r = self.cols(i);
        r.contains(&j).then(|| self.row_start[i] + j - r.start)
    }
    #[inline]
    fn get(&self, i: usize, j: usize) -> i32 {
        self.slot(i, j).map_or(INF, |s| self.c[s])
    }
}

/// Lowest energy any loop with n1 x n2 unpaired nt can have, over all sequence contexts.
/// Lets the DP skip loops that cannot beat the current best without changing its result.
fn loop_floor() -> &'static [[i32; MAXLOOP + 1]; MAXLOOP + 1] {
    static FLOOR: std::sync::OnceLock<[[i32; MAXLOOP + 1]; MAXLOOP + 1]> = std::sync::OnceLock::new();
    FLOOR.get_or_init(|| {
        let min_of = |it: &mut dyn Iterator<Item = i32>| it.min().unwrap();
        let stack = min_of(&mut t::STACK.iter().flatten().copied());
        let int11 = min_of(&mut t::INT11.iter().flatten().flatten().flatten().copied());
        let int21 = min_of(&mut t::INT21.iter().flatten().flatten().flatten().flatten().copied());
        let int22 = min_of(&mut t::INT22.iter().flatten().flatten().flatten().copied());
        let mm23 = min_of(&mut t::MM23.iter().flatten().flatten().copied());
        let mmg = min_of(&mut t::MM_GENERIC.iter().flatten().copied());
        let mut f = [[INF; MAXLOOP + 1]; MAXLOOP + 1];
        for n1 in 0..=MAXLOOP {
            for n2 in 0..=MAXLOOP - n1 {
                f[n1][n2] = match (n1, n2) {
                    (0, 0) => stack,
                    (0, n) | (n, 0) if n == 1 => t::BULGE[1] + stack,
                    (0, n) | (n, 0) => t::BULGE[n],
                    (1, 1) => int11,
                    (1, 2) | (2, 1) => int21,
                    (2, 2) => int22,
                    _ => {
                        let asym = (t::ASYM_PER_NT * (n1 as i32 - n2 as i32).abs()).min(t::ASYM_CAP);
                        let mm = if n1.min(n2) == 1 {
                            2 * t::MM_1XN
                        } else if (n1.min(n2), n1.max(n2)) == (2, 3) {
                            2 * mm23
                        } else {
                            2 * mmg
                        };
                        t::INT_INIT[n1 + n2] + asym + mm.min(0)
                    }
                };
            }
        }
        f
    })
}

/// Visit every pair (k, l) that can close a loop of <= MAXLOOP nt onto (i, j), in duplex.c's order
/// (k = i-1 downwards, then l = j+1 upwards), with the energy of reaching (i, j) through it.
/// Stops early when `visit` returns true.
fn for_each_loop(
    s1: &[u8],
    s2: &[u8],
    c: &Band,
    i: usize,
    j: usize,
    mut visit: impl FnMut(usize, usize, i32) -> bool,
) {
    let n2 = s2.len();
    for k in (i.saturating_sub(MAXLOOP + 1)..i).rev() {
        for l in j + 1..n2 {
            if (i - k - 1) + (l - j - 1) > MAXLOOP {
                break;
            }
            let ckl = c.get(k, l);
            if ckl < INF && visit(k, l, ckl + loop_energy(s1, s2, k, l, i, j)) {
                return;
            }
        }
    }
}

/// Minimum free energy duplex of strand 1 (5'->3') with strand 2 (5'->3'), as `RNA.duplexfold`.
pub fn duplexfold(strand1: &[u8], strand2: &[u8]) -> Option<Duplex> {
    duplexfold_band(strand1, strand2, 0, (strand1.len() + strand2.len()).saturating_sub(2))
}

/// The band that contains a search hit's alignment path, widened by `margin`. `gaps` is the hit's
/// bulge count: each bulge moves the path one step off the anti-diagonal, so the path stays within
/// `gaps` of the diagonals through its two ends. Works for arms in either orientation.
pub fn path_band(len1: usize, len2: usize, gaps: usize, margin: usize) -> (usize, usize) {
    let (a, b) = (len1.min(len2).saturating_sub(1), len1.max(len2).saturating_sub(1));
    (a.saturating_sub(gaps + margin), b + gaps + margin)
}

/// `duplexfold` over the band `lo <= i + j <= hi` only (0-based i on strand 1, j on strand 2).
/// Identical to the full fold whenever the full fold's MFE duplex lies inside the band.
pub fn duplexfold_band(strand1: &[u8], strand2: &[u8], lo: usize, hi: usize) -> Option<Duplex> {
    let s1: Vec<u8> = strand1.iter().map(|&c| base(c)).collect();
    let s2: Vec<u8> = strand2.iter().map(|&c| base(c)).collect();
    let (n1, n2) = (s1.len(), s2.len());
    if n1 == 0 || n2 == 0 || lo > hi {
        return None;
    }
    // c[i][j]: best duplex whose last pair (towards strand 1's 3' end) is (i, j)
    let mut c = Band::new(lo, hi.min(n1 + n2 - 2), n1, n2);
    let (mut emin, mut imin, mut jmin) = (INF, 0, 0);
    for i in 0..n1 {
        for j in c.cols(i).rev() {
            if pair_type(s1[i], s2[j]).is_none() {
                continue;
            }
            let mut v = t::INIT + first_end(&s1, &s2, i, j);
            // Same minimum as visiting every loop, but a loop whose floor cannot beat the current
            // best is skipped before its energy is computed. Traceback below is unchanged.
            let floor = loop_floor();
            for k in (i.saturating_sub(MAXLOOP + 1)..i).rev() {
                let n1l = i - k - 1;
                // row k's band columns, resolved once: l runs over j+1.. within the loop-size limit
                let cols = c.cols(k);
                let (lo, hi) = ((j + 1).max(cols.start), (j + 2 + MAXLOOP - n1l).min(cols.end));
                if lo >= hi {
                    continue;
                }
                let rs = &c.c[c.row_start[k] + lo - cols.start..c.row_start[k] + hi - cols.start];
                let fs = &floor[n1l][lo - j - 1..hi - j - 1];
                // Vectorised quick reject: if no loop in this row can beat v even at its floor, skip
                // the row. Most rows are skipped here; the rest take the exact path below.
                if rs.iter().zip(fs).map(|(a, b)| a + b).min().unwrap_or(INF) >= v {
                    continue;
                }
                for (x, (&ckl, &f)) in rs.iter().zip(fs).enumerate() {
                    if ckl < INF && ckl + f < v {
                        v = v.min(ckl + loop_energy(&s1, &s2, k, lo + x, i, j));
                    }
                }
            }
            let s = c.slot(i, j).unwrap();
            c.c[s] = v;
            let e = v + last_end(&s1, &s2, i, j);
            if e < emin {
                emin = e;
                imin = i;
                jmin = j;
            }
        }
    }
    if emin >= INF {
        return None;
    }

    // Traceback from (imin, jmin) towards strand 1's 5' end; first matching loop wins, as in duplex.c.
    let (mut st1, mut st2) = (vec![b'.'; n1], vec![b'.'; n2]);
    let (mut i, mut j) = (imin, jmin);
    let (mut pairs, mut run, mut longest) = (0, 0, 0);
    loop {
        st1[i] = b'(';
        st2[j] = b')';
        pairs += 1;
        run += 1;
        longest = longest.max(run);
        let target = c.get(i, j);
        let mut next = None;
        for_each_loop(&s1, &s2, &c, i, j, |k, l, e| {
            if e == target {
                next = Some((k, l));
            }
            next.is_some()
        });
        match next {
            Some((k, l)) => {
                if !(k + 1 == i && l == j + 1) {
                    run = 0; // a bulge or interior loop ends the helix
                }
                i = k;
                j = l;
            }
            None => {
                debug_assert_eq!(target - first_end(&s1, &s2, i, j), t::INIT);
                break;
            }
        }
    }
    // one flanking nucleotide per side, as duplexfold reports
    let start1 = i.saturating_sub(1);
    let end1 = (imin + 1).min(n1 - 1);
    let start2 = jmin.saturating_sub(1);
    let end2 = (j + 1).min(n2 - 1);
    let mut structure = String::with_capacity(end1 - start1 + end2 - start2 + 3);
    structure.push_str(std::str::from_utf8(&st1[start1..=end1]).unwrap());
    structure.push('&');
    structure.push_str(std::str::from_utf8(&st2[start2..=end2]).unwrap());
    Some(Duplex {
        energy: emin,
        structure,
        i: end1 + 1,
        j: start2 + 1,
        pairs,
        longest_helix: longest,
    })
}
