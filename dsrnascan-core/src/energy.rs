//! Turner 2004 duplex energies at a given temperature, rescaled from the 37 C free energies and the
//! enthalpies in `turner2004.rs` the way ViennaRNA does: G(T) = H - (H - G37) * (T + 273.15) / 310.15,
//! truncated towards zero; helix-end terms (terminal mismatches, dangles) are capped at 0 and the
//! asymmetry cap is not rescaled. Base index 4 is N, scored as ViennaRNA does (its own entry, the worst
//! of A/C/G/U at 37 C and in enthalpy).

use crate::turner2004 as t;
use std::sync::{Mutex, OnceLock};

pub const MAXLOOP: usize = 30;

pub struct Energy {
    pub celsius: f64,
    pub init: i32,
    pub term_au: i32,
    pub ninio: i32,
    pub ninio_max: i32,
    pub stack: [[i32; 6]; 6],
    pub tmm: [[[i32; 5]; 5]; 6],
    pub d3: [[i32; 5]; 6],
    pub d5: [[i32; 5]; 6],
    pub bulge: [i32; 31],
    pub int_init: [i32; 31],
    pub mmi: [[[i32; 5]; 5]; 6],
    pub mm23i: [[[i32; 5]; 5]; 6],
    pub mm1ni: [[[i32; 5]; 5]; 6],
    pub int11: [[[[i32; 5]; 5]; 6]; 6],
    pub int21: [[[[[i32; 5]; 5]; 5]; 6]; 6],
    pub int22: [[[[i32; 25]; 25]; 6]; 6],
    /// Lowest energy any loop with n1 x n2 unpaired nt can have, over all sequence contexts. Lets
    /// the DP skip loops that cannot beat the current best without changing its result.
    pub floor: [[i32; MAXLOOP + 1]; MAXLOOP + 1],
}

trait Rescale: Sized {
    fn rescale(&self, dh: &Self, tf: f64) -> Self;
    fn min_entry(&self) -> i32;
}

impl Rescale for i32 {
    fn rescale(&self, dh: &i32, tf: f64) -> i32 {
        (*dh as f64 - (*dh - *self) as f64 * tf) as i32 // `as` truncates towards zero, as C does
    }
    fn min_entry(&self) -> i32 {
        *self
    }
}

impl<T: Rescale, const N: usize> Rescale for [T; N] {
    fn rescale(&self, dh: &Self, tf: f64) -> Self {
        std::array::from_fn(|k| self[k].rescale(&dh[k], tf))
    }
    fn min_entry(&self) -> i32 {
        self.iter().map(|x| x.min_entry()).min().unwrap()
    }
}

impl Energy {
    pub fn at(celsius: f64) -> Energy {
        let tf = (celsius + 273.15) / (37.0 + 273.15);
        let cap0 = |d: [[i32; 5]; 6]| d.map(|r| r.map(|v| v.min(0)));
        let cap0_3 = |d: [[[i32; 5]; 5]; 6]| d.map(|m| m.map(|r| r.map(|v| v.min(0))));
        let mut e = Energy {
            celsius,
            init: t::INIT.rescale(&t::INIT_DH, tf),
            term_au: t::TERM_AU.rescale(&t::TERM_AU_DH, tf),
            ninio: t::NINIO.rescale(&t::NINIO_DH, tf),
            ninio_max: t::NINIO_MAX,
            stack: t::STACK.rescale(&t::STACK_DH, tf),
            tmm: cap0_3(t::TMM.rescale(&t::TMM_DH, tf)),
            d3: cap0(t::D3.rescale(&t::D3_DH, tf)),
            d5: cap0(t::D5.rescale(&t::D5_DH, tf)),
            bulge: t::BULGE.rescale(&t::BULGE_DH, tf),
            int_init: t::INT_INIT.rescale(&t::INT_INIT_DH, tf),
            mmi: t::MMI.rescale(&t::MMI_DH, tf),
            mm23i: t::MM23I.rescale(&t::MM23I_DH, tf),
            mm1ni: t::MM1NI.rescale(&t::MM1NI_DH, tf),
            int11: t::INT11.rescale(&t::INT11_DH, tf),
            int21: t::INT21.rescale(&t::INT21_DH, tf),
            int22: t::INT22.rescale(&t::INT22_DH, tf),
            floor: [[0; MAXLOOP + 1]; MAXLOOP + 1],
        };
        let (stack, au) = (e.stack.min_entry(), e.term_au.min(0));
        let (int11, int21, int22) = (e.int11.min_entry(), e.int21.min_entry(), e.int22.min_entry());
        let (mmi, mm23, mm1n) = (e.mmi.min_entry(), e.mm23i.min_entry(), e.mm1ni.min_entry());
        for n1 in 0..=MAXLOOP {
            for n2 in 0..=MAXLOOP - n1 {
                e.floor[n1][n2] = match (n1, n2) {
                    (0, 0) => stack,
                    (0, 1) | (1, 0) => e.bulge[1] + stack,
                    (0, n) | (n, 0) => e.bulge[n] + 2 * au,
                    (1, 1) => int11,
                    (1, 2) | (2, 1) => int21,
                    (2, 2) => int22,
                    _ => {
                        let mm = if n1.min(n2) == 1 {
                            mm1n
                        } else if (n1.min(n2), n1.max(n2)) == (2, 3) {
                            mm23
                        } else {
                            mmi
                        };
                        e.int_init[n1 + n2] + e.asym(n1, n2) + 2 * mm
                    }
                };
            }
        }
        e
    }

    /// Interior-loop asymmetry penalty.
    #[inline]
    pub fn asym(&self, n1: usize, n2: usize) -> i32 {
        (self.ninio * (n1 as i32 - n2 as i32).abs()).min(self.ninio_max)
    }
}

/// The energies at `celsius`, built once per temperature and kept for the run.
pub fn energy(celsius: f64) -> &'static Energy {
    static E37: OnceLock<Energy> = OnceLock::new();
    if celsius == 37.0 {
        return E37.get_or_init(|| Energy::at(37.0)); // the default: no lock per fold
    }
    static CACHE: OnceLock<Mutex<Vec<&'static Energy>>> = OnceLock::new();
    let mut cache = CACHE.get_or_init(|| Mutex::new(Vec::new())).lock().unwrap();
    if let Some(e) = cache.iter().find(|e| e.celsius == celsius) {
        return e;
    }
    let e: &'static Energy = Box::leak(Box::new(Energy::at(celsius)));
    cache.push(e);
    e
}
