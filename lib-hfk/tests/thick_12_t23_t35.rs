//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T2,3#T3,5 ([SS] (v2), Table 1). This type is not in Table 2
//! of [SS]. It is computed here for all seven knots.
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! With the config in `check` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided.
//!
//! The test is slow:
//!
//!     cargo test --release --test thick_12_t23_t35 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T2_3_T3_5};

fn check(k: &K) -> usize {
    k.check(G0CalcConfig { try_fast_ng: true, max_size: 2_000_000, ..Default::default() })
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0() {
    let knots: Vec<_> = KNOTS.iter().collect();
    run_all("g0", &knots, |k| k.name, check);
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "12n242", os: &[5,4,2,1,0,8,7,6,3], xs: &[8,7,6,5,4,3,2,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n472", os: &[8,7,5,3,2,1,10,9,0,6,4], xs: &[10,9,8,7,6,5,4,2,3,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n574", os: &[8,7,6,4,2,1,0,10,9,5,3], xs: &[10,9,8,7,6,5,4,3,1,2,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n679", os: &[9,8,7,5,3,2,11,1,10,0,6,4], xs: &[11,10,9,8,7,6,5,4,2,3,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
    K { name: "12n688", os: &[9,8,6,4,3,2,11,1,7,0,10,5], xs: &[11,10,9,8,7,6,5,3,2,4,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
    K { name: "12n725", os: &[7,5,3,2,1,9,8,0,6,4], xs: &[9,8,7,6,5,4,3,2,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n888", os: &[9,7,5,4,3,11,2,1,8,0,10,6], xs: &[11,10,9,8,7,4,6,5,2,3,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
];
