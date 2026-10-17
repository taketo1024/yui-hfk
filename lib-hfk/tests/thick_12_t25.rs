//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T2,5 ([SS] (v2), Tables 1 and 2).
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! With `config()` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided. Those of this type that need a larger `max_size` (12n549)
//! are in `thick_12_heavy.rs`.
//!
//! The tests are slow:
//!
//!     cargo test --release --test thick_12_t25 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::{run_all, run_tau};
use common::types::{K, T2_5};

fn config() -> G0CalcConfig {
    G0CalcConfig { try_fast_ng: true, max_size: 2_000_000, ..Default::default() }
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn tau_all() {
    let knots: Vec<_> = KNOTS.iter().collect();
    run_tau("tau_all", &knots, config());
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_all() {
    let knots: Vec<_> = KNOTS.iter().collect();
    run_all("g0_all", &knots, |k| k.name, |k| k.check(config()));
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    K { name: "12n16", os: &[9,10,6,7,0,11,8,5,3,4,1,2], xs: &[11,8,9,4,6,3,2,1,10,0,7,5], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n102", os: &[6,0,8,3,2,1,11,7,4,5,10,9], xs: &[2,5,1,0,10,6,4,3,8,9,7,11], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n120", os: &[1,0,2,6,4,3,5,7,8,9,10], xs: &[5,3,9,1,0,10,2,4,6,7,8], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n149", os: &[8,9,5,7,6,0,10,2,1,3,4], xs: &[10,7,8,3,2,5,4,9,6,0,1], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n175", os: &[11,10,9,8,7,6,0,1,4,2,5,3], xs: &[9,7,5,4,11,2,3,8,0,10,1,6], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n179", os: &[11,10,9,7,6,4,0,8,2,3,5,1], xs: &[9,8,5,4,2,11,3,1,7,10,0,6], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n199", os: &[1,0,2,5,4,3,6,7,8,9,10], xs: &[4,3,9,1,0,10,2,5,6,7,8], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n293", os: &[7,5,4,1,2,9,0,6,8,3], xs: &[9,8,6,5,7,3,4,1,2,0], mirror: false, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n321", os: &[7,5,2,3,1,9,0,6,8,4], xs: &[9,8,6,7,4,2,5,1,3,0], mirror: false, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n332", os: &[8,9,5,7,6,0,10,1,3,2,4], xs: &[10,7,8,3,2,5,4,6,0,9,1], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n349", os: &[8,9,4,7,10,6,0,1,2,3,5], xs: &[10,6,8,3,5,2,4,9,7,0,1], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n398", os: &[5,0,9,1,2,3,10,6,4,7,8], xs: &[1,4,3,6,7,0,5,2,8,9,10], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n404", os: &[8,9,6,7,5,10,0,1,2,4,3], xs: &[10,7,8,2,1,4,3,6,5,0,9], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n432", os: &[8,9,5,7,6,0,10,1,4,2,3], xs: &[10,7,8,2,1,4,3,5,9,6,0], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n433", os: &[9,1,0,3,2,4,5,6,8,7,10], xs: &[2,6,4,7,10,1,0,3,5,9,8], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n436", os: &[10,6,8,7,1,3,0,2,4,9,5], xs: &[8,9,4,3,5,10,6,7,1,2,0], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n526", os: &[9,0,10,5,3,4,2,1,7,6,8], xs: &[3,7,2,1,0,8,5,6,4,9,10], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n617", os: &[6,9,0,2,1,4,3,5,8,7,10], xs: &[1,3,5,6,7,0,10,2,4,9,8], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n642", os: &[7,9,8,0,2,10,4,1,5,3,6], xs: &[10,4,1,3,7,5,8,6,0,9,2], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n697", os: &[1,0,10,8,5,2,3,6,4,9,7,11], xs: &[6,4,3,1,0,7,11,2,8,5,10,9], mirror: true, g: 4, g4: 2, tau: -2, g0_type: T2_5 },
    K { name: "12n739", os: &[3,5,7,10,0,1,2,6,4,9,8], xs: &[9,2,1,3,4,5,8,0,7,6,10], mirror: true, g: 4, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n801", os: &[1,7,0,5,8,3,6,2,4,9], xs: &[5,3,4,2,1,9,0,7,8,6], mirror: false, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
];
