//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T2,3 ([SS] (v2), Tables 1 and 2).
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! With `config()` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided. Those of this type that need a larger `max_size` (12n45,
//! 12n382) are in `thick_12_heavy.rs`.
//!
//! The tests are slow:
//!
//!     cargo test --release --test thick_12_t23 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::{run_all, run_tau};
use common::types::{K, T2_3};

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
    K { name: "12n54", os: &[1,2,0,4,6,3,5,7,8,9,10], xs: &[5,9,3,1,2,10,0,4,6,7,8], mirror: true, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n63", os: &[11,8,9,7,0,4,1,3,2,10,5,6], xs: &[9,10,5,3,6,11,7,0,8,4,1,2], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n82", os: &[11,8,9,4,1,3,0,2,6,10,7,5], xs: &[9,10,6,7,5,8,4,11,1,2,3,0], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n115", os: &[1,0,6,10,4,11,2,3,5,7,8,9], xs: &[4,3,1,2,0,5,6,7,8,9,10,11], mirror: true, g: 4, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n119", os: &[10,9,8,6,5,0,3,7,2,1,4], xs: &[8,7,5,3,2,4,10,1,0,6,9], mirror: false, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n121", os: &[2,4,0,8,1,5,9,3,6,7], xs: &[5,1,3,2,6,0,4,7,8,9], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n143", os: &[8,7,4,0,10,2,1,9,3,6,5], xs: &[10,9,8,6,3,7,5,4,0,1,2], mirror: false, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n184", os: &[6,5,0,2,7,8,4,1,3,9,10,11], xs: &[2,1,4,10,0,3,11,6,7,5,8,9], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n192", os: &[1,0,7,10,5,11,2,3,4,6,8,9], xs: &[5,4,1,2,0,3,6,7,8,9,10,11], mirror: true, g: 4, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n195", os: &[9,10,6,8,3,11,2,0,1,7,4,5], xs: &[11,8,9,4,7,5,6,3,10,2,0,1], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n200", os: &[8,9,6,7,2,4,0,10,1,5,3], xs: &[10,7,8,4,5,1,3,2,6,0,9], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n218", os: &[6,7,1,0,2,5,4,3,8,9,10], xs: &[9,2,4,3,10,1,0,6,5,7,8], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n225", os: &[6,5,0,10,7,2,4,1,3,9,8,11], xs: &[2,1,4,3,0,8,11,6,7,5,10,9], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n230", os: &[9,10,6,7,0,11,8,5,3,4,1,2], xs: &[11,8,9,4,6,3,1,2,10,0,5,7], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n231", os: &[9,10,6,7,0,11,8,1,3,4,5,2], xs: &[11,8,9,4,5,3,2,6,10,0,1,7], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n232", os: &[9,10,6,7,11,4,8,3,5,0,1,2], xs: &[11,8,9,3,5,1,2,6,10,4,7,0], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n260", os: &[8,9,4,5,0,7,6,10,1,3,2], xs: &[10,6,8,1,4,3,2,5,9,0,7], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n267", os: &[11,10,9,4,6,1,0,7,3,8,5,2], xs: &[9,7,5,8,3,4,2,1,11,0,10,6], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n273", os: &[0,9,4,3,2,1,5,7,8,6,10], xs: &[5,2,1,0,7,6,10,3,4,9,8], mirror: false, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n282", os: &[8,9,5,6,2,0,10,3,7,1,4], xs: &[10,7,8,3,5,4,1,0,2,6,9], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n310", os: &[8,9,5,4,1,10,7,0,6,2,3], xs: &[10,7,8,6,5,2,3,4,1,9,0], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n336", os: &[2,1,0,4,5,8,3,6,7,9,10,11], xs: &[7,6,3,10,2,1,11,0,4,5,8,9], mirror: true, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n346", os: &[9,8,6,4,7,11,1,2,0,3,10,5], xs: &[11,10,9,8,3,5,6,7,4,1,2,0], mirror: true, g: 4, g4: 1, tau: -1, g0_type: T2_3 },
    K { name: "12n355", os: &[8,9,3,2,0,10,4,6,5,7,1], xs: &[10,6,8,7,5,1,0,2,9,3,4], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n370", os: &[7,5,3,0,1,2,9,4,6,8], xs: &[9,8,6,5,4,7,3,0,1,2], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n371", os: &[0,6,3,2,5,4,1,7,8,9], xs: &[4,2,1,8,0,9,6,3,5,7], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n403", os: &[6,8,0,1,9,3,4,2,5,7], xs: &[3,4,5,7,2,0,1,6,8,9], mirror: true, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n411", os: &[8,9,4,3,6,10,0,7,1,2,5], xs: &[10,6,8,7,1,5,3,2,4,9,0], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n438", os: &[6,5,7,1,9,0,8,4,3,2], xs: &[9,8,4,6,3,5,2,1,0,7], mirror: false, g: 3, g4: 2, tau: 1, g0_type: T2_3 },
    K { name: "12n439", os: &[7,5,3,6,9,0,8,1,2,4], xs: &[9,8,7,2,4,5,3,6,0,1], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n446", os: &[6,0,9,5,1,3,4,2,8,7,10], xs: &[1,4,3,2,7,0,10,6,5,9,8], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n457", os: &[5,0,9,2,1,3,10,6,4,7,8], xs: &[1,4,3,7,6,0,5,2,8,9,10], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n464", os: &[8,6,4,7,0,10,5,3,2,9,1], xs: &[10,9,8,2,6,1,0,7,4,3,5], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n487", os: &[7,5,3,1,0,2,9,6,4,8], xs: &[9,8,6,5,4,7,3,0,1,2], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n495", os: &[11,7,9,8,3,4,2,0,1,10,5,6], xs: &[9,10,5,4,7,11,6,3,8,2,0,1], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n499", os: &[3,5,4,10,2,0,1,6,7,9,8], xs: &[9,2,1,3,7,5,8,0,4,6,10], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n523", os: &[10,7,8,3,4,0,2,1,5,9,6], xs: &[8,9,5,7,1,3,10,6,0,4,2], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n535", os: &[8,9,4,3,10,6,0,7,1,2,5], xs: &[10,6,8,7,5,1,3,2,4,9,0], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n554", os: &[9,10,5,4,7,11,6,0,8,2,3,1], xs: &[11,7,9,8,3,4,2,5,1,10,0,6], mirror: true, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n558", os: &[11,10,9,7,6,1,8,0,3,4,5,2], xs: &[9,8,5,3,11,4,2,7,10,0,1,6], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n601", os: &[0,8,5,3,2,1,6,7,4,9,10,11], xs: &[6,2,1,0,10,4,11,3,8,5,7,9], mirror: true, g: 3, g4: 1, tau: -1, g0_type: T2_3 },
    K { name: "12n629", os: &[8,9,4,3,5,0,10,6,2,7,1], xs: &[10,6,8,7,2,3,1,0,9,4,5], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n665", os: &[6,0,1,5,2,4,3,7,8,9,10,11], xs: &[2,3,10,0,8,1,6,11,4,5,7,9], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n672", os: &[11,7,9,8,3,4,2,0,1,5,10,6], xs: &[9,10,5,4,7,11,6,3,8,0,2,1], mirror: true, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n673", os: &[1,0,10,6,4,11,2,3,5,7,8,9], xs: &[4,3,2,1,0,5,6,7,8,9,10,11], mirror: false, g: 4, g4: 1, tau: -1, g0_type: T2_3 },
    K { name: "12n808", os: &[8,9,3,7,5,0,10,6,4,2,1], xs: &[10,6,8,4,2,3,1,0,9,7,5], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
    K { name: "12n809", os: &[6,8,3,0,5,9,7,4,2,1], xs: &[9,4,7,6,2,3,1,0,8,5], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
];
