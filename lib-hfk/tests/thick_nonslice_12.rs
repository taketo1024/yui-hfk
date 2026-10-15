//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, with `tau != 0`, against the `G0`-types of Tables 1 and 2 in [SS] (v2).
//! Those of `G0`-type `O` (`tau = 0`) are in `tau0_12.rs`.
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `genus`, `g4` and `|tau|` are also from KnotInfo.
//!
//! `regions` is the `G0`-type of the side with `tau > 0`, as in Table 1. Each knot is
//! tested on the side with the smaller linear systems: where `mirror` is set, the test
//! uses the mirror of KnotInfo's diagram, and `tau` is that of the diagram tested. Where
//! `tau < 0`, the expected `G0` is the dual type ([SS] Cor. A.7).
//!
//! With the config in `check` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided, except:
//!
//! - 12n45, 12n220, 12n228, 12n382, 12n549, 12n693: fully verified with `max_size = 5_000_000`.
//! - 12n696: fully verified with `max_size = 7_000_000`.
//! - 12n59, 12n96, 12n203: not decided yet.
//!
//! All tests are slow:
//!
//!     cargo test --release --test thick_nonslice_12 -- --ignored --nocapture --test-threads=1

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T2_3, T2_3_2_1, T2_3_T3_5, T2_5, T2_7, T3_4, T3_5};

fn check(k: &K) -> usize {
    k.check(G0CalcConfig { try_fast_ng: true, max_size: 2_000_000, ..Default::default() })
}

fn run_grid(label: &str, n: usize) {
    let selected: Vec<_> = KNOTS.iter().filter(|k| k.grid_number() == n).collect();
    run_all(label, &selected, |k| k.name, check);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_grid9() {
    run_grid("g0_grid9", 9);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_grid10() {
    run_grid("g0_grid10", 10);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_grid11() {
    run_grid("g0_grid11", 11);
}

#[test]
#[ignore = "slow: run with --release --ignored"]
fn g0_grid12() {
    run_grid("g0_grid12", 12);
}

#[rustfmt::skip]
static KNOTS: &[K] = &[
    // T2,3
    K { name: "12n45", os: &[10,9,8,2,0,1,5,3,7,4,6], xs: &[8,7,5,6,3,4,2,10,0,9,1], mirror: false, g: 3, g4: 1, tau: 1, g0_type: T2_3 },
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
    K { name: "12n382", os: &[8,7,4,0,10,1,2,9,3,5,6], xs: &[10,9,8,6,3,7,5,4,0,1,2], mirror: false, g: 2, g4: 1, tau: 1, g0_type: T2_3 },
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

    // T2,5
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
    K { name: "12n549", os: &[8,7,5,3,10,1,2,0,9,4,6], xs: &[10,9,8,6,4,5,7,3,1,0,2], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n617", os: &[6,9,0,2,1,4,3,5,8,7,10], xs: &[1,3,5,6,7,0,10,2,4,9,8], mirror: true, g: 3, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n642", os: &[7,9,8,0,2,10,4,1,5,3,6], xs: &[10,4,1,3,7,5,8,6,0,9,2], mirror: true, g: 2, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n697", os: &[1,0,10,8,5,2,3,6,4,9,7,11], xs: &[6,4,3,1,0,7,11,2,8,5,10,9], mirror: true, g: 4, g4: 2, tau: -2, g0_type: T2_5 },
    K { name: "12n739", os: &[3,5,7,10,0,1,2,6,4,9,8], xs: &[9,2,1,3,4,5,8,0,7,6,10], mirror: true, g: 4, g4: 2, tau: 2, g0_type: T2_5 },
    K { name: "12n801", os: &[1,7,0,5,8,3,6,2,4,9], xs: &[5,3,4,2,1,9,0,7,8,6], mirror: false, g: 2, g4: 2, tau: 2, g0_type: T2_5 },

    // T2,3;2,1
    K { name: "12n156", os: &[8,9,5,6,1,0,10,7,4,3,2], xs: &[10,7,8,4,5,3,2,1,0,9,6], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "12n425", os: &[10,9,8,2,3,4,5,7,0,1,6], xs: &[8,7,5,6,10,1,2,3,4,9,0], mirror: true, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "12n749", os: &[3,1,0,8,9,4,6,2,5,7], xs: &[6,5,4,2,3,1,0,7,8,9], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_3_2_1 },

    // T3,4
    K { name: "12n59", os: &[9,10,5,3,2,1,11,0,6,8,4,7], xs: &[11,8,9,7,6,5,4,2,1,3,10,0], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n67", os: &[10,6,7,0,8,3,1,2,4,5,9], xs: &[8,9,4,5,2,10,6,7,0,1,3], mirror: true, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n68", os: &[2,1,3,7,0,10,6,5,4,9,8], xs: &[6,5,9,2,4,3,1,0,8,7,10], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n89", os: &[10,9,8,5,0,7,6,3,2,1,4], xs: &[8,6,3,2,4,1,10,0,9,5,7], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n93", os: &[11,10,9,3,6,2,0,8,7,5,4,1], xs: &[9,8,5,7,11,4,3,2,1,10,0,6], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n110", os: &[1,0,10,6,7,3,5,2,4,9,8,11], xs: &[5,4,3,1,2,11,0,6,8,7,10,9], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n134", os: &[10,9,8,5,0,7,6,1,4,3,2], xs: &[8,6,4,1,3,2,10,9,0,7,5], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n138", os: &[9,8,5,4,2,0,11,1,7,10,3,6], xs: &[11,10,9,7,6,4,3,8,2,5,0,1], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n172", os: &[2,1,0,9,5,4,3,6,8,7,10], xs: &[6,4,3,2,1,0,7,10,5,9,8], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n207", os: &[9,10,6,4,2,1,8,0,11,5,7,3], xs: &[11,8,9,7,6,5,4,3,1,0,2,10], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n217", os: &[2,1,10,0,7,5,6,3,4,9,8,11], xs: &[5,4,2,3,1,11,0,7,8,6,10,9], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n220", os: &[9,10,5,3,2,1,0,11,6,8,4,7], xs: &[11,8,9,7,6,5,4,2,0,3,10,1], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n228", os: &[9,10,6,4,2,1,11,0,8,5,7,3], xs: &[11,8,9,7,6,5,3,4,1,0,2,10], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n229", os: &[1,0,7,6,10,3,5,4,2,9,8], xs: &[9,5,3,2,4,8,1,0,7,6,10], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n251", os: &[2,1,0,10,5,4,3,7,6,9,8], xs: &[9,4,3,2,1,0,6,5,8,7,10], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n329", os: &[0,3,1,5,2,4,6,7,8,9,10], xs: &[4,9,6,0,10,1,3,2,5,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n366", os: &[0,2,1,4,3,5,7,6,8,9,10], xs: &[3,9,6,0,10,1,4,2,5,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n368", os: &[1,0,9,5,4,6,3,2,8,7,10], xs: &[6,3,2,1,7,0,10,5,4,9,8], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n387", os: &[1,0,9,5,4,3,2,6,8,7,10], xs: &[6,3,2,1,0,7,5,10,4,9,8], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n418", os: &[10,9,8,2,5,0,7,6,4,3,1], xs: &[8,7,4,6,10,3,2,1,0,9,5], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n419", os: &[8,9,4,1,10,0,7,6,5,3,2], xs: &[10,7,8,5,3,4,2,1,0,9,6], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n503", os: &[2,1,0,10,7,4,3,6,5,9,8], xs: &[9,6,3,2,1,0,5,4,8,7,10], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n528", os: &[0,3,1,6,2,4,7,5,8,9,10], xs: &[5,9,7,0,10,1,3,2,4,6,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n594", os: &[3,2,1,0,7,10,4,6,5,9,8], xs: &[9,6,4,5,2,3,0,1,8,7,10], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n603", os: &[3,2,1,0,9,5,4,6,8,7], xs: &[8,6,4,3,2,1,7,0,5,9], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n624", os: &[8,9,5,3,1,0,10,7,6,4,2], xs: &[10,7,8,6,5,4,2,3,0,1,9], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n639", os: &[8,7,5,2,1,10,6,0,9,4,3], xs: &[10,9,8,7,4,3,2,5,1,0,6], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n648", os: &[8,6,4,2,0,10,7,3,9,1,5], xs: &[10,9,8,7,6,5,1,0,2,4,3], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n655", os: &[1,9,0,6,5,4,2,3,8,7,10], xs: &[4,2,3,1,0,10,6,7,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n660", os: &[0,3,1,6,2,4,8,5,7,9,10], xs: &[5,7,9,0,10,1,3,2,4,6,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n693", os: &[9,10,5,3,11,1,8,0,7,6,4,2], xs: &[11,8,9,6,4,5,3,2,1,0,10,7], mirror: false, g: 4, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n696", os: &[9,10,5,3,1,0,11,6,8,4,2,7], xs: &[11,8,9,7,6,5,4,2,3,1,10,0], mirror: true, g: 4, g4: 3, tau: -3, g0_type: T3_4 },
    K { name: "12n750", os: &[7,5,3,1,2,0,9,6,8,4], xs: &[9,8,6,4,7,5,3,1,2,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n830", os: &[7,5,4,1,2,0,9,6,8,3], xs: &[9,8,6,5,7,4,3,1,2,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n851", os: &[1,9,0,5,6,3,4,2,8,7,10], xs: &[5,3,4,2,1,0,10,7,6,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },

    // T2,7
    K { name: "12n96", os: &[1,0,7,10,6,2,5,3,4,9,8,11], xs: &[5,3,4,2,1,11,0,7,8,6,10,9], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n148", os: &[0,2,1,5,3,6,4,7,8,9,10], xs: &[4,9,6,0,10,2,1,3,5,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n203", os: &[1,0,7,10,6,2,5,3,4,9,8,11], xs: &[5,4,2,3,1,11,0,7,8,6,10,9], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n276", os: &[1,3,0,2,6,4,7,5,8,9,10], xs: &[5,9,7,10,1,0,3,2,4,6,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n402", os: &[0,2,1,4,3,5,6,7,8,9,10], xs: &[3,9,5,0,10,2,1,4,6,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n638", os: &[0,6,9,4,2,5,3,1,8,7,10], xs: &[4,3,2,1,10,0,7,6,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n644", os: &[1,0,6,9,5,4,2,3,8,7,10], xs: &[4,3,1,2,0,10,6,7,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n764", os: &[0,6,9,5,3,4,2,1,8,7,10], xs: &[3,2,4,1,10,0,7,6,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },

    // T3,5
    K { name: "12n91", os: &[9,8,6,4,3,1,11,2,0,7,10,5], xs: &[11,10,9,8,7,5,3,6,4,1,2,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n105", os: &[9,8,7,5,3,2,11,1,0,6,10,4], xs: &[11,10,9,8,7,6,4,5,2,1,3,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n136", os: &[9,8,6,4,2,11,3,1,0,7,10,5], xs: &[11,10,9,8,7,5,6,4,2,1,3,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n153", os: &[8,7,5,10,3,2,1,0,9,6,4], xs: &[10,9,8,6,7,5,4,3,2,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n187", os: &[9,8,6,4,3,11,2,1,7,0,10,5], xs: &[11,10,9,8,7,5,6,3,2,4,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n243", os: &[3,2,1,0,9,6,5,4,8,7], xs: &[8,5,4,3,2,1,0,7,6,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n244", os: &[8,7,6,3,10,2,1,0,9,5,4], xs: &[10,9,8,7,5,6,4,3,2,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n292", os: &[1,0,7,10,6,2,5,4,3,9,8,11], xs: &[5,4,3,2,1,11,0,7,8,6,10,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n328", os: &[7,5,3,2,1,9,0,8,6,4], xs: &[9,8,7,6,4,2,5,3,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n374", os: &[0,9,6,4,5,2,3,1,8,7,10], xs: &[4,3,2,1,0,7,10,6,5,9,8], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n386", os: &[8,7,4,10,2,1,0,9,6,5,3], xs: &[10,9,8,6,7,5,4,3,1,2,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n417", os: &[8,6,4,3,10,2,1,7,0,9,5], xs: &[10,9,8,7,5,6,4,2,3,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n426", os: &[7,5,3,2,9,1,8,0,6,4], xs: &[9,8,7,6,4,5,3,2,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n473", os: &[1,0,9,6,5,4,2,3,8,7,10], xs: &[4,3,2,1,0,10,6,7,5,9,8], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n502", os: &[2,1,0,9,6,5,4,3,8,7], xs: &[8,4,3,2,1,0,7,6,5,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n518", os: &[9,7,5,4,3,11,2,1,8,0,10,6], xs: &[11,10,9,8,6,4,7,3,2,5,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n575", os: &[1,0,9,5,6,3,4,2,8,7,10], xs: &[5,4,3,2,1,0,10,7,6,9,8], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n591", os: &[5,4,2,1,8,0,7,6,3], xs: &[8,7,6,5,3,4,2,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n640", os: &[8,7,5,3,2,10,1,0,6,9,4], xs: &[10,9,8,7,6,4,5,3,1,2,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n647", os: &[8,7,5,3,2,10,1,9,0,6,4], xs: &[10,9,8,7,6,4,5,2,3,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n680", os: &[1,0,10,7,6,3,5,2,4,9,8,11], xs: &[5,4,3,2,1,11,0,6,8,7,10,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n689", os: &[2,1,10,0,7,6,5,3,4,9,8,11], xs: &[5,4,2,3,1,0,11,7,8,6,10,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n691", os: &[1,0,7,10,5,6,3,4,2,9,8,11], xs: &[5,4,3,1,0,2,8,11,7,6,10,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n692", os: &[1,0,10,7,5,6,3,4,2,9,8,11], xs: &[5,3,4,1,0,2,8,11,7,6,10,9], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n694", os: &[8,7,5,2,3,1,10,0,9,6,4], xs: &[10,9,8,7,6,4,2,5,3,1,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },
    K { name: "12n850", os: &[8,7,6,4,2,1,10,0,9,5,3], xs: &[10,9,8,7,6,5,3,4,1,2,0], mirror: false, g: 4, g4: 4, tau: 4, g0_type: T3_5 },

    // T2,3#T3,5
    K { name: "12n242", os: &[5,4,2,1,0,8,7,6,3], xs: &[8,7,6,5,4,3,2,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n472", os: &[8,7,5,3,2,1,10,9,0,6,4], xs: &[10,9,8,7,6,5,4,2,3,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n574", os: &[8,7,6,4,2,1,0,10,9,5,3], xs: &[10,9,8,7,6,5,4,3,1,2,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n679", os: &[9,8,7,5,3,2,11,1,10,0,6,4], xs: &[11,10,9,8,7,6,5,4,2,3,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
    K { name: "12n688", os: &[9,8,6,4,3,2,11,1,7,0,10,5], xs: &[11,10,9,8,7,6,5,3,2,4,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
    K { name: "12n725", os: &[7,5,3,2,1,9,8,0,6,4], xs: &[9,8,7,6,5,4,3,2,1,0], mirror: false, g: 5, g4: 5, tau: 5, g0_type: T2_3_T3_5 },
    K { name: "12n888", os: &[9,7,5,4,3,11,2,1,8,0,10,6], xs: &[11,10,9,8,7,4,6,5,2,3,1,0], mirror: true, g: 5, g4: 5, tau: -5, g0_type: T2_3_T3_5 },
];
