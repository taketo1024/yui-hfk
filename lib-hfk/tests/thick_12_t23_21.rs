//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T2,3;2,1 ([SS] (v2), Tables 1 and 2).
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
//!     cargo test --release --test thick_12_t23_21 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T2_3_2_1};

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
    K { name: "12n156", os: &[8,9,5,6,1,0,10,7,4,3,2], xs: &[10,7,8,4,5,3,2,1,0,9,6], mirror: false, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "12n425", os: &[10,9,8,2,3,4,5,7,0,1,6], xs: &[8,7,5,6,10,1,2,3,4,9,0], mirror: true, g: 4, g4: 2, tau: 2, g0_type: T2_3_2_1 },
    K { name: "12n749", os: &[3,1,0,8,9,4,6,2,5,7], xs: &[6,5,4,2,3,1,0,7,8,9], mirror: false, g: 3, g4: 2, tau: 2, g0_type: T2_3_2_1 },
];
