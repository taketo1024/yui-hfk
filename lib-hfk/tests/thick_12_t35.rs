//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T3,5 ([SS] (v2), Tables 1 and 2).
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
//!     cargo test --release --test thick_12_t35 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T3_5};

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
];
