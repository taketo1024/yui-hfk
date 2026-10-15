//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T2,7 ([SS] (v2), Tables 1 and 2).
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! With the config in `check` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided. Those of this type that need a larger `max_size` (12n96,
//! 12n203) are in `thick_12_heavy.rs`.
//!
//! The test is slow:
//!
//!     cargo test --release --test thick_12_t27 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T2_7};

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
    K { name: "12n148", os: &[0,2,1,5,3,6,4,7,8,9,10], xs: &[4,9,6,0,10,2,1,3,5,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n276", os: &[1,3,0,2,6,4,7,5,8,9,10], xs: &[5,9,7,10,1,0,3,2,4,6,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n402", os: &[0,2,1,4,3,5,6,7,8,9,10], xs: &[3,9,5,0,10,2,1,4,6,7,8], mirror: true, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n638", os: &[0,6,9,4,2,5,3,1,8,7,10], xs: &[4,3,2,1,10,0,7,6,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n644", os: &[1,0,6,9,5,4,2,3,8,7,10], xs: &[4,3,1,2,0,10,6,7,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
    K { name: "12n764", os: &[0,6,9,5,3,4,2,1,8,7,10], xs: &[3,2,4,1,10,0,7,6,5,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 },
];
