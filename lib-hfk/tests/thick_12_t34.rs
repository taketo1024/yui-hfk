//! `G0(K)` for the prime knots with 12 crossings that are homologically thick and
//! non-slice, of `G0`-type T3,4 ([SS] (v2), Tables 1 and 2).
//!
//! The grid diagrams are KnotInfo's. Their grid numbers equal the arc index, so they
//! cannot be simplified. `g`, `g4` and `|tau|` are also from KnotInfo. Each knot is
//! tested on the side with the smaller linear systems (see `K`).
//!
//! With the config in `check` (`try_fast_ng`, `max_size = 2_000_000`), every region of
//! every knot is decided. Those of this type that need a larger `max_size` (12n59,
//! 12n220, 12n228, 12n693, 12n696) are in `thick_12_heavy.rs`.
//!
//! The test is slow:
//!
//!     cargo test --release --test thick_12_t34 -- --ignored --nocapture

mod common;

use yui_hfk::G0CalcConfig;
use common::run_all;
use common::types::{K, T3_4};

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
    K { name: "12n750", os: &[7,5,3,1,2,0,9,6,8,4], xs: &[9,8,6,4,7,5,3,1,2,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n830", os: &[7,5,4,1,2,0,9,6,8,3], xs: &[9,8,6,5,7,4,3,1,2,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
    K { name: "12n851", os: &[1,9,0,5,6,3,4,2,8,7,10], xs: &[5,3,4,2,1,0,10,7,6,9,8], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T3_4 },
];
