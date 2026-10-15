//! Tests of the helpers in `common`, on small knots.

mod common;

use yui_hfk::G0CalcConfig;
use yui_hfk::g_zero::Region;
use common::{check_via_mirror, mirror_tests};
use common::types::{K, Regions, T2_7, to_regions};

// KnotInfo's grid diagram of 7_1 = T(2,7), of G0-type T2,7.
#[rustfmt::skip]
static KNOT_7_1: K = K { name: "7_1", os: &[6,5,4,3,2,1,0,8,7], xs: &[8,7,6,5,4,3,2,1,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 };

// The region of the dual type of T2,7.
const S: &str = "R[(-3, 0), (-2, -1), (-1, -2), (0, -3)]";

#[test]
fn expected_is_the_dual_type_where_tau_is_negative() {
    let t: Regions = T2_7.iter().map(|r| r.iter().copied().collect()).collect();
    assert_eq!(KNOT_7_1.expected(), t);

    let mirror = K { mirror: true, tau: -3, ..KNOT_7_1 };
    let s: Regions = [[(-3, 0), (-2, -1), (-1, -2), (0, -3)].into()].into();
    assert_eq!(mirror.expected(), s);
}

// With `max_size = 1`, the expected R(1,2) and R(2,1) are left undecided, which is tolerated.
#[test]
fn check_g0_counts_undecided() {
    let config = G0CalcConfig { max_size: 1, ..Default::default() };
    assert_eq!(KNOT_7_1.check(config), 2);
}

// With `max_size = 1`, `R(0,3)` is accepted (z lies in it) and `R(1,2)` is left undecided. It is
// decided on the mirror by rejecting R[(-3,0),(0,-3)], as in `mirror_tests_to_reject`.
#[test]
fn t2_7_via_mirror() {
    let config = G0CalcConfig { max_size: 1, ..Default::default() };
    let mirror_config = G0CalcConfig { try_fast_ng: true, ..Default::default() };
    check_via_mirror(&KNOT_7_1, config, mirror_config);
}

// R[(0,2),(1,1),(2,0)] is left undecided besides T2,7 (as for 12n96): accept S on the mirror.
#[test]
fn mirror_tests_to_accept() {
    let t = t2_7();
    let u = [region("R[(0, 2), (1, 1), (2, 0)]")];
    assert_eq!(mirror_tests(&t, &u, &t), (vec![region(S)], vec![]));
}

// R(1,2) and R(2,1) are left undecided: reject R[(-3,0),(0,-3)], the dual of {R(0,3), R(3,0)}.
#[test]
fn mirror_tests_to_reject() {
    let a = [region("R(0, 3)"), region("R(3, 0)")];
    let u = [region("R(1, 2)"), region("R(2, 1)")];
    assert_eq!(mirror_tests(&a, &u, &t2_7()), (vec![], vec![region("R[(-3, 0), (0, -3)]")]));
}

// An undecided region that contains accepted ones changes nothing.
#[test]
fn mirror_tests_none() {
    let t = t2_7();
    let u = [region("R[(0, 3), (1, 2), (2, 1), (3, 0)]")];
    assert_eq!(mirror_tests(&t, &u, &t), (vec![], vec![]));
}

fn t2_7() -> Vec<Region> {
    to_regions(&KNOT_7_1.expected())
}

fn region(s: &str) -> Region {
    s.parse().unwrap()
}
