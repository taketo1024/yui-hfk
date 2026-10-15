//! Tests of the helpers in `common`, on small knots.

mod common;

use yui_hfk::G0CalcConfig;
use common::check_via_mirror;
use common::types::{K, T2_7};

// KnotInfo's grid diagram.
#[rustfmt::skip]
static KNOT_7_1: K = K { name: "7_1", os: &[6,5,4,3,2,1,0,8,7], xs: &[8,7,6,5,4,3,2,1,0], mirror: false, g: 3, g4: 3, tau: 3, g0_type: T2_7 };

// With `max_size = 1`, `R(0,3)` is accepted (z lies in it) and `R(1,2)` is left undecided. It is
// decided on the mirror by rejecting R[(-3,0),(0,-3)], the dual of {R(0,3), R(3,0)}.
#[test]
fn t2_7_via_mirror() {
    let config = G0CalcConfig { max_size: 1, ..Default::default() };
    let mirror_config = G0CalcConfig { try_fast_ng: true, ..Default::default() };
    check_via_mirror(&KNOT_7_1, config, mirror_config);
}
