//! Hard-coded grid diagrams for knots up to 6 crossings, used only in
//! `cargo test` (and downstream crates that opt into the `test-utils`
//! feature). The OX-encoding here is the one already in use by yui-hfk;
//! it does not need to agree with any external database.

use crate::GridDiagram;

impl GridDiagram {
    pub fn test_data(name: &str) -> Result<Self, String> {
        let (os, xs): (&[u8], &[u8]) = match name {
            "0_1" => (&[0, 1], &[1, 0]),
            "3_1" => (&[2, 1, 0, 4, 3], &[4, 3, 2, 1, 0]),
            "4_1" => (&[0, 1, 4, 5, 3, 2], &[3, 5, 0, 2, 1, 4]),
            "5_1" => (&[5, 4, 3, 2, 1, 0, 6], &[0, 6, 5, 4, 3, 2, 1]),
            "5_2" => (&[4, 2, 3, 1, 0, 6, 5], &[0, 5, 6, 4, 2, 1, 3]),
            "6_1" => (&[6, 4, 5, 1, 3, 7, 2, 0], &[3, 7, 2, 6, 0, 1, 4, 5]),
            "6_2" => (&[1, 3, 2, 6, 0, 7, 5, 4], &[5, 7, 4, 1, 3, 2, 0, 6]),
            "6_3" => (&[1, 4, 6, 7, 5, 3, 0, 2], &[3, 7, 0, 2, 1, 6, 4, 5]),
            _ => return Err(format!("no test data for `{name}`")),
        };
        Ok(GridDiagram::from_sequences(os, xs))
    }
}
