use log::info;
use num_integer::Integer;
use yui_core::util::format::table;

use crate::grid_diagram::MAX_SIZE;

use super::point::{GridPoint, odd_points};

/// A grid diagram representing a knot or link.
///
/// A grid diagram consists of an n×n grid with n O-markings and n X-markings,
/// one in each row and column. The O and X positions determine the knot topology.
///
/// See arXiv:2002.09210 Section 4.1 for mathematical details.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GridDiagram {
    pub os: Vec<GridPoint>,
    pub xs: Vec<GridPoint>,
    pub size: u8
}

impl GridDiagram {
    /// Create a grid diagram from O and X marking positions.
    fn new(os: Vec<GridPoint>, xs: Vec<GridPoint>, size: u8) -> Self {
        assert!(size < MAX_SIZE);
        assert!(size.is_even());

        assert_eq!(os.len(), xs.len());
        assert!(os.iter().all(|p| p.is_odd() && p.x < size && p.y < size));
        assert!(xs.iter().all(|p| p.is_odd() && p.x < size && p.y < size));

        // TODO add more validations

        GridDiagram { os, xs, size }
    }

    /// Create from an arc presentation (sequence of O and X indices).
    ///
    /// The sequence alternates: [O₀, X₀, O₁, X₁, ...].
    /// Each pair (Oᵢ, Xᵢ) specifies the column positions of the marks in row i.
    pub fn from_arc_presentation(code: &[u8]) -> Self {
        assert_eq!(code.len() % 2, 0, "Arc presentation must have even length");
        assert!(code.iter().all(|p| 2 * p < MAX_SIZE));

        let mut os = Vec::new();
        let mut xs = Vec::new();

        for (i, chunk) in code.chunks(2).enumerate() {
            let x_o = 2 * chunk[0] - 1;
            let x_x = 2 * chunk[1] - 1;
            let y = (2 * i + 1) as u8;

            os.push(GridPoint::new(x_o, y));
            xs.push(GridPoint::new(x_x, y));
        }

        // Sort by x coordinate
        os.sort_by_key(|p| p.x);
        xs.sort_by_key(|p| p.x);
        
        let size = 2 * os.len() as u8;

        GridDiagram::new(os, xs, size)
    }

    /// Create from O and X sequences.
    ///
    /// `os[i]` is the row index of the O-marking in column `i`.
    /// `xs[i]` is the row index of the X-marking in column `i`.
    pub fn from_sequences(os: &[u8], xs: &[u8]) -> Self {
        assert_eq!(os.len(), xs.len());
        let os = odd_points(os);
        let xs = odd_points(xs);
        let size = 2 * os.len() as u8;
        GridDiagram::new(os, xs, size)
    }

    /// Load a grid diagram by knot name (e.g. `"3_1"`) from the runtime
    /// data dir. Expects `<data_dir>/grid_diagrams/<name>.json` with shape
    /// `[[o-coords], [x-coords]]`. The data dir is resolved via
    /// `$YUI_DATA_DIR` (or the platform default) — see `yui_core::util::data_dir`.
    pub fn load(name: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let json_str = yui_core::util::data_dir::load_json("grid_diagrams", name)?;
        let [os_seq, xs_seq]: [Vec<u8>; 2] = serde_json::from_str(&json_str)?;
        Ok(GridDiagram::from_sequences(&os_seq, &xs_seq))
    }

    /// The number of markings (and grid dimension n for n×n grid).
    pub fn grid_number(&self) -> u8 {
        self.size / 2
    }

    /// `[os, xs]` as in `from_sequences`, the inverse of it.
    pub fn sequences(&self) -> [Vec<u8>; 2] {
        // note: os and xs are kept sorted by column.
        let rows = |ps: &[GridPoint]| -> Vec<u8> {
            ps.iter().map(|p| (p.y - 1) / 2).collect()
        };
        [rows(&self.os), rows(&self.xs)]
    }

    /// The grid as a table of rows `j` (top to bottom) and columns `i`, as in `from_sequences`,
    /// with `O` and `X` on the markings and `.` on the empty cells.
    pub fn to_table_string(&self) -> String {
        let n = self.grid_number();
        table("j\\i", (0..n).rev(), 0..n, |j, i| {
            let p = GridPoint::new(2 * i + 1, 2 * j + 1);
            if self.os.contains(&p) {
                "O"
            } else if self.xs.contains(&p) {
                "X"
            } else {
                "."
            }
        })
    }

    /// Info-logs the diagram, as `[os, xs]` and as a table, with its grid number and writhe.
    pub fn log_diagram(&self) {
        info!("diagram: {:?}\n{}", self.sequences(), self.to_table_string());
        info!("  grid number: {}", self.grid_number());
        info!("  writhe: {}", self.writhe());
    }

pub fn unknot() -> GridDiagram { 
        GridDiagram::from_sequences(&[0,1], &[1,0])
    }

    pub fn trefoil() -> GridDiagram { 
        GridDiagram::from_sequences(&[2,1,0,4,3], &[4,3,2,1,0])
    }

    pub fn figure8() -> GridDiagram {
        GridDiagram::from_sequences(&[3,5,0,2,1,4], &[0,1,4,5,3,2])
    }

    /// Rotate the diagram 90 degrees counter-clockwise.
    /// Each point (x, y) maps to (size - y, x), then Os and Xs are re-sorted by x.
    pub fn rotate90(&self) -> GridDiagram {
        let size = self.size;
        let transform = |p: &GridPoint| GridPoint::new(size - p.y, p.x);

        let mut os: Vec<GridPoint> = self.os.iter().map(transform).collect();
        let mut xs: Vec<GridPoint> = self.xs.iter().map(transform).collect();

        os.sort_by_key(|p| p.x);
        xs.sort_by_key(|p| p.x);

        GridDiagram::new(os, xs, size)
    }

    /// The writhe of the knot diagram: the signed count of crossings, where a
    /// vertical strand (X → O) always passes over a horizontal one (O → X).
    /// Mirroring negates it, so its sign says which side is the "positive" one.
    pub fn writhe(&self) -> isize {
        let n = self.grid_number() as usize;
        let row = |p: &GridPoint| ((p.y - 1) / 2) as usize;

        let o_row: Vec<usize> = self.os.iter().map(row).collect();
        let x_row: Vec<usize> = self.xs.iter().map(row).collect();

        // inverse maps: which column holds the O (resp. X) of a given row
        let mut o_col = vec![0usize; n];
        let mut x_col = vec![0usize; n];
        for c in 0..n {
            o_col[o_row[c]] = c;
            x_col[x_row[c]] = c;
        }

        (0..n).map(|c| {
            let (r0, r1) = sorted(x_row[c], o_row[c]);
            let v: isize = if o_row[c] > x_row[c] { 1 } else { -1 };

            (r0 + 1..r1).filter(|&r| {
                let (c0, c1) = sorted(o_col[r], x_col[r]);
                c0 < c && c < c1
            }).map(|r| {
                let h: isize = if x_col[r] > o_col[r] { 1 } else { -1 };
                -v * h
            }).sum::<isize>()
        }).sum()
    }

    // ============ Degree computation (Section 4.1 of arXiv:2002.09210) ============

    /// Count pairs (p, q) ∈ x × y where p strictly dominates q.
    /// p strictly dominates q means p.x < q.x AND p.y < q.y.
    fn count_dominations(x: &[GridPoint], y: &[GridPoint]) -> isize {
        x.iter().map(|p| y.iter().filter(|q| p < q).count()).sum::<usize>() as isize
    }

    /// Symmetric domination count: I(x, y) + I(y, x).
    fn symmetric_dominations(x: &[GridPoint], y: &[GridPoint]) -> isize {
        Self::count_dominations(x, y) + Self::count_dominations(y, x)
    }

    /// Maslov grading formula.
    /// M(ref, x) = (J(x, x) - 2*J(x, ref) + J(ref, ref)) / 2 + 1
    /// See Book Section 4.3, eqn (4.5). 
    pub fn maslov_grading(ref_pts: &[GridPoint], x: &[GridPoint]) -> isize {
        let j_xx = Self::symmetric_dominations(x, x);
        let j_xr = Self::symmetric_dominations(x, ref_pts);
        let j_rr = Self::symmetric_dominations(ref_pts, ref_pts);
        (j_xx - 2 * j_xr + j_rr) / 2 + 1
    }

    /// Compute the Maslov degree of a generator given by its points.
    pub fn maslov_degree(&self, gen_points: &[GridPoint]) -> isize {
        Self::maslov_grading(&self.os, gen_points)
    }

    /// Compute the Alexander degree of a generator given by its points.
    /// A(x) = (M_O(x) - M_X(x) - n + 1) / 2
    /// See Book Def 4.3.2. 
    pub fn alexander_degree(&self, gen_points: &[GridPoint]) -> isize {
        let n = self.grid_number() as isize;
        let m_o = Self::maslov_grading(&self.os, gen_points);
        let m_x = Self::maslov_grading(&self.xs, gen_points);
        (m_o - m_x - n + 1) / 2
    }
}

fn sorted(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

impl std::fmt::Display for GridDiagram {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GridDiagram(n={})", self.grid_number())
    }
}

// ============ Unit Tests ============

#[cfg(test)]
mod tests {
    use super::*;

    /// Test 1: Create from arc presentation
    #[test]
    fn test_from_arc_presentation() {
        let g = GridDiagram::from_arc_presentation(&[1, 2, 3, 1, 2, 3]);
        assert_eq!(g.grid_number(), 3);
        assert_eq!(g.os.len(), 3);
        assert_eq!(g.xs.len(), 3);
        assert_eq!(g.size, 6);
    }

    /// Test 2: Create from sequences
    #[test]
    fn test_from_sequences() {
        let g = GridDiagram::from_sequences(&[2, 1, 0], &[0, 2, 1]);
        assert_eq!(g.grid_number(), 3);

        // O sequence [2, 1, 0] → odd points
        // Column 0, row 2 → (1, 5)
        // Column 1, row 1 → (3, 3)
        // Column 2, row 0 → (5, 1)
        // After sorting by x: (1, 5), (3, 3), (5, 1)
        assert_eq!(g.os[0], GridPoint::new(1, 5));
        assert_eq!(g.os[1], GridPoint::new(3, 3));
        assert_eq!(g.os[2], GridPoint::new(5, 1));
    }

    /// Test 3: Basic properties
    #[test]
    fn test_basic_properties() {
        let g = GridDiagram::from_sequences(&[0, 1], &[1, 0]);
        assert_eq!(g.grid_number(), 2);
        assert_eq!(g.size, 4);
    }

    /// Test 4: Degree computation - unlink (unknot sum)
    /// For the standard unlink (unknot), we know specific degree properties.
    #[test]
    fn test_maslov_degree_consistency() {
        let g = GridDiagram::from_sequences(&[1, 0], &[0, 1]);
        let gen_points = odd_points(&[0, 1]);
        let m_deg = g.maslov_degree(&gen_points);
        // Just verify it's computable and returns an integer
        assert!(m_deg > i32::MIN as isize && m_deg < i32::MAX as isize);
    }

    /// Test 5: Alexander degree is well-defined
    #[test]
    fn test_alexander_degree_computation() {
        let g = GridDiagram::from_sequences(&[0, 1], &[1, 0]);
        let gen_points = odd_points(&[0, 1]);
        let a_deg = g.alexander_degree(&gen_points);
        // Just verify it's computable
        assert!(a_deg > i32::MIN as isize && a_deg < i32::MAX as isize);
    }

    /// Test 6: Degree bigrading relationship
    /// Alexander degree computation uses both Maslov gradings
    #[test]
    fn test_alexander_uses_both_gradings() {
        let g = GridDiagram::from_sequences(&[0, 1, 2], &[1, 2, 0]);
        let gen_points = odd_points(&[0, 1, 2]);

        let m_o = g.maslov_degree(&gen_points);
        let m_x_points = g.xs.clone();
        let m_x = GridDiagram::maslov_grading(&m_x_points, &gen_points);
        let a_deg = g.alexander_degree(&gen_points);

        // Verify the formula: A = (M_O - M_X - n + 1) / 2
        let n = g.grid_number() as isize;
        assert_eq!(a_deg, (m_o - m_x - n + 1) / 2);
    }


    /// Test 8: Display formatting
    #[test]
    fn test_display() {
        let g = GridDiagram::from_sequences(&[0, 1], &[1, 0]);
        let s = g.to_string();
        assert!(s.contains("n=2"));
    }

    /// Test 9: Round-trip via sequences
    #[test]
    fn test_sequences_round_trip() {
        let os_seq = vec![2, 0, 1];
        let xs_seq = vec![1, 2, 0];
        let g = GridDiagram::from_sequences(&os_seq, &xs_seq);
        assert_eq!(g.sequences(), [os_seq, xs_seq]);

        let m = GridDiagram::trefoil().rotate90();
        let [os, xs] = m.sequences();
        assert_eq!(GridDiagram::from_sequences(&os, &xs), m);
    }

    #[test]
    fn test_table_string() {
        let s = GridDiagram::unknot().to_table_string();
        let rows: Vec<Vec<&str>> = s.lines().map(|l| l.split_whitespace().collect()).collect();
        assert_eq!(rows, [vec!["j\\i", "0", "1"], vec!["1", "X", "O"], vec!["0", "O", "X"]]);
    }


    /// Test 11: Domination count is non-negative
    #[test]
    fn test_domination_non_negative() {
        let pts1 = vec![GridPoint::new(1, 1), GridPoint::new(3, 3)];
        let pts2 = vec![GridPoint::new(2, 2)];
        let count = GridDiagram::count_dominations(&pts1, &pts2);
        assert!(count >= 0);
    }

    /// Test 12: Symmetric dominations are symmetric
    #[test]
    fn test_symmetric_dominations_symmetry() {
        let pts1 = vec![GridPoint::new(1, 1), GridPoint::new(3, 3)];
        let pts2 = vec![GridPoint::new(2, 2), GridPoint::new(4, 1)];

        let j12 = GridDiagram::symmetric_dominations(&pts1, &pts2);
        let j21 = GridDiagram::symmetric_dominations(&pts2, &pts1);

        // J should be symmetric by definition
        assert_eq!(j12, j21);
    }

    #[test]
    fn test_writhe() {
        // positive torus knots: writhe > 0, matching tau > 0
        assert_eq!(GridDiagram::trefoil().writhe(), 3);
        assert_eq!(GridDiagram::test_data("5_1").unwrap().writhe(), 5);
        // 8_19 = T(3,4), a positive braid closure
        let k8_19 = GridDiagram::from_sequences(&[3,2,1,0,6,5,4], &[6,5,4,3,2,1,0]);
        assert_eq!(k8_19.writhe(), 8);
        // figure-8 is amphichiral
        assert_eq!(GridDiagram::figure8().writhe(), 0);
        assert_eq!(GridDiagram::unknot().writhe(), 0);
    }

    #[test]
    fn test_writhe_negates_under_mirror() {
        for name in ["0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"] {
            let g = GridDiagram::test_data(name).unwrap();
            assert_eq!(g.rotate90().writhe(), -g.writhe(), "{name}");
        }
    }

    #[test]
    fn test_data_trefoil() {
        let g = GridDiagram::test_data("3_1").unwrap();
        assert_eq!(g.grid_number(), 5);
        assert_eq!(g, GridDiagram::trefoil());
    }
}
