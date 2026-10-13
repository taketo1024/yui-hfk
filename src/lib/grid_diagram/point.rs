use super::MAX_SIZE;
use std::fmt;
use num_integer::Integer;

/// A point on the grid with coordinates (x, y).
/// Used for O-markings, X-markings, and generator positions.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct GridPoint {
    pub x: u8,
    pub y: u8,
}

impl GridPoint {
    pub fn new(x: u8, y: u8) -> Self {
        assert!(x < MAX_SIZE, "x={x} must be < MAX_SIZE={MAX_SIZE}");
        assert!(y < MAX_SIZE, "y={y} must be < MAX_SIZE={MAX_SIZE}");
        GridPoint { x, y }
    }

    /// Shift the point by (dx, dy).
    pub fn shift(&self, dx: i8, dy: i8) -> Self {
        assert!(dx >= 0 || -dx <= self.x as i8);
        assert!(dy >= 0 || -dy <= self.y as i8);
        GridPoint {
            x: (self.x as i8 + dx) as u8,
            y: (self.y as i8 + dy) as u8,
        }
    }

    /// The four corners of a unit square centered at this point.
    pub fn corners(&self) -> [GridPoint; 4] {
        [
            self.shift(1, 1),
            self.shift(-1, 1),
            self.shift(-1, -1),
            self.shift(1, -1),
        ]
    }

    pub fn is_even(&self) -> bool { 
        self.x.is_even() && self.y.is_even()
    }

    pub fn is_odd(&self) -> bool { 
        self.x.is_odd() && self.y.is_odd()
    }
}

impl fmt::Display for GridPoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

impl PartialOrd for GridPoint {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        use std::cmp::Ordering::*;
        match (self.x.cmp(&other.x), self.y.cmp(&other.y)) {
            (Less,    Less)    => Some(Less),
            (Equal,   Equal)   => Some(Equal),
            (Greater, Greater) => Some(Greater),
            _                  => None,
        }
    }
}

/// Convert a sequence to even points: `(2i, 2*seq[i])`.
pub fn even_points(seq: &[u8]) -> Vec<GridPoint> {
    seq.iter()
        .enumerate()
        .map(|(i, &j)| GridPoint::new(2 * i as u8, 2 * j))
        .collect()
}

/// Convert a sequence to odd points: `(2i+1, 2*seq[i]+1)`.
pub fn odd_points(seq: &[u8]) -> Vec<GridPoint> {
    seq.iter()
        .enumerate()
        .map(|(i, &j)| GridPoint::new(2 * i as u8 + 1, 2 * j + 1))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let p = GridPoint::new(3, 5);
        assert_eq!(p.x, 3);
        assert_eq!(p.y, 5);
    }

    #[test]
    #[should_panic]
    fn test_new_x_out_of_bounds() {
        GridPoint::new(MAX_SIZE, 0);
    }

    #[test]
    #[should_panic]
    fn test_new_y_out_of_bounds() {
        GridPoint::new(0, MAX_SIZE);
    }

    #[test]
    fn test_shift() {
        let p = GridPoint::new(3, 5);
        assert_eq!(p.shift(2, -1), GridPoint::new(5, 4));
        assert_eq!(p.shift(-3, -5), GridPoint::new(0, 0));
    }

    #[test]
    fn test_corners() {
        let p = GridPoint::new(3, 5);
        let c = p.corners();
        assert_eq!(c[0], GridPoint::new(4, 6));
        assert_eq!(c[1], GridPoint::new(2, 6));
        assert_eq!(c[2], GridPoint::new(2, 4));
        assert_eq!(c[3], GridPoint::new(4, 4));
    }

    #[test]
    fn test_display() {
        assert_eq!(GridPoint::new(3, 5).to_string(), "(3, 5)");
    }

    #[test]
    fn test_partial_ord() {
        let p = GridPoint::new(1, 2);
        let q = GridPoint::new(2, 3);
        let r = GridPoint::new(2, 1);
        assert!(p < q);           // both components less
        assert!(!(p < r));        // x: 1 < 2 but y: 2 > 1 → incomparable
        assert!(!(r < p));
        assert_eq!(p.partial_cmp(&r), None);
    }
}
