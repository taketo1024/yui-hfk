use super::{GridPoint, MAX_SIZE};
use std::fmt;

/// A rectangle on the toroidal grid.
/// origin: the lower-left corner
/// size: (width, height)
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GridRect {
    pub origin: GridPoint,
    pub size: (u8, u8),
}

impl GridRect {
    pub fn new(origin: GridPoint, size: (u8, u8)) -> Self {
        assert!(size.0 <= MAX_SIZE);
        assert!(size.1 <= MAX_SIZE);

        GridRect { origin, size }
    }

    /// Construct the positive rectangle from p to q (wrapping around the torus).
    pub fn between(p: GridPoint, q: GridPoint, grid_size: u8) -> Self {
        assert!(grid_size <= MAX_SIZE);
        assert!(p <= q || p < q.shift(grid_size as i8, grid_size as i8)); 

        let l = grid_size;
        let size = if p <= q { 
            (q.x - p.x, q.y - p.y)
        } else { 
            (q.x + l - p.x, q.y + l - p.y)
        };
        GridRect::new(p, size)
    }

    /// The destination point (where the rectangle ends).
    /// May exceed grid_size; not a valid grid coordinate.
    pub fn end(&self) -> GridPoint {
        GridPoint {
            x: self.origin.x + self.size.0,
            y: self.origin.y + self.size.1,
        }
    }

    /// Check if a point is contained in this rectangle.
    /// interior: if true, check strict containment (exclude boundary)
    pub fn contains(&self, p: GridPoint, grid_size: u8, interior: bool) -> bool {
        let ox = self.origin.x;
        let oy = self.origin.y;
        let bx = ox + self.size.0;
        let by = oy + self.size.1;

        in_range(p.x, ox, bx, grid_size, interior) && 
        in_range(p.y, oy, by, grid_size, interior)
    }

    /// Check if any point in the list is in this rectangle.
    pub fn intersects(&self, points: &[GridPoint], grid_size: u8, interior: bool) -> bool {
        points.iter().any(|&p| self.contains(p, grid_size, interior))
    }

    /// Count how many points in the list are in this rectangle.
    pub fn count_intersections(&self, points: &[GridPoint], grid_size: u8, interior: bool) -> usize {
        points.iter().filter(|&&p| self.contains(p, grid_size, interior)).count()
    }
}

impl fmt::Display for GridRect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[origin: {}, size: ({}, {})]", self.origin, self.size.0, self.size.1)
    }
}

fn in_range(p: u8, a: u8, b: u8, grid_size: u8, interior: bool) -> bool {
    if interior { 
        if a < p { p < b } else { p + grid_size < b }
    } else { 
        if a <= p { p <= b } else { p + grid_size <= b }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contains_no_wrap() {
        let rect = GridRect::new(GridPoint::new(2, 2), (4, 4));
        // Rectangle from (2,2) to (6,6), no wrapping
        assert!(rect.contains(GridPoint::new(3, 3), 8, false));
        assert!(rect.contains(GridPoint::new(2, 2), 8, false));
        assert!(rect.contains(GridPoint::new(6, 6), 8, false));
        assert!(!rect.contains(GridPoint::new(1, 3), 8, false));
        assert!(!rect.contains(GridPoint::new(7, 3), 8, false));
    }

    #[test]
    fn test_between() {
        let p = GridPoint::new(2, 2);
        let q = GridPoint::new(5, 5);
        let rect = GridRect::between(p, q, 8);
        // Should have size (3, 3)
        assert_eq!(rect.size, (3, 3));
        assert_eq!(rect.end(), GridPoint::new(5, 5));
    }

    #[test]
    fn test_between_wrap() {
        let p = GridPoint::new(6, 6);
        let q = GridPoint::new(2, 2);
        let rect = GridRect::between(p, q, 8);
        // Should wrap: size = (8 + 2 - 6) % 8 = 4
        assert_eq!(rect.size, (4, 4));
    }

    #[test]
    fn test_end() {
        let rect = GridRect::new(GridPoint::new(2, 3), (4, 5));
        assert_eq!(rect.end(), GridPoint { x: 6, y: 8 });
    }

    #[test]
    fn test_contains_interior() {
        let rect = GridRect::new(GridPoint::new(2, 2), (4, 4));
        // boundary points are not strictly inside
        assert!(!rect.contains(GridPoint::new(2, 2), 8, true));
        assert!(!rect.contains(GridPoint::new(6, 6), 8, true));
        // interior point is strictly inside
        assert!(rect.contains(GridPoint::new(3, 3), 8, true));
    }

    #[test]
    fn test_contains_wrap() {
        // rect from (6,6), size (4,4): destination (10,10), wraps on an 8-grid
        let rect = GridRect::between(GridPoint::new(6, 6), GridPoint::new(2, 2), 8);
        // (7,7): no wrap needed, inside
        assert!(rect.contains(GridPoint::new(7, 7), 8, false));
        // (1,1): wraps around, on boundary
        assert!(rect.contains(GridPoint::new(1, 1), 8, false));
        // (2,2): on far boundary
        assert!(rect.contains(GridPoint::new(2, 2), 8, false));
        // (2,2) strictly inside? 2+8=10, not < 10 → false
        assert!(!rect.contains(GridPoint::new(2, 2), 8, true));
        // (5,5): outside
        assert!(!rect.contains(GridPoint::new(5, 5), 8, false));
    }

    #[test]
    fn test_intersects() {
        let rect = GridRect::new(GridPoint::new(2, 2), (4, 4));
        let inside = vec![GridPoint::new(3, 3), GridPoint::new(5, 5)];
        let outside = vec![GridPoint::new(1, 1), GridPoint::new(7, 7)];
        assert!(rect.intersects(&inside, 8, false));
        assert!(!rect.intersects(&outside, 8, false));
    }

    #[test]
    fn test_count_intersections() {
        let rect = GridRect::new(GridPoint::new(2, 2), (4, 4));
        let points = vec![
            GridPoint::new(3, 3), // inside
            GridPoint::new(5, 5), // inside
            GridPoint::new(1, 1), // outside
            GridPoint::new(7, 7), // outside
        ];
        assert_eq!(rect.count_intersections(&points, 8, false), 2);
        // boundary excluded in interior mode: (2,2) and (6,6) are boundary
        let boundary = vec![GridPoint::new(2, 2), GridPoint::new(6, 6)];
        assert_eq!(rect.count_intersections(&boundary, 8, false), 2);
        assert_eq!(rect.count_intersections(&boundary, 8, true), 0);
    }
}
