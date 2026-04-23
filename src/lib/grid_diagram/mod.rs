pub mod point;
pub mod rect;
pub mod diagram;

pub use point::GridPoint;
pub use rect::GridRect;
pub use diagram::GridDiagram;

pub const MAX_SIZE: u8 = 128;
