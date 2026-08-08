//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

mod bcurve;
mod circle;
mod common;
mod ellipse;
mod line;

pub use bcurve::{BCurve, BCurveData, CurveError};
pub use circle::Circle;
pub use common::Curve;
pub use ellipse::Ellipse;
pub use line::Line;
