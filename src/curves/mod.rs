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
mod hyperbola;
mod line;
mod parabola;

pub use bcurve::{BCurve, BCurveData, CurveError};
pub use circle::Circle;
pub use common::Curve;
pub use ellipse::Ellipse;
pub use hyperbola::Hyperbola;
pub use line::Line;
pub use parabola::Parabola;
