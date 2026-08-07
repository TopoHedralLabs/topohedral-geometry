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
mod common;
mod line;

pub use bcurve::{BCurve, BCurveData, CurveError};
pub use common::Curve;
pub use line::Line;
