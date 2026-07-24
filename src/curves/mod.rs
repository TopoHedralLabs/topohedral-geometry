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

pub use bcurve::{BCurve, BCurveData, CurveError};
pub use common::Curve;
