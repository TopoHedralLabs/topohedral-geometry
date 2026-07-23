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

mod common;
mod curves;
mod splines;

pub use curves::{BCurve, BCurveData, Curve, CurveError};
pub use splines::{
    eval as spline_eval, eval_diff as spline_eval_diff, eval_diff_all as spline_eval_diff_all,
    multiplicites as spline_multiplicities,
};
