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
mod utils;

pub use curves::{BCurve, BCurveData, Circle, Curve, CurveError, Line};
pub use splines::{
    eval as spline_eval, eval_diff as spline_eval_diff, eval_diff_all as spline_eval_diff_all,
    knot_multiplicity, multiplicites as knot_multiplicities,
};
