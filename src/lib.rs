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
mod surfaces;
mod utils;

pub use common::{Vec3, Vec4};
pub use curves::{
    BCurve, BCurveData, Circle, Curve, CurveError, Ellipse, Hyperbola, Line, Parabola,
};
pub use splines::{
    eval as spline_eval, eval_diff as spline_eval_diff, eval_diff_all as spline_eval_diff_all,
    knot_multiplicity, multiplicites as knot_multiplicities,
};
pub use surfaces::{Cylinder, Plane, Surface};
