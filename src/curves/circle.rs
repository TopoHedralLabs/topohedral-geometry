//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::common::Vec3;
use crate::utils::orthogonal_vector;
use crate::Curve;
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Circle {
    center: Vec3,
    radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Circle {
    pub fn new(
        center: Vec3,
        normal: Vec3,
        radius: f64,
    ) -> Self {
        assert!(normal.norm() > 0.0);
        assert!(radius > 0.0);
        let x_axis = orthogonal_vector(&normal);
        let y_axis = x_axis.cross(&normal);
        Self {
            center,
            radius,
            x_axis,
            y_axis,
        }
    }
}

impl Curve for Circle {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            radius: r,
            x_axis,
            y_axis,
        } = self;

        let out = (c + r * u.cos() * x_axis + r * u.sin() * y_axis).into();
        out
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            return self.eval_point(u);
        }

        let Self {
            center: _,
            radius: r,
            x_axis,
            y_axis,
        } = self;
        if m % 2 == 0 {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            let deriv = (sign * (r * u.cos() * x_axis + r * u.sin() * y_axis)).into();
            deriv
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            let deriv = (sign * (-r * u.sin() * x_axis + r * u.cos() * y_axis)).into();
            deriv
        }
    }
}
