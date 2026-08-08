//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::common::Vec3;
use crate::utils::{is_unit_vector, orthogonal_vector};
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
        assert!(is_unit_vector(&normal));
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
            x_axis: ex,
            y_axis: ey,
        } = self;

        (c + r * u.cos() * ex + r * u.sin() * ey).into()
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
            x_axis: ex,
            y_axis: ey,
        } = self;
        if m.is_multiple_of(2) {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (r * u.cos() * ex + r * u.sin() * ey)).into()
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (-r * u.sin() * ex + r * u.cos() * ey)).into()
        }
    }
}
