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

pub struct Ellipse {
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    x_radius: f64,
    y_radius: f64,
}

impl Ellipse {
    pub fn new(
        center: Vec3,
        normal: Vec3,
        x_radius: f64,
        y_radius: f64,
    ) -> Self {
        assert!(is_unit_vector(&normal));
        assert!(x_radius > 0.0);
        assert!(y_radius > 0.0);

        let x_axis = orthogonal_vector(&normal);
        let y_axis = x_axis.cross(&normal);
        Self {
            center,
            x_axis,
            y_axis,
            x_radius,
            y_radius,
        }
    }
}

impl Curve for Ellipse {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            x_axis,
            y_axis,
            x_radius: rx,
            y_radius: ry,
        } = self;

        (c + rx * u.cos() * x_axis + ry * u.sin() * y_axis).into()
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        let Self {
            center: _,
            x_axis,
            y_axis,
            x_radius: rx,
            y_radius: ry,
        } = self;

        if m == 0 {
            return self.eval_point(u);
        }

        if m.is_multiple_of(2) {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (rx * u.cos() * x_axis + ry * u.sin() * y_axis)).into()
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (-rx * u.sin() * x_axis + ry * u.cos() * y_axis)).into()
        }
    }
}
