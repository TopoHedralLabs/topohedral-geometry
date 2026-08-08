//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::common::Vector;
use crate::utils::orthogonal_vector;
use crate::Curve;
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Ellipse {
    center: Vector,
    x_axis: Vector,
    y_axis: Vector,
    x_radius: f64,
    y_radius: f64,
}

impl Ellipse {
    pub fn new(
        center: Vector,
        normal: Vector,
        x_radius: f64,
        y_radius: f64,
    ) -> Self {
        assert_eq!(center.len(), normal.len());
        assert_eq!(center.len(), 3);
        assert!(normal.norm() > 0.0);
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
    fn dimension(&self) -> usize {
        self.x_axis.len()
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vector {
        let Self {
            center: c,
            x_axis,
            y_axis,
            x_radius: rx,
            y_radius: ry,
        } = self;

        let x = (c + rx * u.cos() * x_axis + ry * u.sin() * y_axis).into();
        x
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vector {
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

        if m % 2 == 0 {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            let deriv = (sign * (rx * u.cos() * x_axis + ry * u.sin() * y_axis)).into();
            deriv
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            let deriv = (sign * (-rx * u.sin() * x_axis + ry * u.cos() * y_axis)).into();
            deriv
        }
    }
}
