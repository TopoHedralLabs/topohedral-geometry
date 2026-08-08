//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    utils::{is_unit_vector, orthogonal_vector},
    Surface, Vec3,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Cylinder {
    center: Vec3,
    axis: Vec3,
    radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Cylinder {
    pub fn new(
        center: Vec3,
        axis: Vec3,
        radius: f64,
    ) -> Self {
        assert!(radius > 0.0);
        assert!(is_unit_vector(&axis));

        let x_axis = orthogonal_vector(&axis);
        let y_axis = x_axis.cross(&axis);

        Self {
            center,
            axis,
            radius,
            x_axis,
            y_axis,
        }
    }
}

impl Surface for Cylinder {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            axis: ez,
            radius: r,
            x_axis: ex,
            y_axis: ey,
        } = self;
        (c + r * u.cos() * ex + r * u.sin() * ey + v * ez).into()
    }

    fn eval_derivative(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
    ) -> Vec3 {
        if m == 0 && n == 0 {
            return self.eval_point(u, v);
        }

        if m == 0 && n == 1 {
            return self.axis;
        }

        if m > 0 && n == 0 {
            let Self {
                center: _,
                axis: _,
                radius: r,
                x_axis: ex,
                y_axis: ey,
            } = self;

            if m.is_multiple_of(2) {
                let is_odd = ((m / 2) % 2) == 1;
                let sign = if is_odd { -1.0 } else { 1.0 };
                return (sign * (r * u.cos() * ex + r * u.sin() * ey)).into();
            } else {
                let is_odd = ((m / 2) % 2) == 1;
                let sign = if is_odd { -1.0 } else { 1.0 };
                return (sign * (-r * u.sin() * ex + r * u.cos() * ey)).into();
            }
        }
        // all other configurations, mixed derivatives, v derivative > 1 all zero
        Vec3::zeros()
    }
}
