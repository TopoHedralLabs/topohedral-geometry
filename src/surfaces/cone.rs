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

use topohedral_linalg::VectorOps;

use crate::{
    utils::{is_unit_vector, orthogonal_vector},
    Surface, Vec3,
};

pub struct Cone {
    center: Vec3,
    axis: Vec3,
    radius_bottom: f64,
    tan_alpha: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Cone {
    pub fn new(
        center: Vec3,
        axis: Vec3,
        height: f64,
        radius_bottom: f64,
        radius_top: f64,
    ) -> Self {
        assert!(is_unit_vector(&axis));
        assert!(radius_top > 0.0);
        assert!(height > 0.0);
        assert!(radius_bottom > radius_top);
        let tan_alpha = (radius_bottom - radius_top) / height;
        let x_axis = orthogonal_vector(&axis);
        let y_axis = x_axis.cross(&axis);
        Self {
            center,
            axis,
            radius_bottom,
            tan_alpha,
            x_axis,
            y_axis,
        }
    }

    fn radius(
        &self,
        v: f64,
    ) -> f64 {
        self.radius_bottom - v * self.tan_alpha
    }
}

impl Surface for Cone {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            axis: ez,
            radius_bottom: _,
            tan_alpha: _,
            x_axis: ex,
            y_axis: ey,
        } = self;
        let r = self.radius(v);
        (c + r * (u.cos() * ex + u.sin() * ey) + v * ez).into()
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

        if n == 0 || n == 1 {
            let Self {
                center: _,
                axis: _,
                radius_bottom: _,
                tan_alpha: _,
                x_axis: ex,
                y_axis: ey,
            } = self;

            let v_const = if n == 0 {
                self.radius(v)
            } else {
                -self.tan_alpha
            };

            if m.is_multiple_of(2) {
                let is_odd = ((m / 2) % 2) == 1;
                let sign = if is_odd { -1.0 } else { 1.0 };
                return (sign * v_const * (u.cos() * ex + u.sin() * ey)).into();
            } else {
                let is_odd = (((m - 1) / 2) % 2) == 1;
                let sign = if is_odd { -1.0 } else { 1.0 };
                return (sign * v_const * (-u.sin() * ex + u.cos() * ey)).into();
            }
        }

        Vec3::zeros()
    }
}
