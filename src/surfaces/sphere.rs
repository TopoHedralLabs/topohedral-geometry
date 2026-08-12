//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use topohedral_linalg::VectorOps;

use crate::Surface;
//{{{ crate imports
use crate::common::Vec3;
use crate::utils::{cos_derivative, is_unit_vector, orthogonal_vector, sin_derivative};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Sphere {
    center: Vec3,
    radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
    z_axis: Vec3,
}

impl Sphere {
    pub fn new(
        center: Vec3,
        radius: f64,
        polar_axis: Vec3,
    ) -> Self {
        assert!(is_unit_vector(&polar_axis));
        let x_axis = orthogonal_vector(&polar_axis);
        let y_axis = x_axis.cross(&polar_axis);
        Self {
            center,
            radius,
            x_axis,
            y_axis,
            z_axis: polar_axis,
        }
    }

    fn eval_latitude(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        self.x_axis * cos_derivative(u, m) + self.y_axis * sin_derivative(u, m)
    }
}

impl Surface for Sphere {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            radius,
            x_axis: ex,
            y_axis: ey,
            z_axis: ez,
        } = self;
        let r = *radius;
        (c + r * v.cos() * self.eval_latitude(u, 0) + r * v.sin() * ez).into()
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
        let r = self.radius;
        let ez = self.z_axis;
        r * (cos_derivative(v, n) * self.eval_latitude(u, m) + sin_derivative(v, n) * ez)
    }
}
