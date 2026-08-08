//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{common::Vec3, utils::is_unit_vector, Curve};
//}}}
//{{{ std imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Parabola {
    origin: Vec3,
    normal: Vec3,
    direction: Vec3,
    slope: f64,
}

impl Parabola {
    pub fn new(
        origin: Vec3,
        normal: Vec3,
        direction: Vec3,
        slope: f64,
    ) -> Self {
        assert!(is_unit_vector(&normal));
        assert!(is_unit_vector(&direction));
        assert!(slope > 0.0);
        Self {
            origin,
            normal,
            direction,
            slope,
        }
    }
}

impl Curve for Parabola {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            origin: o,
            normal: n,
            direction: dir_x,
            slope: a,
        } = self;
        let dir_y = dir_x.cross(n);

        (o + 4.0 * a * u * dir_x + u * u * dir_y).into()
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            return self.eval_point(u);
        }

        if m > 2 {
            return Vec3::zeros();
        }

        let Self {
            origin: _,
            normal: n,
            direction: dir_x,
            slope: a,
        } = self;
        let dir_y = dir_x.cross(n);

        if m == 1 {
            let dx_du = (4.0 * a * dir_x + 2.0 * u * dir_y).into();
            return dx_du;
        }

        2.0 * dir_y
    }
}
