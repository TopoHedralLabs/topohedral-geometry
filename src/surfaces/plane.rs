//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use approx::assert_relative_eq;
use topohedral_linalg::VectorOps;

//{{{ crate imports
use crate::common::Vec3;
use crate::surfaces::common::Surface;
use crate::utils::orthogonal_vector;
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Plane {
    origin: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Plane {
    pub fn new(
        origin: Vec3,
        normal: Vec3,
    ) -> Self {
        assert_relative_eq!(normal.norm(), 1.0);
        let x_axis = orthogonal_vector(&normal);
        let y_axis = x_axis.cross(&normal);

        Self {
            origin,
            x_axis,
            y_axis,
        }
    }
}

impl Surface for Plane {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            origin: o,
            x_axis,
            y_axis,
        } = self;
        (o + u * x_axis + v * y_axis).into()
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

        if m == 1 && n == 0 {
            return self.x_axis;
        }

        if m == 0 && n == 1 {
            return self.y_axis;
        }

        Vec3::zeros()
    }
}
