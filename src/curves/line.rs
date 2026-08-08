//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{common::Vec3, utils::is_unit_vector, Curve};
//}}}
//{{{ std imports
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Line {
    origin: Vec3,
    direction: Vec3,
}

impl Line {
    pub fn new(
        origin: Vec3,
        direction: Vec3,
    ) -> Self {
        assert!(is_unit_vector(&direction));
        Self { origin, direction }
    }
}

impl Curve for Line {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        self.origin + u * self.direction
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            self.eval_point(u)
        } else if m == 1 {
            self.direction
        } else {
            Vec3::zeros()
        }
    }
}
