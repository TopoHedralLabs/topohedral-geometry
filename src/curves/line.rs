//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::{zero_vec, Vec3},
    Curve,
};
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
        Self { origin, direction }
    }
}

impl Curve for Line {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let c = (&self.origin + u * &self.direction).into();
        c
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            self.eval_point(u)
        } else if m == 1 {
            self.direction.clone()
        } else {
            zero_vec()
        }
    }
}
