//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::{zero_vec, Vector},
    Curve,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Line {
    origin: Vector,
    direction: Vector,
}

impl Line {
    pub fn new(
        origin: Vector,
        direction: Vector,
    ) -> Self {
        assert_eq!(origin.len(), direction.len());
        Self { origin, direction }
    }
}

impl Curve for Line {
    fn dimension(&self) -> usize {
        self.origin.len()
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vector {
        let c = (&self.origin + u * &self.direction).into();
        c
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vector {
        if m == 0 {
            self.eval_point(u)
        } else if m == 1 {
            self.direction.clone()
        } else {
            zero_vec(self.dimension())
        }
    }
}
