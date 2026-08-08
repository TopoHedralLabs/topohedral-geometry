//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use topohedral_linalg::VectorOps;

//{{{ crate imports
use crate::common::{zero_vec, Vec3};
//}}}
//{{{ std imports
//}}}
//--------------------------------------------------------------------------------------------------

pub(crate) fn orthogonal_vector(vec: &Vec3) -> Vec3 {
    let v0 = vec[0];
    let s = v0.signum();
    let beta = 1.0 / (1.0 + v0.abs());

    let mut orth_vec = zero_vec();
    orth_vec[0] = -s * vec[1];
    orth_vec[1] = 1.0 - beta * vec[1] * vec[1];
    for i in 2..3 {
        orth_vec[i] = -beta * vec[1] * vec[i];
    }
    orth_vec
}

pub(crate) fn is_unit_vector(vec: &Vec3) -> bool {
    let norm = vec.norm();
    (norm - 1.0).abs() < 1e-12
}
