//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::common::{zero_vec, Vector};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub(crate) fn orthogonal_vector(vec: &Vector) -> Vector {
    let v0 = vec[0];
    let s = v0.signum();
    let beta = 1.0 / (1.0 + v0.abs());

    let mut orth_vec = zero_vec(vec.len());
    orth_vec[0] = -s * vec[1];
    orth_vec[1] = 1.0 - beta * vec[1] * vec[1];
    for i in 2..vec.len() {
        orth_vec[i] = -beta * vec[1] * vec[i];
    }
    orth_vec
}
