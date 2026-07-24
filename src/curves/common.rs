//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

use core::panic;
use std::todo;

use topohedral_linalg::VectorOps;

//{{{ crate imports
use crate::common::{vec_from_values, zero_vec, Vector};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

pub trait Curve {
    fn dimension(&self) -> usize;

    fn eval_point(
        &self,
        u: f64,
    ) -> Vector;

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vector;

    fn eval_all_derivatives(
        &self,
        u: f64,
        m: usize,
        derivatives: &mut [Vector],
    );

    fn eval_tangent(
        &self,
        u: f64,
        normalize: bool,
    ) -> Vector {
        let tangent = self.eval_derivative(u, 1);
        if normalize {
            return tangent.normalize();
        }
        tangent
    }

    fn eval_normal(
        &self,
        u: f64,
        normalize: bool,
    ) -> Vector {
        match self.dimension() {
            2 => {
                let tangent = self.eval_tangent(u, normalize);
                let normal = vec_from_values(&[-tangent[0], tangent[1]]);
                normal
            }
            3 => {
                let mut curve_dirivatives: [Vector; 3] = std::array::from_fn(|_| zero_vec(3));
                self.eval_all_derivatives(u, 2, &mut curve_dirivatives);

                let velocity = &curve_dirivatives[1];
                let acceleration = &curve_dirivatives[2];

                let b = acceleration.cross(velocity);
                let b_norm = b.norm();
                let mut normal = velocity.cross(&b) * b_norm;

                if normalize {
                    return normal.normalize();
                }
                normal
            }
            _ => {
                panic!("Unsupported dimensionality")
            }
        }
    }

    fn eval_binormal(
        &self,
        u: f64,
        normalize: bool,
    ) -> Vector {
        match self.dimension() {
            2 => zero_vec(2),
            3 => {
                let tangent = self.eval_tangent(u, false);
                let normal = self.eval_normal(u, false);
                let binormal = tangent.cross(&normal);
                if normalize {
                    return binormal.normalize();
                }
                binormal
            }
            _ => {
                panic!("Unsupported dimensionality")
            }
        }
    }
}
