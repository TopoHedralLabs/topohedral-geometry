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
//{{{ dep imports
//}}}
//--------------------------------------------------------------------------------------------------

pub trait Curve {
    fn dimension(&self) -> usize {
        3
    }

    fn eval_point(
        &self,
        _u: f64,
    ) -> Vec3;

    fn eval_derivative(
        &self,
        _u: f64,
        m: usize,
    ) -> Vec3;

    fn eval_all_derivatives(
        &self,
        u: f64,
        m: usize,
        derivatives: &mut [Vec3],
    ) {
        assert_eq!(derivatives.len(), m + 1);
        for (i, derivative) in derivatives.iter_mut().enumerate().take(m + 1) {
            derivative.copy_from(&self.eval_derivative(u, i));
        }
    }

    fn eval_tangent(
        &self,
        u: f64,
        normalize: bool,
    ) -> Vec3 {
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
    ) -> Vec3 {
        let mut curve_dirivatives: [Vec3; 3] = std::array::from_fn(|_| zero_vec());
        self.eval_all_derivatives(u, 2, &mut curve_dirivatives);

        let velocity = &curve_dirivatives[1];
        let acceleration = &curve_dirivatives[2];

        let b = acceleration.cross(velocity);
        let b_norm = b.norm();
        let normal = velocity.cross(&b) * b_norm;

        if normalize {
            return normal.normalize();
        }
        normal
    }

    fn eval_binormal(
        &self,
        u: f64,
        normalize: bool,
    ) -> Vec3 {
        let tangent = self.eval_tangent(u, false);
        let normal = self.eval_normal(u, false);
        let binormal = tangent.cross(&normal);
        if normalize {
            return binormal.normalize();
        }
        binormal
    }

    fn eval_curvature(
        &self,
        _u: f64,
    ) -> f64 {
        0.0
    }

    fn eval_torsion(
        &self,
        _u: f64,
    ) -> f64 {
        0.0
    }

    fn eval_arclen(
        &self,
        _ustart: f64,
        _uend: f64,
    ) -> f64 {
        0.0
    }
}
