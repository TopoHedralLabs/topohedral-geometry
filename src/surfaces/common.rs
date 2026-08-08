//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{common::IndexHelper, Vec3};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub trait Surface {
    fn dimension(&self) -> usize {
        3
    }

    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3;

    fn eval_derivative(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
    ) -> Vec3;

    fn eval_all_derivatives(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
        derivatives: &mut [Vec3],
    ) {
        assert!(derivatives.len() == (m + 1) * (n + 1));
        let idx_helper = IndexHelper::<2>::new(&[m + 1, n + 1]);
        for i in 0..=m {
            for j in 0..=n {
                let idx = idx_helper.lin_index(&[i, j]);
                derivatives[idx].copy_from(self.eval_derivative(u, v, i, j));
            }
        }
    }

    fn eval_derivative_u(
        &self,
        u: f64,
        v: f64,
        m: usize,
    ) -> Vec3 {
        self.eval_derivative(u, v, m, 0)
    }

    fn eval_derivative_v(
        &self,
        u: f64,
        v: f64,
        n: usize,
    ) -> Vec3 {
        self.eval_derivative(u, v, 0, n)
    }

    fn eval_normal(
        &self,
        u: f64,
        v: f64,
        normalize: bool,
    ) -> Vec3 {
        let deriv_u = self.eval_derivative_u(u, v, 1);
        let deriv_v = self.eval_derivative_v(u, v, 1);
        let normal = deriv_u.cross(&deriv_v);
        if normalize {
            return normal.normalize();
        }
        normal
    }
}
