//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::{homog, inv_homog, DynArray, StatArray},
    splines::{self, is_member, PMAX, PMAX_PLUS_ONE, PMAX_PLUS_ONE_SQR, PMAX_SQR},
    utils::WEIGHT_ULPS,
    Surface, Vec3, Vec4,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::ulps_eq;
use topohedral_linalg::TransformOps;
use topohedral_tracing::trace_fn;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ struct: BSurface
pub struct BSurface {
    order_u: usize,
    order_v: usize,
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
    points_w: DynArray<Vec4, 2>,
}
//}}}
//{{{ impl: BSurface
impl BSurface {
    //{{{ fun: new
    #[trace_fn]
    pub fn new(
        order_u: usize,
        order_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        points: Vec<Vec3>,
        weights: Vec<f64>,
    ) -> Self {
        assert!(order_u <= PMAX, "order_u too high");
        assert!(order_u <= PMAX, "order v too high");
        assert!(knots_u.is_sorted(), "knots u not sorted");
        assert!(knots_v.is_sorted(), "knots v not sorted");
        assert!(
            weights.iter().all(|&weight| weight >= 0.0),
            "Weights must be positive"
        );
        assert_eq!(
            points.len(),
            weights.len(),
            "Points and weights must have same size"
        );
        assert_eq!(
            (knots_u.len() + order_u + 1) * (knots_v.len() + order_v + 1),
            weights.len(),
            "Knots have wrong size"
        );

        let num_knots_u = knots_u.len();
        let num_knots_v = knots_v.len();
        let num_points_u = num_knots_u - (order_u + 1);
        let num_points_v = num_knots_v - (order_v + 1);

        let mut points_w = DynArray::<Vec4, 2>::new(Vec4::zeros(), &[num_knots_u, num_knots_v]);

        let mut idx = 0;
        for i in 0..num_points_u {
            for j in 0..num_points_v {
                let point_ij = points[idx];
                let weight_ij = weights[idx];
                let point_w = homog(&point_ij, weight_ij);
                points_w[&[i, j]].copy_from(point_w);
                idx += 1;
            }
        }
        Self {
            order_u,
            order_v,
            knots_u,
            knots_v,
            points_w,
        }
    }
    //}}}
    //{{{ fun: is_rational
    pub fn is_rational(&self) -> bool {
        let w00 = self.points_w[&[0, 0]][3];
        for vec in self.points_w.as_slice().iter().take(1) {
            if !ulps_eq!(vec[3], w00, max_ulps = WEIGHT_ULPS) {
                return false;
            }
        }
        true
    }
    //}}}
}
//}}}
//{{{ impl: Surface for BSurface
impl Surface for BSurface {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        assert!(is_member(&self.knots_u, u));
        assert!(is_member(&self.knots_v, v));

        let Self {
            order_u,
            order_v,
            knots_u,
            knots_v,
            points_w,
        } = self;

        let p = *order_u;
        let q = *order_v;

        let mut point_w = Vec4::zeros();
        let (start_u, end_u, _) = splines::non_zero_basis(knots_u, u, p);
        let (start_v, end_v, _) = splines::non_zero_basis(knots_v, v, q);

        let mut basis_funs_u = [0.0; PMAX];
        splines::eval(knots_u, u, p, &mut basis_funs_u);
        let mut basis_funs_v = [0.0; PMAX];
        splines::eval(knots_v, v, q, &mut basis_funs_v);

        for i in start_u..end_u {
            let bas_i = basis_funs_u[i - start_u];
            for j in start_v..end_v {
                let bas_j = basis_funs_v[j - start_v];
                point_w += bas_i * bas_j * points_w[&[i, j]];
            }
        }
        inv_homog(&point_w)
    }

    fn eval_derivative(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
    ) -> Vec3 {
        let mut derivatives = StatArray::<Vec3, PMAX_SQR, 2>::new(Vec3::zeros(), &[m, n]);
        self.eval_all_derivatives(u, v, m, n, derivatives.as_slice_mut());
        return derivatives[&[m, n]];
    }

    fn eval_all_derivatives(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
        derivatives: &mut [Vec3],
    ) {
        assert!(m <= self.order_u);
        assert!(n <= self.order_v);
        assert!(is_member(&self.knots_u, u));
        assert!(is_member(&self.knots_v, v));

        let Self {
            order_u,
            order_v,
            knots_u,
            knots_v,
            points_w,
        } = self;

        let p = *order_u;
        let q = *order_v;

        let (start_u, end_u, _) = splines::non_zero_basis(knots_u, u, p);
        let (start_v, end_v, _) = splines::non_zero_basis(knots_v, v, q);
        let mut basis_ders_u = StatArray::<f64, PMAX_PLUS_ONE_SQR, 2>::new(0.0, &[m + 1, p + 1]);
        splines::eval_diff_all(knots_u, u, p, m, basis_ders_u.as_slice_mut());
        let mut basis_ders_v = StatArray::<f64, PMAX_PLUS_ONE_SQR, 2>::new(0.0, &[m + 1, p + 1]);
        splines::eval_diff_all(knots_v, v, q, n, basis_ders_v.as_slice_mut());

        let mut vv = Vec4::zeros();
        let mut derivatives =
            StatArray::<Vec4, PMAX_PLUS_ONE_SQR, 2>::new(Vec4::zeros(), &[m + 1, n + 1]);
        let mut temp = [Vec4::zeros(); PMAX_PLUS_ONE];

        todo!()
    }
}
//}}}
