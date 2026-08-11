//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::{binom_coeff, homog, inv_homog, DynArray, StatArray},
    splines::{self, is_member, PMAX, PMAX_PLUS_ONE, PMAX_PLUS_ONE_SQR},
    utils::WEIGHT_ULPS,
    Surface, Vec3, Vec4,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::ulps_eq;
use topohedral_linalg::SubViewable;
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
        assert!(order_v <= PMAX, "order_v too high");
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

        let num_points_u = knots_u
            .len()
            .checked_sub(order_u + 1)
            .expect("Knots u have wrong size");
        let num_points_v = knots_v
            .len()
            .checked_sub(order_v + 1)
            .expect("Knots v have wrong size");
        assert_eq!(
            num_points_u * num_points_v,
            weights.len(),
            "Knots have wrong size"
        );
        assert!(!points.is_empty(), "At least one point is required");

        let mut points_w = DynArray::<Vec4, 2>::new(Vec4::zeros(), &[num_points_u, num_points_v]);

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
        self.points_w
            .as_slice()
            .iter()
            .skip(1)
            .any(|vec| !ulps_eq!(vec[3], w00, max_ulps = WEIGHT_ULPS))
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

        let mut basis_funs_u = [0.0; PMAX_PLUS_ONE];
        splines::eval(knots_u, u, p, &mut basis_funs_u);
        let mut basis_funs_v = [0.0; PMAX_PLUS_ONE];
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
        let num_derivatives = (m + 1) * (n + 1);
        let mut derivatives =
            StatArray::<Vec3, PMAX_PLUS_ONE_SQR, 2>::new(Vec3::zeros(), &[m + 1, n + 1]);
        self.eval_all_derivatives(
            u,
            v,
            m,
            n,
            &mut derivatives.as_slice_mut()[..num_derivatives],
        );
        derivatives[&[m, n]]
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
        assert!(derivatives.len() == (m + 1) * (n + 1));

        let Self {
            order_u,
            order_v,
            knots_u,
            knots_v,
            points_w,
        } = self;

        let p = *order_u;
        let q = *order_v;

        let (start_u, end_u, num_basis_u) = splines::non_zero_basis(knots_u, u, p);
        let (start_v, end_v, num_basis_v) = splines::non_zero_basis(knots_v, v, q);
        let mut basis_ders_u =
            StatArray::<f64, PMAX_PLUS_ONE_SQR, 2>::new(0.0, &[num_basis_u, m + 1]);
        splines::eval_diff_all(knots_u, u, p, m, basis_ders_u.as_slice_mut());
        let mut basis_ders_v =
            StatArray::<f64, PMAX_PLUS_ONE_SQR, 2>::new(0.0, &[num_basis_v, n + 1]);
        splines::eval_diff_all(knots_v, v, q, n, basis_ders_v.as_slice_mut());

        let mut derivatives_arr =
            StatArray::<Vec4, PMAX_PLUS_ONE_SQR, 2>::new(Vec4::zeros(), &[m + 1, n + 1]);

        for k in 0..m + 1
        // derivatives of u
        {
            for l in 0..n + 1
            // derivatives of v
            {
                for i in start_u..end_u {
                    let i0 = i - start_u;
                    for j in start_v..end_v {
                        let j0 = j - start_v;
                        derivatives_arr[&[k, l]] +=
                            basis_ders_u[&[i0, k]] * basis_ders_v[&[j0, l]] * points_w[&[i, j]]
                    }
                }
            }
        }

        if !self.is_rational() {
            let weight = derivatives_arr[&[0, 0]][3];
            for (derivative, derivative_w) in derivatives.iter_mut().zip(derivatives_arr.as_slice())
            {
                derivative.copy_from(derivative_w.rows_range(0, 2) / weight);
            }
            return;
        }

        let mut surface_derivatives =
            StatArray::<Vec3, PMAX_PLUS_ONE_SQR, 2>::new(Vec3::zeros(), &[m + 1, n + 1]);
        let binom = binom_coeff::<PMAX_PLUS_ONE_SQR>(m.max(n));
        let weight = derivatives_arr[&[0, 0]][3];

        for k in 0..=m
        // derivatives of u
        {
            for l in 0..=n
            // subtract derivatives of v
            {
                let mut derivative = Vec3::zeros();
                derivative.copy_from(derivatives_arr[&[k, l]].rows_range(0, 2));

                for j in 1..=l {
                    derivative -= binom[&[l, j]]
                        * derivatives_arr[&[0, j]][3]
                        * surface_derivatives[&[k, l - j]];
                }

                for i in 1..=k {
                    derivative -= binom[&[k, i]]
                        * derivatives_arr[&[i, 0]][3]
                        * surface_derivatives[&[k - i, l]];

                    let mut mixed_terms = Vec3::zeros();
                    for j in 1..=l {
                        mixed_terms += binom[&[l, j]]
                            * derivatives_arr[&[i, j]][3]
                            * surface_derivatives[&[k - i, l - j]];
                    }
                    derivative -= binom[&[k, i]] * mixed_terms;
                }

                surface_derivatives[&[k, l]].copy_from(derivative / weight);
            }
        }

        derivatives.copy_from_slice(&surface_derivatives.as_slice()[..derivatives.len()]);
    }
}
//}}}
