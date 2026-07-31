//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use super::common::Curve;
use crate::common::{homog, inv_homog, zero_vec, MultDimArray, Vector};
use crate::knot_multiplicity;
use crate::splines::{self, find_span, is_member, PMAX};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use thiserror::Error;
use topohedral_linalg::{SubViewable, TransformOps, VectorOps};
use topohedral_tracing::*;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ collection: constants
pub const BCURVE_DER_MAX: usize = 5;
pub const BCURVE_DER_MAX_SQRD: usize = 25;
//}}}
//{{{ enum: CurveErrro
#[derive(Debug, Error, Clone)]
pub enum CurveError {
    #[error("Construction error {0}")]
    ConstructionError(String),
}
//}}}
//{{{ struct BCurveData
#[derive(Debug, Clone)]
pub struct BCurveData {
    pub p: usize,
    pub knots: Vec<f64>,
    pub points: Vec<Vector>,
    pub weights: Vec<f64>,
}
//}}}
//{{{ struct BCurve
#[derive(Debug, Clone, Default)]
pub struct BCurve {
    dimension: usize,
    order: usize,
    knots: Vec<f64>,
    points_w: Vec<Vector>,
}
//}}}
//{{{ impl: BCurveData
impl BCurveData {
    #[trace_fn]
    fn validate(&self) -> Result<(), CurveError> {
        if self.p > PMAX {
            return Err(CurveError::ConstructionError("Order too high".to_string()));
        }
        if !self.knots.is_sorted() {
            return Err(CurveError::ConstructionError(
                "Knots not sorted".to_string(),
            ));
        }
        if !self.weights.iter().all(|&x| x >= 0.0) {
            return Err(CurveError::ConstructionError(
                "Weights must by positive".to_string(),
            ));
        }
        if self.points.len() != self.weights.len() {
            return Err(CurveError::ConstructionError(
                "Points and weights must have same size".to_string(),
            ));
        }
        if self.knots.len() != (self.weights.len() + self.p + 1) {
            return Err(CurveError::ConstructionError(
                "Knots have wrong size".to_string(),
            ));
        }
        let dim = self.points[0].len();
        if !self.points.iter().all(|x| x.len() == dim) {
            return Err(CurveError::ConstructionError(
                "All points must have same length".to_string(),
            ));
        }
        Ok(())
    }

    fn dim(&self) -> usize {
        self.points[0].len()
    }
}
//}}}
//{{{ ufn: binom_coeff
fn binom_coeff(n: usize) -> MultDimArray<f64, BCURVE_DER_MAX_SQRD, 2> {
    let mut binom = MultDimArray::<f64, BCURVE_DER_MAX_SQRD, 2>::new(0.0, &[n + 1, n + 1]);
    for i in 0..n + 1 {
        binom[&[i, i]] = 1.0;
        binom[&[i, 0]] = 1.0;
    }
    for n2 in 2..n + 1 {
        for k2 in 1..n2 {
            binom[&[n2, k2]] = binom[&[n2 - 1, k2 - 1]] + binom[&[n2 - 1, k2]];
        }
    }
    binom
}
//}}}
//{{{ impl: BCurve
impl BCurve {
    //{{{ fun: new
    #[trace_fn]
    pub fn new(data: BCurveData) -> Result<Self, CurveError> {
        data.validate()?;
        let dim = data.dim();
        let n = data.weights.len();
        let mut points_w = vec![zero_vec(dim + 1); n];
        for (point_w, (point, &weight)) in points_w
            .iter_mut()
            .zip(data.points.iter().zip(&data.weights))
        {
            point_w.copy_from(&homog(point, weight));
        }
        Ok(Self {
            dimension: dim,
            order: data.p,
            knots: data.knots,
            points_w,
        })
    }
    //}}}
    //{{{ fun: insert_knot
    #[trace_fn]
    pub fn insert_knot(
        &mut self,
        u: f64,
        num_insertions: usize,
    ) {
        //{{{ trace
        trace!("Inserting knot {u} {num_insertions} times");
        //}}}
        let Self {
            dimension,
            order,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

        let dim = *dimension;
        let ord = *order;
        let knot_mult = knot_multiplicity(knots_before, u);
        assert!(num_insertions >= 1 && num_insertions <= (ord - knot_mult));
        assert!(is_member(knots_before, u));

        let span = find_span(knots_before, u, ord);
        let num_knots_before = knots_before.len();
        let num_knots_after = num_knots_before + num_insertions;
        let num_points_before = points_w_before.len();
        let num_points_after = num_points_before + num_insertions;

        //{{{ trace
        trace!("knot_mul = {knot_mult} span = {span} m_before = {num_knots_before} m_after = {num_knots_after}");
        //}}}

        let mut knots_after = vec![0.0; num_knots_after];
        let mut points_w_after = vec![zero_vec(dim + 1); num_points_after];
        let mut points_auxiliary = vec![zero_vec(dim + 1); ord + 1];

        knots_after[0..=span].copy_from_slice(&knots_before[0..=span]);
        knots_after[(span + 1)..=(span + num_insertions)].fill(u);
        knots_after[(span + num_insertions + 1)..num_knots_after]
            .copy_from_slice(&knots_before[(span + 1)..num_knots_before]);

        for i in 0..(span - ord + 1) {
            points_w_after[i].copy_from(&points_w_before[i]);
        }
        for i in (span - knot_mult)..num_points_before {
            points_w_after[i + num_insertions].copy_from(&points_w_before[i]);
        }
        for i in 0..(ord - knot_mult + 1) {
            points_auxiliary[i].copy_from(&points_w_before[span - ord + i]);
        }

        let mut p1 = zero_vec(dim + 1);
        let mut p2 = zero_vec(dim + 1);

        let mut l = span;
        for j in 1..(num_insertions + 1) {
            // insert the knot num_insertions times
            l = span - ord + j;
            for i in 0..(ord - j - knot_mult + 1) {
                let numer = u - knots_before[l + i];
                let denom = knots_before[i + span + 1] - knots_before[l + i];
                let alpha = numer / denom;
                p1.copy_from(&points_auxiliary[i + 1]);
                p2.copy_from(&points_auxiliary[i]);
                points_auxiliary[i].copy_from(alpha * &p1 + (1.0 - alpha) * &p2);
            }
            points_w_after[l].copy_from(&points_auxiliary[0]);
            points_w_after[span + num_insertions - j - knot_mult]
                .copy_from(&points_auxiliary[ord - j - knot_mult]);
        }
        for i in (l + 1)..(span - knot_mult) {
            points_w_after[i].copy_from(&points_auxiliary[i - l]);
        }
        *knots_before = knots_after;
        *points_w_before = points_w_after;
    }
    //}}}
    //{{{ fun: with_inserted_knot
    #[trace_fn]
    pub fn with_inserted_knot(
        &self,
        u: f64,
        r: usize,
    ) -> Self {
        let mut new_curve = self.clone();
        new_curve.insert_knot(u, r);
        new_curve
    }
    //}}}
    //{{{ fun: refine_knots
    #[trace_fn]
    pub fn refine_knots(
        &mut self,
        u: &[f64],
    ) {
        debug_assert!(!u.is_empty());
        debug_assert!(u.is_sorted());
        debug_assert!(is_member(&self.knots, u[0]));
        debug_assert!(is_member(&self.knots, u[u.len() - 1]));

        //{{{ trace
        trace!("Inserting knots {u:?}");
        //}}}

        let Self {
            dimension,
            order,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

        let dim = *dimension;
        let ord = *order;

        let n_new_knots = u.len();
        let n_knots_before = knots_before.len();
        let n_points_before = points_w_before.len();

        let n_knots_after = n_knots_before + n_new_knots + 1;
        let n_points_after = n_points_before + n_new_knots + 1;
        let mut knots_after = vec![0.0; n_knots_after];
        let mut points_w_after = vec![zero_vec(dim + 1); n_points_after];

        let span_start = find_span(knots_before, *u.first().unwrap(), ord);
        let span_end = find_span(knots_before, *u.last().unwrap(), ord) + 1;
        //{{{ trace
        trace!("span_start = {span_start} span_end = {span_end}");
        //}}}

        // copy knots before affected_range
        knots_after[0..=span_start].copy_from_slice(&knots_before[0..=span_start]);
        // copy knots after affected range
        knots_after[(span_end + n_new_knots + ord + 1)..n_knots_after]
            .copy_from_slice(&knots_before[(span_end + ord + 1)..n_knots_before]);

        // copy points before affected range
        for j in 0..=(span_start - ord) {
            points_w_after[j].copy_from(&points_w_before[j]);
        }
        // copy points after affected range
        for j in (span_end - 1)..n_points_after {
            points_w_after[j + n_new_knots + 1].copy_from(&points_w_before[j]);
        }

        let mut low_lim = span_end + ord - 1;
        let mut high_lim = span_end + ord + n_new_knots;

        for j in (0..=n_new_knots).rev() {
            //{{{ trace
            trace_scope!("loop1");
            trace!("...................... j = {j}");
            //}}}
            while u[j] <= knots_before[low_lim] && low_lim > span_start {
                let idx_after = high_lim - ord - 1;
                let idx_before = low_lim - ord - 1;
                //{{{ trace
                trace_scope!("loop2");
                trace!(".... low_lim = {low_lim} high_lim = {high_lim}");
                trace!("Setting new_point {idx_after} from {idx_before}");
                //}}}
                points_w_after[idx_after].copy_from(&points_w_before[idx_before]);
                knots_after[high_lim] = knots_before[low_lim];
                low_lim += 1;
                high_lim -= 1;
            }
            points_w_after[high_lim - ord - 1].copy_from(&points_w_before[low_lim]);
            for l in 1..=ord {
                //{{{ trace
                trace_scope!("loop3");
                trace!("...................... j = {j}");
                //}}}
                let ind = high_lim - ord + l;
                let mut alpha = knots_after[high_lim + l] - u[j];
                if alpha.abs() == 0.0 {
                    let tmp = points_w_after[ind].clone();
                    points_w_after[ind - 1].copy_from(&tmp);
                } else {
                    alpha /= (knots_after[high_lim + l] - knots_before[low_lim - ord + l]);
                    let tmp1 = points_w_after[ind - 1].clone();
                    let tmp2 = points_w_after[ind].clone();
                    points_w_after[ind - 1].copy_from(alpha * &tmp1 + (1.0 - alpha) * tmp2);
                }
            }
            knots_after[high_lim] = knots_before[low_lim];
            high_lim -= 1;
        }
        *knots_before = knots_after;
        *points_w_before = points_w_after;
    }
    //}}}
    //{{{ fun: remove_knot
    pub fn remove_knot(
        &mut self,
        u: f64,
    ) {
    }
    //}}}
    //{{{ fun: elevate_degree
    pub fn elevate_degree(
        &mut self,
        u: f64,
    ) {
    }
    //}}}
    //{{{ fun: reduce_degree
    pub fn reduce_degree(
        &mut self,
        u: f64,
    ) {
    }
    //}}}
}
//}}}
//{{{ impl Curve for Bcurve
impl Curve for BCurve {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vector {
        assert!(is_member(&self.knots, u));

        let mut point_w = zero_vec(self.dimension + 1);
        let (start, end, _) = splines::non_zero_basis(&self.knots, u, self.order);

        let mut bases_funs = [0.0; PMAX];
        splines::eval(&self.knots, u, self.order, &mut bases_funs);
        for i in start..end {
            point_w += bases_funs[i - start] * &self.points_w[i];
        }
        inv_homog(&point_w)
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vector {
        assert!(m < BCURVE_DER_MAX, "derivative order exceeds maximum");

        let mut derivatives: [Vector; BCURVE_DER_MAX] =
            std::array::from_fn(|_| zero_vec(self.dimension));
        self.eval_all_derivatives(u, m, &mut derivatives[..=m]);
        derivatives[m].clone()
    }

    fn eval_all_derivatives(
        &self,
        u: f64,
        derivative_order: usize,
        derivatives: &mut [Vector],
    ) {
        assert!(is_member(&self.knots, u));
        assert!(derivatives.len() > derivative_order);
        assert!(
            derivative_order < BCURVE_DER_MAX,
            "derivative order exceeds maximum"
        );

        if derivative_order == 0 {
            derivatives[0].copy_from(self.eval_point(u));
        } else {
            let p = self.order;
            let knots = &self.knots;
            let (start, _, num_basis) = splines::non_zero_basis(knots, u, p);
            let num_derivatives = derivative_order + 1;
            let dim = self.dimension;
            let mut curve_derivatives_w: [Vector; BCURVE_DER_MAX] =
                std::array::from_fn(|_| zero_vec(dim + 1));
            let mut basis_ders_arr = MultDimArray::<f64, BCURVE_DER_MAX_SQRD, 2>::new(
                0.0,
                &[num_basis, num_derivatives],
            );
            splines::eval_diff_all(knots, u, p, derivative_order, basis_ders_arr.as_slice_mut());

            for m in 0..num_derivatives {
                for j in 0..num_basis {
                    let dnj_dxm = basis_ders_arr[&[j, m]];
                    let point_j = &self.points_w[start + j];
                    curve_derivatives_w[m] += dnj_dxm * point_j
                }
            }

            let binom = binom_coeff(derivative_order);
            let mut ders_loc: [Vector; BCURVE_DER_MAX] = std::array::from_fn(|_| zero_vec(dim));
            let w0 = curve_derivatives_w[0][dim];
            let mut v = zero_vec(dim);

            for m in 0..num_derivatives {
                v.fill(0.0);
                v.copy_from(&curve_derivatives_w[m].rows_range(0, dim - 1));

                for j in 1..m + 1 {
                    let wj = curve_derivatives_w[j][dim];
                    let bmj = binom[&[m, j]];
                    v -= bmj * wj * &ders_loc[m - j];
                }
                ders_loc[m].copy_from(&v / w0);
                derivatives[m].copy_from(&ders_loc[m]);
            }
        }
    }
}
//}}}
