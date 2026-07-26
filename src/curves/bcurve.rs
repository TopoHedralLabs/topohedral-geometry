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
use topohedral_linalg::{tan, SubViewable, TransformOps, VectorOps};
use topohedral_tracing::*;
//}}}
//--------------------------------------------------------------------------------------------------

pub const BCURVE_DER_MAX: usize = 5;
pub const BCURVE_DER_MAX_SQRD: usize = 25;

#[derive(Debug, Error, Clone)]
pub enum CurveError {
    #[error("Construction error {0}")]
    ConstructionError(String),
}

#[derive(Debug, Clone)]
pub struct BCurveData {
    pub p: usize,
    pub knots: Vec<f64>,
    pub points: Vec<Vector>,
    pub weights: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct BCurve {
    dimension: usize,
    order: usize,
    knots: Vec<f64>,
    points_w: Vec<Vector>,
}

impl BCurveData {
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

impl BCurve {
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

    #[trace_fn]
    pub fn insert_knot(
        &mut self,
        u: f64,
        r: usize,
    ) {
        //{{{ trace
        trace!("Inserting knot {u} {r} times into {self:?}");
        //}}}
        let Self {
            dimension: dim,
            order: p,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

        let knot_mul = knot_multiplicity(knots_before, u);
        assert!(r >= 1 && r <= (*p - knot_mul));

        let span = find_span(knots_before, u, *p);
        let m_before = knots_before.len();
        let m_after = m_before + r;
        let n_before = points_w_before.len();
        let n_after = n_before + r;

        //{{{ trace
        trace!("knot_mul = {knot_mul} span = {span} m_before = {m_before} m_after = {m_after}");
        //}}}

        let mut knots_after = vec![0.0; m_after];
        let mut points_w_after = vec![zero_vec(*dim + 1); n_after];
        let mut point_auxiliary = vec![zero_vec(*dim + 1); *p + 1];

        knots_after[0..=span].copy_from_slice(&knots_before[0..=span]);
        knots_after[(span + 1)..=(span + r)].fill(u);
        knots_after[(span + r + 1)..m_after].copy_from_slice(&knots_before[(span + 1)..m_before]);

        for i in 0..(span - *p + 1) {
            points_w_after[i].copy_from(&points_w_before[i]);
        }
        for i in (span - knot_mul)..n_before {
            points_w_after[i + r].copy_from(&points_w_before[i]);
        }
        for i in 0..(*p - knot_mul + 1) {
            point_auxiliary[i].copy_from(&points_w_after[span - *p + i]);
        }
    }

    pub fn with_inserted_knot(
        &self,
        u: f64,
        r: usize,
    ) -> Self {
        let mut new_curve = self.clone();
        new_curve.insert_knot(u, r);
        new_curve
    }

    pub fn insert_knots(
        &mut self,
        u: &[f64],
    ) {
    }

    pub fn remove_knot(
        &mut self,
        u: f64,
    ) {
    }

    pub fn elevate_degree(
        &mut self,
        u: f64,
    ) {
    }

    pub fn reduce_degree(
        &mut self,
        u: f64,
    ) {
    }
}

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
