//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use super::common::Curve;
use crate::common::{binom_coeff, homog, inv_homog, zero_vec, zero_vec4, MultDimArray, Vec3, Vec4};
use crate::knot_multiplicity;
use crate::splines::{self, find_span, is_member, knot_eq, PMAX, PMAX_SQR};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::TransformOps;
use topohedral_tracing::*;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ collection: constants
pub const BCURVE_DER_MAX: usize = 5;
pub const BCURVE_DER_MAX_SQRD: usize = 25;
//}}}
//{{{ struct BCurve
#[derive(Debug, Clone, Default)]
pub struct BCurve {
    order: usize,
    knots: Vec<f64>,
    points_w: Vec<Vec4>,
}
//}}}
//{{{ impl: BCurve
impl BCurve {
    //{{{ fun: new
    #[trace_fn]
    pub fn new(
        order: usize,
        knots: Vec<f64>,
        points: Vec<Vec3>,
        weights: Vec<f64>,
    ) -> Self {
        assert!(order <= PMAX, "Order too high");
        assert!(knots.is_sorted(), "Knots not sorted");
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
            knots.len(),
            weights.len() + order + 1,
            "Knots have wrong size"
        );
        assert!(!points.is_empty(), "At least one point is required");

        let n = weights.len();
        let mut points_w = vec![zero_vec4(); n];
        for (point_w, (point, &weight)) in points_w.iter_mut().zip(points.iter().zip(&weights)) {
            point_w.copy_from(&homog(point, weight));
        }
        Self {
            order,
            knots,
            points_w,
        }
    }
    //}}}
    //{{{ fun: dimension
    pub fn dimension(&self) -> usize {
        3
    }
    //}}}
    //{{{ fun: order
    pub fn order(&self) -> usize {
        self.order
    }
    //}}}
    //{{{ fun: knots
    pub fn knots(&self) -> &Vec<f64> {
        &self.knots
    }
    //}}}
    //{{{ fun: points_w
    pub fn points_w(&self) -> &Vec<Vec4> {
        &self.points_w
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
            order,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

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
        let mut points_w_after = vec![zero_vec4(); num_points_after];
        let mut points_auxiliary = vec![zero_vec4(); ord + 1];

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

        let mut p1 = zero_vec4();
        let mut p2 = zero_vec4();

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
                points_auxiliary[i].copy_from(alpha * p1 + (1.0 - alpha) * p2);
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
            order,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

        let ord = *order;

        let n_new_knots = u.len();
        let n_knots_before = knots_before.len();
        let n_points_before = points_w_before.len();

        let n_knots_after = n_knots_before + n_new_knots;
        let n_points_after = n_points_before + n_new_knots;
        let mut knots_after = vec![0.0; n_knots_after];
        let mut points_w_after = vec![zero_vec4(); n_points_after];

        let span_start = find_span(knots_before, *u.first().unwrap(), ord);
        let span_end = find_span(knots_before, *u.last().unwrap(), ord) + 1;
        //{{{ trace
        trace!("span_start = {span_start} span_end = {span_end}");
        //}}}

        // copy knots before affected_range
        let r1 = (0, span_start);
        //{{{ trace
        trace!(
            "Copying knots before affected range in ({}, {})",
            r1.0,
            r1.1
        );
        //}}}
        knots_after[r1.0..=r1.1].copy_from_slice(&knots_before[r1.0..=r1.1]);

        // copy knots after affected range
        let r2 = (span_end + ord, n_knots_before);
        //{{{ trace
        trace!("Copying knots after affected range in ({}, {})", r2.0, r2.1);
        //}}}
        knots_after[(r2.0 + n_new_knots)..n_knots_after].copy_from_slice(&knots_before[r2.0..r2.1]);

        // copy points before affected range
        let r3 = (0, span_start - ord);
        //{{{ trace
        trace!("Copying points before affected range ({} {})", r3.0, r3.1);
        //}}}
        for j in r3.0..=r3.1 {
            points_w_after[j].copy_from(&points_w_before[j]);
        }
        // copy points after affected range
        let r4 = (span_end - 1, n_points_before);
        //{{{ trace
        trace!("Copying points after affected range ({} {})", r4.0, r4.1);
        //}}}
        for j in r4.0..r4.1 {
            points_w_after[j + n_new_knots].copy_from(&points_w_before[j]);
        }

        let mut low_lim = span_end + ord - 1;
        let mut high_lim = span_end + ord + n_new_knots - 1;

        for j in (0..n_new_knots).rev() {
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
                low_lim -= 1;
                high_lim -= 1;
            }
            let tmp = points_w_after[high_lim - ord];
            points_w_after[high_lim - ord - 1].copy_from(&tmp);
            for l in 1..=ord {
                //{{{ trace
                trace_scope!("loop3");
                trace!("...................... j = {j}");
                //}}}
                let ind = high_lim - ord + l;
                let mut alpha = knots_after[high_lim + l] - u[j];
                if alpha.abs() == 0.0 {
                    let tmp = points_w_after[ind];
                    points_w_after[ind - 1].copy_from(&tmp);
                } else {
                    alpha /= knots_after[high_lim + l] - knots_before[low_lim - ord + l];
                    let tmp1 = points_w_after[ind - 1];
                    let tmp2 = points_w_after[ind];
                    points_w_after[ind - 1].copy_from(alpha * tmp1 + (1.0 - alpha) * tmp2);
                }
            }
            knots_after[high_lim] = u[j];
            high_lim -= 1;
        }
        *knots_before = knots_after;
        *points_w_before = points_w_after;
    }
    //}}}
    //{{{ fun: with_refined_knots
    #[trace_fn]
    pub fn with_refined_knots(
        &self,
        u: &[f64],
    ) -> BCurve {
        let mut new_curve = self.clone();
        new_curve.refine_knots(u);
        new_curve
    }
    //}}}
    //{{{ fun: elevate_degree
    pub fn elevate_degree(
        &mut self,
        q: usize,
    ) {
        if q == 0 {
            return;
        }

        let Self {
            order,
            knots: knots_before,
            points_w: points_w_before,
        } = self;

        let ord_before = *order;
        let ord_after = ord_before + q;
        assert!(ord_after <= PMAX, "Elevated order exceeds maximum");

        let n_knots_before = knots_before.len();
        // `last_knot` is the `m` index from Algorithm A5.9 in The NURBS Book.
        let last_knot = n_knots_before - 1;

        let comb = binom_coeff::<PMAX_SQR>(ord_after);

        let mut bezalfs =
            MultDimArray::<f64, PMAX_SQR, 2>::new(0.0, &[ord_after + 1, ord_before + 1]);

        bezalfs[&[0, 0]] = 1.0;
        bezalfs[&[ord_after, ord_before]] = 1.0;

        for i in 1..((ord_after / 2) + 1) {
            let inv = 1.0 / comb[&[ord_after, i]];
            let mpi = ord_before.min(i);
            let start = (0i64.max(i as i64 - q as i64)) as usize;
            for j in start..=mpi {
                bezalfs[&[i, j]] = inv * comb[&[ord_before, j]] * comb[&[q, i - j]];
            }
        }
        for i in ((ord_after / 2) + 1)..ord_after {
            let mpi = ord_before.min(i);
            let start = (0i64.max(i as i64 - q as i64)) as usize;
            for j in start..=mpi {
                bezalfs[&[i, j]] = bezalfs[&[ord_after - i, ord_before - j]]
            }
        }

        let n_knots_after_ub = (n_knots_before + 1) * (q + 1) * ord_after + 1;
        let n_points_after_ub = (n_knots_before + 1) * (q + 1);
        let mut knots_after = vec![0.0; n_knots_after_ub];
        let mut point_w_after = vec![zero_vec4(); n_points_after_ub];

        let mut mh = ord_after;
        let mut kind = ord_after + 1;
        let mut r = -1isize;
        let mut a = ord_before;
        let mut b = ord_before + 1;
        let mut cind = 1;
        let mut ua = *knots_before.first().unwrap();
        point_w_after[0].copy_from(&points_w_before[0]);
        knots_after[0..(ord_after + 1)].fill(ua);

        let mut bpts: Vec<Vec4> = points_w_before[0..(ord_before + 1)].to_vec();

        let mut ebpts = vec![zero_vec4(); ord_after + 1];
        let mut next_bpts = vec![zero_vec4(); ord_before.saturating_sub(1)];

        while b < last_knot {
            let i = b;
            while b < last_knot && knot_eq(knots_before[b], knots_before[b + 1]) {
                b += 1;
            }

            let mul = b - i + 1;
            mh += mul + q;
            let ub = knots_before[b];
            let oldr = r;
            r = ord_before as isize - mul as isize;
            let lbz = if oldr > 0 {
                ((oldr + 2) / 2) as usize
            } else {
                1
            };
            let rbz = if r > 0 {
                (ord_after as isize - ((r + 1) / 2)) as usize
            } else {
                ord_after
            };

            if r > 0 {
                let numer = ub - ua;
                let mut alphas = vec![0.0; ord_before.saturating_sub(1)];

                for k in (mul + 1..=ord_before).rev() {
                    alphas[k - mul - 1] = numer / (knots_before[a + k] - ua);
                }
                for j in 1..=(r as usize) {
                    let save = r as usize - j;
                    let s = mul + j;
                    for k in (s..=ord_before).rev() {
                        let alpha = alphas[k - s];
                        let old_k = bpts[k];
                        let old_k_minus_one = bpts[k - 1];
                        let tmp = alpha * old_k + (1.0 - alpha) * old_k_minus_one;
                        bpts[k].copy_from(&tmp);
                    }
                    next_bpts[save].copy_from(&bpts[ord_before]);
                }
            }

            for i in lbz..=ord_after {
                ebpts[i].fill(0.0);
                let mpi = ord_before.min(i);
                let start = i.saturating_sub(q);
                for j in start..=mpi {
                    let tmp = bezalfs[&[i, j]] * bpts[j];
                    ebpts[i] += tmp;
                }
            }

            if oldr > 1 {
                let mut first = kind - 2;
                let den = ub - ua;
                let bet = (ub - knots_after[kind - 1]) / den;

                for (last, tr) in (kind..).zip(1..(oldr as usize)) {
                    let mut i = first;
                    let mut j = last;
                    let mut kj = j - kind + 1;
                    while j - i > tr {
                        if i < cind {
                            let alf = (ub - knots_after[i]) / (ua - knots_after[i]);
                            let old_i = point_w_after[i];
                            let old_i_minus_one = point_w_after[i - 1];
                            let tmp = alf * old_i + (1.0 - alf) * old_i_minus_one;
                            point_w_after[i].copy_from(&tmp);
                        }
                        if j >= lbz {
                            let threshold = kind as isize - ord_after as isize + oldr;
                            if (j - tr) as isize <= threshold {
                                let gam = (ub - knots_after[j - tr]) / den;
                                let old_kj = ebpts[kj];
                                let old_kj_plus_one = ebpts[kj + 1];
                                let tmp = gam * old_kj + (1.0 - gam) * old_kj_plus_one;
                                ebpts[kj].copy_from(&tmp);
                            } else {
                                let old_kj = ebpts[kj];
                                let old_kj_plus_one = ebpts[kj + 1];
                                let tmp = bet * old_kj + (1.0 - bet) * old_kj_plus_one;
                                ebpts[kj].copy_from(&tmp);
                            }
                        }
                        i += 1;
                        j -= 1;
                        kj -= 1;
                    }
                    first -= 1;
                }
            }

            if a != ord_before {
                for _ in 0..(ord_after - oldr as usize) {
                    knots_after[kind] = ua;
                    kind += 1;
                }
            }
            for point in ebpts.iter().take(rbz + 1).skip(lbz) {
                point_w_after[cind].copy_from(point);
                cind += 1;
            }

            if b < last_knot {
                for j in 0..(r as usize) {
                    bpts[j].copy_from(&next_bpts[j]);
                }
                for j in (r as usize)..=ord_before {
                    bpts[j].copy_from(&points_w_before[b - ord_before + j]);
                }
                a = b;
                b += 1;
                ua = ub;
            } else {
                for i in 0..=ord_after {
                    knots_after[kind + i] = ub;
                }
            }
        }

        let nh = mh - ord_after - 1;
        *order = ord_after;
        *knots_before = knots_after[..=mh].to_vec();
        *points_w_before = point_w_after[..=nh].to_vec();
    }
    //}}}
}
//}}}
//{{{ impl Curve for Bcurve
impl Curve for BCurve {
    fn dimension(&self) -> usize {
        3
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        assert!(is_member(&self.knots, u));

        let mut point_w = zero_vec4();
        let (start, end, _) = splines::non_zero_basis(&self.knots, u, self.order);

        let mut bases_funs = [0.0; PMAX];
        splines::eval(&self.knots, u, self.order, &mut bases_funs);
        for i in start..end {
            point_w += bases_funs[i - start] * self.points_w[i];
        }
        inv_homog(&point_w)
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        assert!(m < BCURVE_DER_MAX, "derivative order exceeds maximum");

        let mut derivatives: [Vec3; BCURVE_DER_MAX] = std::array::from_fn(|_| zero_vec());
        self.eval_all_derivatives(u, m, &mut derivatives[..=m]);
        derivatives[m]
    }

    fn eval_all_derivatives(
        &self,
        u: f64,
        derivative_order: usize,
        derivatives: &mut [Vec3],
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
            let mut curve_derivatives_w: [Vec4; BCURVE_DER_MAX] =
                std::array::from_fn(|_| zero_vec4());
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

            let binom = binom_coeff::<BCURVE_DER_MAX_SQRD>(derivative_order);
            let mut ders_loc: [Vec3; BCURVE_DER_MAX] = std::array::from_fn(|_| zero_vec());
            let w0 = curve_derivatives_w[0][3];
            let mut v = zero_vec();

            for m in 0..num_derivatives {
                v.fill(0.0);
                v[0] = curve_derivatives_w[m][0];
                v[1] = curve_derivatives_w[m][1];
                v[2] = curve_derivatives_w[m][2];

                for j in 1..m + 1 {
                    let wj = curve_derivatives_w[j][3];
                    let bmj = binom[&[m, j]];
                    v -= bmj * wj * ders_loc[m - j];
                }
                ders_loc[m].copy_from(v / w0);
                derivatives[m].copy_from(&ders_loc[m]);
            }
        }
    }
}
//}}}
