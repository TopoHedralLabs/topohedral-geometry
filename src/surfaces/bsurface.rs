//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::{binom_coeff, homog, inv_homog, DynArray, StatArray},
    curves::BCurveVtkOptions,
    splines::{self, is_member, non_zero_basis, PMAX, PMAX_PLUS_ONE, PMAX_PLUS_ONE_SQR},
    surfaces::common::add_surface_to_vtk,
    utils::WEIGHT_ULPS,
    vtk::{ToVtk, VtkUgridFile},
    BCurve, Surface, Vec3, Vec4,
};

//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::ulps_eq;
use topohedral_linalg::SubViewable;
use topohedral_tracing::trace_fn;
use vtkio::model::CellType;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ struct: BSurface
#[derive(Debug, Clone, Default)]
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
    //{{{ fun: order_u
    pub fn order_u(&self) -> usize {
        self.order_u
    }
    //}}}
    //{{{ fun: order_v
    pub fn order_v(&self) -> usize {
        self.order_v
    }
    //}}}
    //{{{ fun: knots_u
    pub fn knots_u(&self) -> &Vec<f64> {
        &self.knots_u
    }
    //}}}
    //{{{ fun: knots_v
    pub fn knots_v(&self) -> &Vec<f64> {
        &self.knots_v
    }
    //}}}
    //{{{ fun: num_points_u
    /// Number of control points in the `u`-direction of the control net.
    pub fn num_points_u(&self) -> usize {
        self.knots_u.len().saturating_sub(self.order_u + 1)
    }
    //}}}
    //{{{ fun: num_points_v
    /// Number of control points in the `v`-direction of the control net.
    pub fn num_points_v(&self) -> usize {
        self.knots_v.len().saturating_sub(self.order_v + 1)
    }
    //}}}
    //{{{ fun: map_u_curves
    /// Applies `op` to every iso-`v` row of the control net, viewed as a curve in `u`.
    ///
    /// Each `u`-direction operation on a tensor-product surface is the corresponding curve
    /// operation applied to each row of the control net over the shared `u`-knot vector: the
    /// recursion coefficients depend only on the knots, so every row yields the same new order
    /// and the same new knot vector. The `v`-direction data is untouched.
    fn map_u_curves<F>(
        &mut self,
        mut op: F,
    ) where
        F: FnMut(&mut BCurve),
    {
        let num_points_u = self.num_points_u();
        let num_points_v = self.num_points_v();

        let mut curves = Vec::with_capacity(num_points_v);
        for j in 0..num_points_v {
            let points_w = (0..num_points_u).map(|i| self.points_w[&[i, j]]).collect();
            let mut curve = BCurve::from_points_w(self.order_u, self.knots_u.clone(), points_w);
            op(&mut curve);
            curves.push(curve);
        }

        let num_points_u_after = curves[0].points_w().len();
        let mut points_w_after =
            DynArray::<Vec4, 2>::new(Vec4::zeros(), &[num_points_u_after, num_points_v]);
        for (j, curve) in curves.iter().enumerate() {
            for (i, point_w) in curve.points_w().iter().enumerate() {
                points_w_after[&[i, j]].copy_from(point_w);
            }
        }

        self.order_u = curves[0].order();
        self.knots_u = curves[0].knots().clone();
        self.points_w = points_w_after;
    }
    //}}}
    //{{{ fun: map_v_curves
    /// Applies `op` to every iso-`u` column of the control net, viewed as a curve in `v`.
    ///
    /// This is the `v`-direction counterpart of [`BSurface::map_u_curves`].
    fn map_v_curves<F>(
        &mut self,
        mut op: F,
    ) where
        F: FnMut(&mut BCurve),
    {
        let num_points_u = self.num_points_u();
        let num_points_v = self.num_points_v();

        let mut curves = Vec::with_capacity(num_points_u);
        for i in 0..num_points_u {
            let points_w = (0..num_points_v).map(|j| self.points_w[&[i, j]]).collect();
            let mut curve = BCurve::from_points_w(self.order_v, self.knots_v.clone(), points_w);
            op(&mut curve);
            curves.push(curve);
        }

        let num_points_v_after = curves[0].points_w().len();
        let mut points_w_after =
            DynArray::<Vec4, 2>::new(Vec4::zeros(), &[num_points_u, num_points_v_after]);
        for (i, curve) in curves.iter().enumerate() {
            for (j, point_w) in curve.points_w().iter().enumerate() {
                points_w_after[&[i, j]].copy_from(point_w);
            }
        }

        self.order_v = curves[0].order();
        self.knots_v = curves[0].knots().clone();
        self.points_w = points_w_after;
    }
    //}}}
    //{{{ fun: insert_knot_u
    /// Inserts the knot `u` into the `u`-knot vector `num_insertions` times.
    ///
    /// The surface is geometrically unchanged, it simply gains `num_insertions` control points in
    /// the `u`-direction.
    #[trace_fn]
    pub fn insert_knot_u(
        &mut self,
        u: f64,
        num_insertions: usize,
    ) {
        self.map_u_curves(|curve| curve.insert_knot(u, num_insertions));
    }
    //}}}
    //{{{ fun: with_inserted_knot_u
    pub fn with_inserted_knot_u(
        &self,
        u: f64,
        num_insertions: usize,
    ) -> Self {
        let mut new = self.clone();
        new.insert_knot_u(u, num_insertions);
        new
    }
    //}}}
    //{{{ fun: refine_knots_u
    /// Inserts all of the knots in `u`, which must be sorted, into the `u`-knot vector.
    #[trace_fn]
    pub fn refine_knots_u(
        &mut self,
        u: &[f64],
    ) {
        self.map_u_curves(|curve| curve.refine_knots(u));
    }
    //}}}
    //{{{ fun: with_refined_knots_u
    pub fn with_refined_knots_u(
        &self,
        u: &[f64],
    ) -> Self {
        let mut new = self.clone();
        new.refine_knots_u(u);
        new
    }
    //}}}
    //{{{ fun: elevate_degree_u
    /// Raises the order of the surface in the `u`-direction by `q`.
    #[trace_fn]
    pub fn elevate_degree_u(
        &mut self,
        q: usize,
    ) {
        self.map_u_curves(|curve| curve.elevate_degree(q));
    }
    //}}}
    //{{{ fun: with_elevated_degree_u
    pub fn with_elevated_degree_u(
        &self,
        q: usize,
    ) -> Self {
        let mut new = self.clone();
        new.elevate_degree_u(q);
        new
    }
    //}}}
    //{{{ fun: insert_knot_v
    /// Inserts the knot `v` into the `v`-knot vector `num_insertions` times.
    ///
    /// The surface is geometrically unchanged, it simply gains `num_insertions` control points in
    /// the `v`-direction.
    #[trace_fn]
    pub fn insert_knot_v(
        &mut self,
        v: f64,
        num_insertions: usize,
    ) {
        self.map_v_curves(|curve| curve.insert_knot(v, num_insertions));
    }
    //}}}
    //{{{ fun: with_inserted_knot_v
    pub fn with_inserted_knot_v(
        &self,
        v: f64,
        num_insertions: usize,
    ) -> Self {
        let mut new = self.clone();
        new.insert_knot_v(v, num_insertions);
        new
    }
    //}}}
    //{{{ fun: refine_knots_v
    /// Inserts all of the knots in `v`, which must be sorted, into the `v`-knot vector.
    #[trace_fn]
    pub fn refine_knots_v(
        &mut self,
        v: &[f64],
    ) {
        self.map_v_curves(|curve| curve.refine_knots(v));
    }
    //}}}
    //{{{ fun: with_refined_knots_v
    pub fn with_refined_knots_v(
        &self,
        v: &[f64],
    ) -> Self {
        let mut new = self.clone();
        new.refine_knots_v(v);
        new
    }
    //}}}
    //{{{ fun: elevate_degree_v
    /// Raises the order of the surface in the `v`-direction by `q`.
    #[trace_fn]
    pub fn elevate_degree_v(
        &mut self,
        q: usize,
    ) {
        self.map_v_curves(|curve| curve.elevate_degree(q));
    }
    //}}}
    //{{{ fun: with_elevated_degree_v
    pub fn with_elevated_degree_v(
        &self,
        q: usize,
    ) -> Self {
        let mut new = self.clone();
        new.elevate_degree_v(q);
        new
    }
    //}}}
    //{{{ fun: extract_isocurve_at_u
    pub fn extract_isocurve_at_u(
        &self,
        u: f64,
    ) -> BCurve {
        let (start_u, end_u, _) = non_zero_basis(&self.knots_u, u, self.order_u);
        let mut basis_funs_u = [0.0; PMAX_PLUS_ONE];
        splines::eval(&self.knots_u, u, self.order_u, &mut basis_funs_u);

        let p = self.order_v;
        let knots = self.knots_v.clone();
        let num_points_v = self.num_points_v();
        let mut points_w = Vec::<Vec4>::with_capacity(num_points_v);
        for i in 0..num_points_v {
            let mut new_point_w = Vec4::zeros();
            for j in start_u..end_u {
                new_point_w += basis_funs_u[j - start_u] * self.points_w[&[j, i]]
            }
            points_w.push(new_point_w);
        }
        BCurve::from_points_w(p, knots, points_w)
    }
    //}}}
    //{{{ fun: extract_isocurve_v
    pub fn extract_isocurve_at_v(
        &self,
        v: f64,
    ) -> BCurve {
        let (start_v, end_v, _) = non_zero_basis(&self.knots_v, v, self.order_v);
        let mut basis_funs_v = [0.0; PMAX_PLUS_ONE];
        splines::eval(&self.knots_v, v, self.order_v, &mut basis_funs_v);

        let p = self.order_u;
        let knots = self.knots_u.clone();
        let num_points_u = self.num_points_u();
        let mut points_w = Vec::<Vec4>::with_capacity(num_points_u);
        for i in 0..num_points_u {
            let mut new_point_w = Vec4::zeros();
            for j in start_v..end_v {
                new_point_w += basis_funs_v[j - start_v] * self.points_w[&[i, j]]
            }
            points_w.push(new_point_w);
        }
        BCurve::from_points_w(p, knots, points_w)
    }
    //}}}
}
//}}}
//{{{ impl: Surface for BSurface
impl Surface for BSurface {
    //{{{ fun: eval_point
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
    //}}}
    //{{{ fun: eval_derivative
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
    //}}}
    //{{{ fun: eval_all_derivatives
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
    //}}}
}
//}}}
//{{{ struct: VtkOptions
pub(crate) struct VtkOptions {
    u_range: (f64, f64),
    v_range: (f64, f64),
    num_points: (usize, usize),
    show_control_points: bool,
    isolines: Option<(Vec<f64>, Vec<f64>)>,
}
//}}}
//{{{ impl: ToVtk for BSurface
impl ToVtk for BSurface {
    type Options = VtkOptions;

    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        let VtkOptions {
            u_range,
            v_range,
            num_points,
            show_control_points,
            isolines,
        } = opts;

        add_surface_to_vtk(self, u_range, v_range, num_points, vtk_file);

        if show_control_points {
            let np_init = vtk_file.num_points() as u64;
            for (i, point_w) in self.points_w.iter().enumerate() {
                let point = inv_homog(point_w);
                vtk_file.add_point(point);
                vtk_file.add_cell(vec![np_init + i as u64], CellType::Vertex);
            }
            // do the u-directions
            for j in 0..self.num_points_v() {
                for i in 0..self.num_points_u() - 1 {
                    let idx1 = np_init + self.points_w.index_helper().lin_index(&[i, j]) as u64;
                    let idx2 = np_init + self.points_w.index_helper().lin_index(&[i + 1, j]) as u64;
                    vtk_file.add_cell(vec![idx1, idx2], CellType::Line);
                }
            }

            // do the v-directions
            for i in 0..self.num_points_u() {
                for j in 0..self.num_points_v() - 1 {
                    let idx1 = np_init + self.points_w.index_helper().lin_index(&[i, j]) as u64;
                    let idx2 = np_init + self.points_w.index_helper().lin_index(&[i, j + 1]) as u64;
                    vtk_file.add_cell(vec![idx1, idx2], CellType::Line);
                }
            }
        }

        if let Some((isolines_u, isolines_v)) = isolines {
            for u in isolines_u {
                let iso_curve_u = self.extract_isocurve_at_u(u);
                iso_curve_u.add_to_vtk(
                    BCurveVtkOptions {
                        ustart: u_range.0,
                        uend: u_range.1,
                        num_points: 10 * num_points.0,
                        show_control_points: false,
                    },
                    vtk_file,
                );
            }
            for v in isolines_v {
                let iso_curve_v = self.extract_isocurve_at_v(v);
                iso_curve_v.add_to_vtk(
                    BCurveVtkOptions {
                        ustart: v_range.0,
                        uend: v_range.1,
                        num_points: 10 * num_points.1,
                        show_control_points: false,
                    },
                    vtk_file,
                );
            }
        }
    }
}
//}}}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::*;

    fn reference_surface() -> BSurface {
        let order_u = 3;
        let order_v = 2;
        let knots_u = vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0];
        let knots_v = vec![0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0];

        let num_points_u = knots_u.len() - order_u - 1;
        let num_points_v = knots_v.len() - order_v - 1;
        let mut points = Vec::with_capacity(num_points_u * num_points_v);
        let mut weights = Vec::with_capacity(num_points_u * num_points_v);
        for i in 0..num_points_u {
            for j in 0..num_points_v {
                points.push(Vec3::from_col_slice(&[
                    i as f64,
                    1.5 * j as f64,
                    ((i * j) % 5) as f64 - 2.0,
                ]));
                weights.push(1.0 + 0.25 * ((i + 2 * j) % 4) as f64);
            }
        }

        BSurface::new(order_u, order_v, knots_u, knots_v, points, weights)
    }

    #[test]
    fn test_vtk_output() {
        let surface = reference_surface();
        let mut vtk_file = VtkUgridFile::default();
        surface.add_to_vtk(
            VtkOptions {
                u_range: (0.0, 1.0),
                v_range: (0.0, 1.0),
                num_points: (100, 100),
                show_control_points: true,
                isolines: Some((vec![0.2, 0.2, 0.6, 0.8], vec![0.2, 0.2, 0.6, 0.8])),
            },
            &mut vtk_file,
        );
        vtk_file
            .save("bsurface.vtu")
            .expect("B-surface VTK output should be written");
    }
}
//}}}
