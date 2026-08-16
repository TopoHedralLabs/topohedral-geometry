//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{common::Vec3, vtk::VtkUgridFile};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
use topohedral_tracing::*;
use vtkio::model::CellType;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ trait: Curve
pub trait Curve {
    fn param_range(&self) -> (f64, f64);

    fn in_range(
        &self,
        u: f64,
    ) -> bool {
        let range = self.param_range();
        u >= range.0 && u <= range.1
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
        let mut curve_dirivatives: [Vec3; 3] = std::array::from_fn(|_| Vec3::zeros());
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
//}}}
//{{{ struct: CurveVtkOptions
pub(crate) struct CurveVtkOptions {
    pub(crate) ustart: f64,
    pub(crate) uend: f64,
    pub(crate) num_points: usize,
}
//}}}

#[trace_fn]
pub(crate) fn add_curve_to_vtk<C: Curve + ?Sized>(
    curve: &C,
    ustart: f64,
    uend: f64,
    num_points: usize,
    vtk_file: &mut VtkUgridFile,
) {
    assert!(num_points >= 2, "at least two points are required");

    let np_init = vtk_file.num_points() as u64;
    let num_segments = (num_points - 1) as u64;
    let du = (uend - ustart) / num_segments as f64;
    let mut ui = ustart;
    for _ in 0..num_points {
        let xi = curve.eval_point(ui);
        vtk_file.add_point(xi);
        ui += du;
    }

    let curve_cell: Vec<u64> = (np_init..(np_init + num_points as u64)).collect();
    vtk_file.add_cell(curve_cell, CellType::PolyLine);
}
