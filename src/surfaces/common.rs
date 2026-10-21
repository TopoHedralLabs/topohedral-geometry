//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{common::IndexHelper, vtk::VtkUgridFile, Vec3};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
use topohedral_tracing::*;
use vtkio::model::CellType;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ trait: Surface
pub trait Surface {
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
//}}}
//{{{ struct: SurfaceVtkOptions
pub(crate) struct SurfaceVtkOptions {
    pub(crate) u_range: (f64, f64),
    pub(crate) v_range: (f64, f64),
    pub(crate) num_points: (usize, usize),
}
//}}}
//{{{ fun: add_surface_to_vtk
#[trace_fn]
pub(crate) fn add_surface_to_vtk<S: Surface + ?Sized>(
    surface: &S,
    u_range: (f64, f64),
    v_range: (f64, f64),
    num_points: (usize, usize),
    vtk_file: &mut VtkUgridFile,
) {
    let np_init = vtk_file.num_points();
    let np_u = num_points.0;
    let np_v = num_points.1;
    let num_seg_u = np_u - 1;
    let num_seg_v = np_v - 1;
    let du = (u_range.1 - u_range.0) / num_seg_u as f64;
    let dv = (v_range.1 - v_range.0) / num_seg_v as f64;

    for j in 0..num_points.1 {
        let v = v_range.0 + j as f64 * dv;
        for i in 0..num_points.0 {
            let u = u_range.0 + i as f64 * du;
            let point = surface.eval_point(u, v);
            vtk_file.add_point(point);
        }
    }

    for j in 0..num_seg_v {
        for i in 0..num_seg_u {
            let v0 = (np_init + np_u * j + i) as u64;
            let v1 = (np_init + np_u * j + i + 1) as u64;
            let v2 = (np_init + np_u * (j + 1) + i + 1) as u64;
            let v3 = (np_init + np_u * (j + 1) + i) as u64;
            let tri0 = vec![v0, v1, v2];
            let tri1 = vec![v0, v2, v3];
            vtk_file.add_cell(tri0, CellType::Triangle);
            vtk_file.add_cell(tri1, CellType::Triangle);
        }
    }
}
//}}}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{vtk::ToVtk, Cone, Cylinder, Plane, Sphere, Torus};

    fn assert_vtk_export<S>(surface: S)
    where
        S: Surface + ToVtk<Options = SurfaceVtkOptions>,
    {
        let mut vtk_file = VtkUgridFile::default();
        surface.add_to_vtk(
            SurfaceVtkOptions {
                u_range: (0.0, 1.0),
                v_range: (0.0, 1.0),
                num_points: (3, 2),
            },
            &mut vtk_file,
        );

        assert_eq!(vtk_file.num_points(), 6);
        assert_eq!(vtk_file.num_cells(), 4);
    }

    #[test]
    fn parametric_surfaces_export_to_vtk() {
        let origin = Vec3::zeros();
        let z_axis = Vec3::from_col_slice(&[0.0, 0.0, 1.0]);

        assert_vtk_export(Sphere::default());
        assert_vtk_export(Torus::default());
        assert_vtk_export(Cone::new(origin, z_axis, 2.0, 2.0, 1.0));
        assert_vtk_export(Cylinder::new(origin, z_axis, 1.0));
        assert_vtk_export(Plane::new(origin, z_axis));
    }
}
//}}}
