//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::common::Vec3;
use crate::curves::common::{add_curve_to_vtk, CurveVtkOptions};
use crate::utils::{is_unit_vector, orthogonal_vector};
use crate::vtk::{ToVtk, VtkUgridFile};
use crate::Curve;
//}}}
//{{{ std imports
use std::f64;
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
use vtkio::model::CellType;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Circle {
    center: Vec3,
    radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Circle {
    pub fn new(
        center: Vec3,
        normal: Vec3,
        radius: f64,
    ) -> Self {
        assert!(is_unit_vector(&normal));
        assert!(radius > 0.0);
        let x_axis = orthogonal_vector(&normal);
        let y_axis = x_axis.cross(&normal);
        Self {
            center,
            radius,
            x_axis,
            y_axis,
        }
    }
}

impl Curve for Circle {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            radius: r,
            x_axis: ex,
            y_axis: ey,
        } = self;

        (c + r * u.cos() * ex + r * u.sin() * ey).into()
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            return self.eval_point(u);
        }

        let Self {
            center: _,
            radius: r,
            x_axis: ex,
            y_axis: ey,
        } = self;
        if m.is_multiple_of(2) {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            (sign * r * (u.cos() * ex + u.sin() * ey)).into()
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };
            (sign * r * (-u.sin() * ex + u.cos() * ey)).into()
        }
    }

    fn param_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

pub(crate) struct VtkOptions {
    pub(crate) curve: CurveVtkOptions,
    pub(crate) show_center: bool,
}

impl ToVtk for Circle {
    type Options = VtkOptions;

    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        let VtkOptions {
            curve:
                CurveVtkOptions {
                    ustart,
                    uend,
                    num_points,
                },
            show_center,
        } = opts;

        add_curve_to_vtk(self, ustart, uend, num_points, vtk_file);
        if show_center {
            let center_id = vtk_file.num_points() as u64;
            vtk_file.add_point(self.center);
            vtk_file.add_cell(vec![center_id], CellType::Vertex);
        }
    }
}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vtk_output() {
        let curve = Circle::new(
            Vec3::from_col_slice(&[0.0, 0.0, 0.0]),
            Vec3::from_col_slice(&[0.0, 0.0, 1.0]),
            1.0,
        );
        let mut vtk_file = VtkUgridFile::default();
        curve.add_to_vtk(
            VtkOptions {
                curve: CurveVtkOptions {
                    ustart: 0.0,
                    uend: std::f64::consts::TAU,
                    num_points: 100,
                },
                show_center: false,
            },
            &mut vtk_file,
        );
        vtk_file
            .save("circle.vtu")
            .expect("circle VTK output should be written");
    }
}
//}}}
