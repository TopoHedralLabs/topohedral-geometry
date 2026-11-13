//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::Vec3,
    curves::common::{add_curve_to_vtk, CurveVtkOptions},
    utils::{is_unit_vector, orthogonal_vector},
    vtk::{ToVtk, VtkUgridFile},
    Curve,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Ellipse {
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    x_radius: f64,
    y_radius: f64,
}

impl Ellipse {
    pub fn new(
        center: Vec3,
        normal: Vec3,
        x_radius: f64,
        y_radius: f64,
    ) -> Self {
        assert!(is_unit_vector(&normal));
        assert!(x_radius > 0.0);
        assert!(y_radius > 0.0);

        let x_axis = orthogonal_vector(&normal);
        let y_axis = x_axis.cross(&normal);
        Self {
            center,
            x_axis,
            y_axis,
            x_radius,
            y_radius,
        }
    }
}

impl Curve for Ellipse {
    fn param_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            x_axis,
            y_axis,
            x_radius: rx,
            y_radius: ry,
        } = self;

        (c + rx * u.cos() * x_axis + ry * u.sin() * y_axis).into()
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        let Self {
            center: _,
            x_axis,
            y_axis,
            x_radius: rx,
            y_radius: ry,
        } = self;

        if m == 0 {
            return self.eval_point(u);
        }

        if m.is_multiple_of(2) {
            let is_odd = ((m / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (rx * u.cos() * x_axis + ry * u.sin() * y_axis)).into()
        } else {
            let is_odd = (((m - 1) / 2) % 2) == 1;
            let sign = if is_odd { -1.0 } else { 1.0 };

            (sign * (-rx * u.sin() * x_axis + ry * u.cos() * y_axis)).into()
        }
    }
}

impl ToVtk for Ellipse {
    type Options = CurveVtkOptions;

    fn add_to_vtk(
        &self,
        CurveVtkOptions {
            ustart,
            uend,
            num_points,
        }: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        add_curve_to_vtk(self, ustart, uend, num_points, vtk_file);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vtk_output() {
        let curve = Ellipse::new(
            Vec3::from_col_slice(&[0.0, 0.0, 0.0]),
            Vec3::from_col_slice(&[0.0, 0.0, 1.0]),
            2.0,
            1.0,
        );
        let mut vtk_file = VtkUgridFile::default();
        curve.add_to_vtk(
            CurveVtkOptions {
                ustart: 0.0,
                uend: std::f64::consts::TAU,
                num_points: 100,
            },
            &mut vtk_file,
        );
        vtk_file
            .save("ellipse.vtu")
            .expect("ellipse VTK output should be written");
    }
}
