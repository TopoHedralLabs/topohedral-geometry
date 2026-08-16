//! A hyperbolic conic curve embedded in three-dimensional space.

use crate::{
    common::Vec3,
    curves::common::add_curve_to_vtk,
    utils::{is_unit_vector, orthogonal_vector},
    vtk::{ToVtk, VtkUgridFile},
    Curve,
};
use topohedral_linalg::VectorOps;

/// A hyperbola parameterized as
/// `center + x_radius * cosh(u) * x_axis + y_radius * sinh(u) * y_axis`.
pub struct Hyperbola {
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    x_radius: f64,
    y_radius: f64,
}

impl Hyperbola {
    /// Creates a hyperbola in the plane defined by `center` and `normal`.
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

impl Curve for Hyperbola {
    fn param_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        let Self {
            center,
            x_axis,
            y_axis,
            x_radius,
            y_radius,
        } = self;

        (center + x_radius * u.cosh() * x_axis + y_radius * u.sinh() * y_axis).into()
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
            x_axis,
            y_axis,
            x_radius,
            y_radius,
            ..
        } = self;

        if m.is_multiple_of(2) {
            (x_radius * u.cosh() * x_axis + y_radius * u.sinh() * y_axis).into()
        } else {
            (x_radius * u.sinh() * x_axis + y_radius * u.cosh() * y_axis).into()
        }
    }
}

pub(crate) struct VtkOptions {
    pub(crate) ustart: f64,
    pub(crate) uend: f64,
    pub(crate) num_points: usize,
}

impl ToVtk for Hyperbola {
    type Options = VtkOptions;

    fn add_to_vtk(
        &self,
        VtkOptions {
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
        let curve = Hyperbola::new(
            Vec3::from_col_slice(&[0.0, 0.0, 0.0]),
            Vec3::from_col_slice(&[0.0, 0.0, 1.0]),
            1.0,
            0.5,
        );
        let mut vtk_file = VtkUgridFile::default();
        curve.add_to_vtk(
            VtkOptions {
                ustart: -1.0,
                uend: 1.0,
                num_points: 100,
            },
            &mut vtk_file,
        );
        vtk_file
            .save("hyperbola.vtu")
            .expect("hyperbola VTK output should be written");
    }
}
