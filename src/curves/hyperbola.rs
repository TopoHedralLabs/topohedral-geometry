//! A hyperbolic conic curve embedded in three-dimensional space.

use crate::common::Vec3;
use crate::utils::{is_unit_vector, orthogonal_vector};
use crate::Curve;
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
