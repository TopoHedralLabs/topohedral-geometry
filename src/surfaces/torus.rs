//! A ring torus with a configurable center, axis, and major and minor radii.
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    surfaces::common::{add_surface_to_vtk, SurfaceVtkOptions},
    utils::{cos_derivative, is_unit_vector, orthogonal_vector, sin_derivative},
    vtk::{ToVtk, VtkUgridFile},
    Surface, Vec3,
};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct Torus {
    center: Vec3,
    major_radius: f64,
    minor_radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
    z_axis: Vec3,
}

impl Default for Torus {
    fn default() -> Self {
        let ex = Vec3::from_col_slice(&[1.0, 0.0, 0.0]);
        let ey = Vec3::from_col_slice(&[0.0, 1.0, 0.0]);
        let ez = Vec3::from_col_slice(&[0.0, 0.0, 1.0]);
        Self {
            center: Vec3::zeros(),
            major_radius: 1.0,
            minor_radius: 0.1,
            x_axis: ex,
            y_axis: ey,
            z_axis: ez,
        }
    }
}

impl Torus {
    /// Creates a ring torus with `major_radius > minor_radius > 0` and a unit axis.
    ///
    /// Both surface parameters are angles in radians with period `2 * PI`:
    /// `u` runs around the axis and `v` runs around the tube.
    pub fn new(
        center: Vec3,
        axis: Vec3,
        major_radius: f64,
        minor_radius: f64,
    ) -> Self {
        assert!(major_radius > 0.0 && minor_radius > 0.0);
        assert!(major_radius > minor_radius);
        assert!(is_unit_vector(&axis));

        let xaxis = orthogonal_vector(&axis);
        let yaxis = xaxis.cross(&axis);
        Self {
            center,
            major_radius,
            minor_radius,
            x_axis: xaxis,
            y_axis: yaxis,
            z_axis: axis,
        }
    }

    fn eval_latitude_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        self.x_axis * cos_derivative(u, m) + self.y_axis * sin_derivative(u, m)
    }
}

impl Surface for Torus {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        self.center
            + (self.major_radius + self.minor_radius * v.cos())
                * self.eval_latitude_derivative(u, 0)
            + self.minor_radius * v.sin() * self.z_axis
    }

    fn eval_derivative(
        &self,
        u: f64,
        v: f64,
        m: usize,
        n: usize,
    ) -> Vec3 {
        if m == 0 && n == 0 {
            return self.eval_point(u, v);
        }

        let radial_factor = if n == 0 {
            self.major_radius + self.minor_radius * v.cos()
        } else {
            self.minor_radius * cos_derivative(v, n)
        };
        let radial = radial_factor * self.eval_latitude_derivative(u, m);
        if m == 0 {
            radial + self.minor_radius * sin_derivative(v, n) * self.z_axis
        } else {
            radial
        }
    }
}

impl ToVtk for Torus {
    type Options = SurfaceVtkOptions;

    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        add_surface_to_vtk(self, opts.u_range, opts.v_range, opts.num_points, vtk_file);
    }
}
