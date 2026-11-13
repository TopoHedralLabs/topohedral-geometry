//! Short Description of module
//!
//! Longer description of module
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

pub struct Cylinder {
    center: Vec3,
    axis: Vec3,
    radius: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}

impl Cylinder {
    pub fn new(
        center: Vec3,
        axis: Vec3,
        radius: f64,
    ) -> Self {
        assert!(radius > 0.0);
        assert!(is_unit_vector(&axis));

        let x_axis = orthogonal_vector(&axis);
        let y_axis = x_axis.cross(&axis);

        Self {
            center,
            axis,
            radius,
            x_axis,
            y_axis,
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

impl Surface for Cylinder {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            axis: ez,
            radius: r,
            x_axis: _,
            y_axis: _,
        } = self;
        (c + *r * self.eval_latitude_derivative(u, 0) + v * ez).into()
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

        if m == 0 && n == 1 {
            return self.axis;
        }

        if m > 0 && n == 0 {
            return self.radius * self.eval_latitude_derivative(u, m);
        }
        // all other configurations, mixed derivatives, v derivative > 1 all zero
        Vec3::zeros()
    }
}

impl ToVtk for Cylinder {
    type Options = SurfaceVtkOptions;

    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        add_surface_to_vtk(self, opts.u_range, opts.v_range, opts.num_points, vtk_file);
    }
}
