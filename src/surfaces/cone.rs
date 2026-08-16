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

//{{{ struct: Cone
pub struct Cone {
    center: Vec3,
    axis: Vec3,
    radius_bottom: f64,
    tan_alpha: f64,
    x_axis: Vec3,
    y_axis: Vec3,
}
//}}}
//{{{ impl: Cone
impl Cone {
    pub fn new(
        center: Vec3,
        axis: Vec3,
        height: f64,
        radius_bottom: f64,
        radius_top: f64,
    ) -> Self {
        assert!(is_unit_vector(&axis));
        assert!(radius_top > 0.0);
        assert!(height > 0.0);
        assert!(radius_bottom > radius_top);
        let tan_alpha = (radius_bottom - radius_top) / height;
        let x_axis = orthogonal_vector(&axis);
        let y_axis = x_axis.cross(&axis);
        Self {
            center,
            axis,
            radius_bottom,
            tan_alpha,
            x_axis,
            y_axis,
        }
    }

    fn radius(
        &self,
        v: f64,
    ) -> f64 {
        self.radius_bottom - v * self.tan_alpha
    }

    fn eval_latitude_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        self.x_axis * cos_derivative(u, m) + self.y_axis * sin_derivative(u, m)
    }
}
//}}}
//{{{ impl: Surface
impl Surface for Cone {
    fn eval_point(
        &self,
        u: f64,
        v: f64,
    ) -> Vec3 {
        let Self {
            center: c,
            axis: ez,
            radius_bottom: _,
            tan_alpha: _,
            x_axis: _,
            y_axis: _,
        } = self;
        let r = self.radius(v);
        (c + r * self.eval_latitude_derivative(u, 0) + v * ez).into()
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

        if n == 0 {
            return self.radius(v) * self.eval_latitude_derivative(u, m);
        }

        if n == 1 {
            let axial_derivative = if m == 0 { self.axis } else { Vec3::zeros() };
            return -self.tan_alpha * self.eval_latitude_derivative(u, m) + axial_derivative;
        }

        Vec3::zeros()
    }
}
//}}}

impl ToVtk for Cone {
    type Options = SurfaceVtkOptions;

    fn add_to_vtk(
        &self,
        opts: Self::Options,
        vtk_file: &mut VtkUgridFile,
    ) {
        add_surface_to_vtk(self, opts.u_range, opts.v_range, opts.num_points, vtk_file);
    }
}
