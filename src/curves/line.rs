//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use crate::{
    common::Vec3,
    curves::common::{add_curve_to_vtk, CurveVtkOptions},
    utils::is_unit_vector,
    vtk::{ToVtk, VtkUgridFile},
    Curve,
};
//}}}
//{{{ std imports
//}}}
//--------------------------------------------------------------------------------------------------

pub struct Line {
    origin: Vec3,
    direction: Vec3,
}

impl Line {
    pub fn new(
        origin: Vec3,
        direction: Vec3,
    ) -> Self {
        assert!(is_unit_vector(&direction));
        Self { origin, direction }
    }
}

impl Curve for Line {
    fn eval_point(
        &self,
        u: f64,
    ) -> Vec3 {
        self.origin + u * self.direction
    }

    fn eval_derivative(
        &self,
        u: f64,
        m: usize,
    ) -> Vec3 {
        if m == 0 {
            self.eval_point(u)
        } else if m == 1 {
            self.direction
        } else {
            Vec3::zeros()
        }
    }

    fn param_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

impl ToVtk for Line {
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
        let curve = Line::new(
            Vec3::from_col_slice(&[0.0, 0.0, 0.0]),
            Vec3::from_col_slice(&[1.0, 0.0, 0.0]),
        );
        let mut vtk_file = VtkUgridFile::default();
        curve.add_to_vtk(
            CurveVtkOptions {
                ustart: -1.0,
                uend: 1.0,
                num_points: 10,
            },
            &mut vtk_file,
        );
        vtk_file
            .save("line.vtu")
            .expect("line VTK output should be written");
    }
}
