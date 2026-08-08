//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Plane, Surface, Vec3};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::assert_abs_diff_eq;
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

const TOLERANCE: f64 = 1.0e-12;

fn vector(values: [f64; 3]) -> Vec3 {
    Vec3::from_col_slice(&values)
}

fn assert_vector_eq(
    actual: &Vec3,
    expected: &[f64],
) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = TOLERANCE);
    }
}

fn reference_plane() -> Plane {
    Plane::new(vector([1.0, -2.0, 3.0]), vector([0.0, 0.0, 1.0]))
}

#[test]
fn dimension_and_point_evaluation() {
    let plane = reference_plane();

    assert_eq!(plane.dimension(), 3);
    assert_vector_eq(&plane.eval_point(0.0, 0.0), &[1.0, -2.0, 3.0]);
    assert_vector_eq(&plane.eval_point(2.0, -3.0), &[-2.0, 0.0, 3.0]);
}

#[test]
fn first_derivatives_are_parameter_axes() {
    let plane = reference_plane();

    assert_vector_eq(&plane.eval_derivative(4.0, -1.5, 1, 0), &[0.0, 1.0, 0.0]);
    assert_vector_eq(&plane.eval_derivative(4.0, -1.5, 0, 1), &[1.0, 0.0, 0.0]);
}

#[test]
fn higher_and_mixed_derivatives_are_zero() {
    let plane = reference_plane();

    assert_vector_eq(&plane.eval_derivative(0.0, 0.0, 2, 0), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&plane.eval_derivative(0.0, 0.0, 0, 2), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&plane.eval_derivative(0.0, 0.0, 1, 1), &[0.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_use_column_major_parameter_order() {
    let plane = reference_plane();
    let mut derivatives = vec![vector([99.0, 99.0, 99.0]); 4];

    plane.eval_all_derivatives(2.0, -3.0, 1, 1, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[-2.0, 0.0, 3.0]);
    assert_vector_eq(&derivatives[1], &[0.0, 1.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[1.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[0.0, 0.0, 0.0]);
}

#[test]
fn normal_is_constant_and_has_plane_orientation() {
    let plane = reference_plane();

    assert_vector_eq(&plane.eval_normal(0.0, 0.0, false), &[0.0, 0.0, -1.0]);
    assert_vector_eq(&plane.eval_normal(2.0, -3.0, true), &[0.0, 0.0, -1.0]);
}

#[test]
#[should_panic]
fn rejects_non_unit_normal() {
    Plane::new(vector([0.0, 0.0, 0.0]), vector([0.0, 0.0, 2.0]));
}
