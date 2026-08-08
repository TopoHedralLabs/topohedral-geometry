//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Cylinder, Surface, Vec3};
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

fn reference_cylinder() -> Cylinder {
    Cylinder::new(vector([1.0, -2.0, 3.0]), vector([0.0, 0.0, 1.0]), 2.0)
}

#[test]
fn dimension_and_cardinal_points() {
    let cylinder = reference_cylinder();

    assert_eq!(cylinder.dimension(), 3);
    assert_vector_eq(&cylinder.eval_point(0.0, 0.0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(
        &cylinder.eval_point(std::f64::consts::FRAC_PI_2, 4.0),
        &[3.0, -2.0, 7.0],
    );
    assert_vector_eq(
        &cylinder.eval_point(std::f64::consts::PI, -1.0),
        &[1.0, -4.0, 2.0],
    );
}

#[test]
fn angular_derivatives_repeat_every_four_orders() {
    let cylinder = reference_cylinder();

    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 0, 0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 1, 0), &[2.0, 0.0, 0.0]);
    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 2, 0), &[0.0, -2.0, 0.0]);
    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 3, 0), &[-2.0, 0.0, 0.0]);
    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 4, 0), &[0.0, 2.0, 0.0]);
    assert_vector_eq(&cylinder.eval_derivative(0.0, 0.0, 5, 0), &[2.0, 0.0, 0.0]);
}

#[test]
fn axial_and_mixed_derivatives_are_correct() {
    let cylinder = reference_cylinder();

    assert_vector_eq(&cylinder.eval_derivative(1.2, -3.0, 0, 1), &[0.0, 0.0, 1.0]);
    assert_vector_eq(&cylinder.eval_derivative(1.2, -3.0, 0, 2), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&cylinder.eval_derivative(1.2, -3.0, 1, 1), &[0.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_use_column_major_parameter_order() {
    let cylinder = reference_cylinder();
    let mut derivatives = vec![vector([99.0, 99.0, 99.0]); 6];

    cylinder.eval_all_derivatives(0.0, 4.0, 2, 1, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[1.0, 0.0, 7.0]);
    assert_vector_eq(&derivatives[1], &[2.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[0.0, -2.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[0.0, 0.0, 1.0]);
    assert_vector_eq(&derivatives[4], &[0.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[5], &[0.0, 0.0, 0.0]);
}

#[test]
fn normal_is_radial_and_constant_along_the_axis() {
    let cylinder = reference_cylinder();

    assert_vector_eq(&cylinder.eval_normal(0.0, 0.0, false), &[0.0, -2.0, 0.0]);
    assert_vector_eq(&cylinder.eval_normal(0.0, 4.0, true), &[0.0, -1.0, 0.0]);
}

#[test]
#[should_panic]
fn rejects_zero_radius() {
    Cylinder::new(vector([0.0, 0.0, 0.0]), vector([0.0, 0.0, 1.0]), 0.0);
}

#[test]
#[should_panic]
fn rejects_negative_radius() {
    Cylinder::new(vector([0.0, 0.0, 0.0]), vector([0.0, 0.0, 1.0]), -1.0);
}
