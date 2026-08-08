//! Integration tests for ellipse curves.

use approx::assert_abs_diff_eq;
use topohedral_geometry::{Curve, Ellipse, Vec3};
use topohedral_linalg::VectorOps;

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

fn reference_ellipse() -> Ellipse {
    Ellipse::new(vector([1.0, -2.0, 3.0]), vector([0.0, 0.0, 1.0]), 2.0, 3.0)
}

#[test]
fn dimension_and_cardinal_points() {
    let ellipse = reference_ellipse();

    assert_eq!(ellipse.dimension(), 3);
    assert_vector_eq(&ellipse.eval_point(0.0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(
        &ellipse.eval_point(std::f64::consts::FRAC_PI_2),
        &[4.0, -2.0, 3.0],
    );
    assert_vector_eq(&ellipse.eval_point(std::f64::consts::PI), &[1.0, -4.0, 3.0]);
    assert_vector_eq(
        &ellipse.eval_point(3.0 * std::f64::consts::FRAC_PI_2),
        &[-2.0, -2.0, 3.0],
    );
}

#[test]
fn derivatives_follow_the_ellipse_periodicity() {
    let ellipse = reference_ellipse();

    assert_vector_eq(&ellipse.eval_derivative(0.0, 0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(&ellipse.eval_derivative(0.0, 1), &[3.0, 0.0, 0.0]);
    assert_vector_eq(&ellipse.eval_derivative(0.0, 2), &[0.0, -2.0, 0.0]);
    assert_vector_eq(&ellipse.eval_derivative(0.0, 3), &[-3.0, 0.0, 0.0]);
    assert_vector_eq(&ellipse.eval_derivative(0.0, 4), &[0.0, 2.0, 0.0]);
}

#[test]
fn all_derivatives_are_written_to_the_output_slice() {
    let ellipse = reference_ellipse();
    let mut derivatives = vec![vector([99.0, 99.0, 99.0]); 5];

    ellipse.eval_all_derivatives(0.0, 4, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[1.0, 0.0, 3.0]);
    assert_vector_eq(&derivatives[1], &[3.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[0.0, -2.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[-3.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[4], &[0.0, 2.0, 0.0]);
}

#[test]
fn tangent_and_frame_follow_the_ellipse_orientation() {
    let ellipse = reference_ellipse();

    assert_vector_eq(&ellipse.eval_tangent(0.0, true), &[1.0, 0.0, 0.0]);
    assert_vector_eq(&ellipse.eval_normal(0.0, true), &[0.0, -1.0, 0.0]);
    assert_vector_eq(&ellipse.eval_binormal(0.0, true), &[0.0, 0.0, -1.0]);
}

#[test]
#[should_panic]
fn rejects_non_positive_radii() {
    Ellipse::new(vector([0.0, 0.0, 0.0]), vector([0.0, 0.0, 1.0]), 0.0, 1.0);
}
