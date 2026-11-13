//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Sphere, Surface, Vec3};
//}}}
//{{{ std imports
use std::f64::consts::{FRAC_PI_2, PI};
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

fn assert_vectors_eq(
    actual: &Vec3,
    expected: &Vec3,
) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected.iter()) {
        assert_abs_diff_eq!(actual, expected, epsilon = TOLERANCE);
    }
}

fn reference_sphere() -> Sphere {
    Sphere::new(vector([1.0, -2.0, 3.0]), 2.0, vector([0.0, 0.0, 1.0]))
}

#[test]
fn dimension_and_cardinal_points() {
    let sphere = reference_sphere();

    let c = vector([1.0, -2.0, 3.0]);
    let r = 2.0;
    let ex = vector([0.0, 1.0, 0.0]);
    let ey = vector([1.0, 0.0, 0.0]);
    let ez = vector([0.0, 0.0, 1.0]);

    let p1: Vec3 = c + r * ex;
    assert_vector_eq(&sphere.eval_point(0.0, 0.0), p1.as_ref());
    let p2: Vec3 = c + r * ey;
    assert_vector_eq(&sphere.eval_point(FRAC_PI_2, 0.0), p2.as_ref());
    let p3: Vec3 = c + r * ez;
    assert_vector_eq(&sphere.eval_point(0.0, FRAC_PI_2), p3.as_ref());
    let p4: Vec3 = c - r * ez;
    assert_vector_eq(&sphere.eval_point(0.0, -FRAC_PI_2), p4.as_ref());
    let p5: Vec3 = c - r * ex;
    assert_vector_eq(&sphere.eval_point(PI, 0.0), p5.as_ref());
}

#[test]
fn eval_derivative_at_zero_order_matches_eval_point() {
    let sphere = reference_sphere();

    assert_vectors_eq(
        &sphere.eval_derivative(0.7, -0.4, 0, 0),
        &sphere.eval_point(0.7, -0.4),
    );
}

#[test]
fn first_derivatives_are_correct() {
    let sphere = reference_sphere();

    assert_vector_eq(&sphere.eval_derivative(0.0, 0.0, 1, 0), &[2.0, 0.0, 0.0]);
    assert_vector_eq(&sphere.eval_derivative(0.0, 0.0, 0, 1), &[0.0, 0.0, 2.0]);
    assert_vector_eq(
        &sphere.eval_derivative(FRAC_PI_2, 0.0, 1, 0),
        &[0.0, -2.0, 0.0],
    );
}

#[test]
fn second_derivatives_are_correct() {
    let sphere = reference_sphere();

    assert_vector_eq(&sphere.eval_derivative(0.0, 0.0, 2, 0), &[0.0, -2.0, 0.0]);
    assert_vector_eq(&sphere.eval_derivative(0.0, 0.0, 0, 2), &[0.0, -2.0, 0.0]);
}

#[test]
fn mixed_derivative_is_correct() {
    let sphere = reference_sphere();

    assert_vector_eq(
        &sphere.eval_derivative(0.0, FRAC_PI_2, 1, 1),
        &[-2.0, 0.0, 0.0],
    );
}

#[test]
fn all_derivatives_use_column_major_parameter_order() {
    let sphere = reference_sphere();
    let mut derivatives = vec![vector([0.0, 0.0, 0.0]); 4];

    sphere.eval_all_derivatives(0.0, 0.0, 1, 1, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[1.0, 0.0, 3.0]);
    assert_vector_eq(&derivatives[1], &[2.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[0.0, 0.0, 2.0]);
    assert_vector_eq(&derivatives[3], &[0.0, 0.0, 2.0]);
}

#[test]
fn normal_is_correct() {
    let sphere = reference_sphere();

    assert_vector_eq(&sphere.eval_normal(0.0, 0.0, false), &[0.0, -4.0, 0.0]);
    assert_vector_eq(&sphere.eval_normal(0.0, 0.0, true), &[0.0, -1.0, 0.0]);
}

#[test]
#[should_panic]
fn rejects_non_unit_polar_axis() {
    Sphere::new(vector([0.0, 0.0, 0.0]), 1.0, vector([0.0, 0.0, 2.0]));
}
