//! Integration tests for parabolic curves.

use approx::assert_abs_diff_eq;
use topohedral_geometry::{Curve, Parabola, Vec3};

const TOLERANCE: f64 = 1.0e-12;

fn point(values: [f64; 3]) -> Vec3 {
    Vec3::from_col_slice(&values)
}

fn assert_vector_eq(
    actual: &Vec3,
    expected: [f64; 3],
) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(*actual, expected, epsilon = TOLERANCE);
    }
}

fn reference_parabola() -> Parabola {
    Parabola::new(
        point([1.0, -2.0, 3.0]),
        point([0.0, 0.0, 1.0]),
        point([1.0, 0.0, 0.0]),
        0.5,
    )
}

#[test]
fn evaluation_follows_the_parabola_parameterization() {
    let parabola = reference_parabola();

    assert_eq!(parabola.dimension(), 3);
    assert_vector_eq(&parabola.eval_point(0.0), [1.0, -2.0, 3.0]);
    assert_vector_eq(&parabola.eval_point(1.0), [3.0, -3.0, 3.0]);
    assert_vector_eq(&parabola.eval_point(2.0), [5.0, -6.0, 3.0]);
}

#[test]
fn derivatives_match_the_quadratic_curve() {
    let parabola = reference_parabola();

    assert_vector_eq(&parabola.eval_derivative(2.0, 0), [5.0, -6.0, 3.0]);
    assert_vector_eq(&parabola.eval_derivative(2.0, 1), [2.0, -4.0, 0.0]);
    assert_vector_eq(&parabola.eval_derivative(2.0, 2), [0.0, -2.0, 0.0]);
    assert_vector_eq(&parabola.eval_derivative(2.0, 3), [0.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_and_frenet_frame_are_consistent() {
    let parabola = reference_parabola();
    let mut derivatives = vec![Vec3::zeros(); 4];

    parabola.eval_all_derivatives(0.0, 3, &mut derivatives);

    assert_vector_eq(&derivatives[0], [1.0, -2.0, 3.0]);
    assert_vector_eq(&derivatives[1], [2.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], [0.0, -2.0, 0.0]);
    assert_vector_eq(&derivatives[3], [0.0, 0.0, 0.0]);
    assert_vector_eq(&parabola.eval_tangent(0.0, true), [1.0, 0.0, 0.0]);
    assert_vector_eq(&parabola.eval_normal(0.0, true), [0.0, -1.0, 0.0]);
    assert_vector_eq(&parabola.eval_binormal(0.0, true), [0.0, 0.0, -1.0]);
}
