//! Integration tests for hyperbolic curves.

use approx::assert_abs_diff_eq;
use topohedral_geometry::{Curve, Hyperbola, Vec3};

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

fn reference_hyperbola() -> Hyperbola {
    Hyperbola::new(point([1.0, -2.0, 3.0]), point([0.0, 0.0, 1.0]), 2.0, 3.0)
}

#[test]
fn evaluation_follows_the_hyperbolic_parameterization() {
    let hyperbola = reference_hyperbola();

    assert_vector_eq(&hyperbola.eval_point(0.0), [1.0, 0.0, 3.0]);
    assert_vector_eq(
        &hyperbola.eval_point(1.0),
        [4.525603580931405, 1.0861612696304874, 3.0],
    );
}

#[test]
fn derivatives_alternate_between_the_two_hyperbolic_forms() {
    let hyperbola = reference_hyperbola();

    assert_vector_eq(&hyperbola.eval_derivative(0.0, 1), [3.0, 0.0, 0.0]);
    assert_vector_eq(&hyperbola.eval_derivative(0.0, 2), [0.0, 2.0, 0.0]);
    assert_vector_eq(
        &hyperbola.eval_derivative(1.0, 3),
        [4.629241904445731, 2.3504023872876028, 0.0],
    );
    assert_vector_eq(
        &hyperbola.eval_derivative(1.0, 4),
        [3.525603580931404, 3.0861612696304874, 0.0],
    );
}

#[test]
fn tangent_and_frame_follow_the_hyperbola_orientation() {
    let hyperbola = reference_hyperbola();

    assert_vector_eq(&hyperbola.eval_tangent(0.0, true), [1.0, 0.0, 0.0]);
    assert_vector_eq(&hyperbola.eval_normal(0.0, true), [0.0, 1.0, 0.0]);
    assert_vector_eq(&hyperbola.eval_binormal(0.0, true), [0.0, 0.0, 1.0]);
}

#[test]
#[should_panic]
fn rejects_non_positive_radii() {
    Hyperbola::new(point([0.0, 0.0, 0.0]), point([0.0, 0.0, 1.0]), 0.0, 1.0);
}
