//! Tests for tensor-product B-spline surface evaluation.

use approx::assert_relative_eq;
use topohedral_geometry::{BSurface, Surface, Vec3};

fn point(values: [f64; 3]) -> Vec3 {
    Vec3::from_col_slice(&values)
}

fn assert_vector_eq(
    actual: &Vec3,
    expected: [f64; 3],
) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert_relative_eq!(*actual, expected, epsilon = 1.0e-12, max_relative = 1.0e-12);
    }
}

fn polynomial_surface() -> BSurface {
    // This biquadratic-linear surface is exactly S(u, v) = (u^2, v, uv).
    BSurface::new(
        2,
        1,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        vec![
            point([0.0, 0.0, 0.0]),
            point([0.0, 1.0, 0.0]),
            point([0.0, 0.0, 0.0]),
            point([0.0, 1.0, 0.5]),
            point([1.0, 0.0, 0.0]),
            point([1.0, 1.0, 1.0]),
        ],
        // A common non-unit weight must not scale points or derivatives.
        vec![2.0; 6],
    )
}

fn rational_surface() -> BSurface {
    // A degree-elevated representation of the rational bilinear surface with weighted numerator
    // (u(3 + v), 2v(1 + u), 4uv) and denominator 1 + 2u + v.
    BSurface::new(
        2,
        2,
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![
            point([0.0, 0.0, 0.0]),
            point([0.0, 2.0 / 3.0, 0.0]),
            point([0.0, 1.0, 0.0]),
            point([0.75, 0.0, 0.0]),
            point([0.7, 0.6, 0.4]),
            point([2.0 / 3.0, 1.0, 2.0 / 3.0]),
            point([1.0, 0.0, 0.0]),
            point([1.0, 4.0 / 7.0, 4.0 / 7.0]),
            point([1.0, 1.0, 1.0]),
        ],
        vec![1.0, 1.5, 2.0, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0],
    )
}

#[test]
fn eval_point_matches_polynomial_and_rational_reference_values() {
    assert_vector_eq(
        &polynomial_surface().eval_point(0.25, 0.4),
        [0.0625, 0.4, 0.1],
    );

    assert_vector_eq(
        &rational_surface().eval_point(0.25, 0.5),
        [0.4375, 0.625, 0.25],
    );
}

#[test]
fn eval_derivatives_match_polynomial_reference_values() {
    let surface = polynomial_surface();

    assert_vector_eq(&surface.eval_derivative(0.25, 0.4, 1, 0), [0.5, 0.0, 0.4]);
    assert_vector_eq(&surface.eval_derivative(0.25, 0.4, 2, 0), [2.0, 0.0, 0.0]);
    assert_vector_eq(&surface.eval_derivative(0.25, 0.4, 0, 1), [0.0, 1.0, 0.25]);
    assert_vector_eq(&surface.eval_derivative(0.25, 0.4, 1, 1), [0.0, 0.0, 1.0]);

    let mut derivatives = vec![Vec3::zeros(); 6];
    surface.eval_all_derivatives(0.25, 0.4, 2, 1, &mut derivatives);
    for (actual, expected) in derivatives.iter().zip([
        [0.0625, 0.4, 0.1],
        [0.5, 0.0, 0.4],
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.25],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],
    ]) {
        assert_vector_eq(actual, expected);
    }
}

#[test]
fn eval_derivatives_match_rational_reference_values() {
    let surface = rational_surface();

    assert_vector_eq(
        &surface.eval_derivative(0.25, 0.5, 2, 0),
        [-2.625, 0.25, -1.5],
    );
    assert_vector_eq(
        &surface.eval_derivative(0.25, 0.5, 0, 2),
        [0.09375, -0.9375, -0.375],
    );

    let mut derivatives = vec![Vec3::zeros(); 4];
    surface.eval_all_derivatives(0.25, 0.5, 1, 1, &mut derivatives);

    for (actual, expected) in derivatives.iter().zip([
        [0.4375, 0.625, 0.25],
        [1.3125, -0.125, 0.75],
        [-0.09375, 0.9375, 0.375],
        [-0.0625, 0.125, 1.25],
    ]) {
        assert_vector_eq(actual, expected);
    }
}
