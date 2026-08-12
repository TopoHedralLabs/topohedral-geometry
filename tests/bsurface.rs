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
#[should_panic(expected = "order_u too high")]
fn new_rejects_u_order_above_maximum() {
    BSurface::new(usize::MAX, 0, vec![], vec![], vec![], vec![]);
}

#[test]
#[should_panic(expected = "order_v too high")]
fn new_rejects_v_order_above_maximum() {
    BSurface::new(0, usize::MAX, vec![], vec![], vec![], vec![]);
}

#[test]
#[should_panic(expected = "knots u not sorted")]
fn new_rejects_unsorted_u_knots() {
    BSurface::new(
        0,
        0,
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![point([0.0, 0.0, 0.0])],
        vec![1.0],
    );
}

#[test]
#[should_panic(expected = "knots v not sorted")]
fn new_rejects_unsorted_v_knots() {
    BSurface::new(
        0,
        0,
        vec![0.0, 1.0],
        vec![1.0, 0.0],
        vec![point([0.0, 0.0, 0.0])],
        vec![1.0],
    );
}

#[test]
#[should_panic(expected = "Weights must be positive")]
fn new_rejects_negative_weights() {
    BSurface::new(
        0,
        0,
        vec![0.0, 1.0],
        vec![0.0, 1.0],
        vec![point([0.0, 0.0, 0.0])],
        vec![-1.0],
    );
}

#[test]
#[should_panic(expected = "Points and weights must have same size")]
fn new_rejects_mismatched_point_and_weight_counts() {
    BSurface::new(
        0,
        0,
        vec![0.0, 1.0],
        vec![0.0, 1.0],
        vec![point([0.0, 0.0, 0.0])],
        vec![1.0, 1.0],
    );
}

#[test]
#[should_panic(expected = "Knots have wrong size")]
fn new_rejects_inconsistent_control_net_size() {
    BSurface::new(
        1,
        1,
        vec![0.0, 0.0, 1.0, 1.0],
        vec![0.0, 0.0, 1.0, 1.0],
        vec![point([0.0, 0.0, 0.0])],
        vec![1.0],
    );
}

#[test]
#[should_panic(expected = "Knots u have wrong size")]
fn new_rejects_too_short_u_knot_vector() {
    BSurface::new(2, 0, vec![0.0, 1.0], vec![0.0, 1.0], vec![], vec![]);
}

#[test]
#[should_panic(expected = "Knots v have wrong size")]
fn new_rejects_too_short_v_knot_vector() {
    BSurface::new(0, 2, vec![0.0, 1.0], vec![0.0, 1.0], vec![], vec![]);
}

#[test]
#[should_panic(expected = "At least one point is required")]
fn new_rejects_empty_control_net() {
    BSurface::new(1, 1, vec![0.0, 1.0], vec![0.0, 1.0], vec![], vec![]);
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
