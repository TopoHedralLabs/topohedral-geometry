//! Tests for the ring torus surface.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use approx::assert_abs_diff_eq;
use topohedral_geometry::{Surface, Torus, Vec3};
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

fn reference_torus() -> Torus {
    Torus::new(vector([1.0, -2.0, 3.0]), vector([0.0, 0.0, 1.0]), 3.0, 1.0)
}

#[test]
fn cardinal_points() {
    let torus = reference_torus();

    assert_vector_eq(&torus.eval_point(0.0, 0.0), &[1.0, 2.0, 3.0]);
    assert_vector_eq(&torus.eval_point(FRAC_PI_2, 0.0), &[5.0, -2.0, 3.0]);
    assert_vector_eq(&torus.eval_point(PI, 0.0), &[1.0, -6.0, 3.0]);
    assert_vector_eq(&torus.eval_point(0.0, PI), &[1.0, 0.0, 3.0]);
    assert_vector_eq(&torus.eval_point(0.0, FRAC_PI_2), &[1.0, 1.0, 4.0]);
    assert_vector_eq(&torus.eval_point(0.0, -FRAC_PI_2), &[1.0, 1.0, 2.0]);
}

#[test]
fn default_uses_standard_axes_and_radii() {
    let torus = Torus::default();

    assert_vector_eq(&torus.eval_point(0.0, 0.0), &[1.1, 0.0, 0.0]);
    assert_vector_eq(&torus.eval_point(FRAC_PI_2, 0.0), &[0.0, 1.1, 0.0]);
    assert_vector_eq(&torus.eval_point(0.0, FRAC_PI_2), &[1.0, 0.0, 0.1]);
}

#[test]
fn zero_order_derivative_matches_point() {
    let torus = reference_torus();

    assert_vector_eq(
        &torus.eval_derivative(0.7, -0.4, 0, 0),
        torus.eval_point(0.7, -0.4).as_ref(),
    );
}

#[test]
fn angular_derivatives_repeat_every_four_orders() {
    let torus = reference_torus();
    let expected_u = [
        [4.0, 0.0, 0.0],
        [0.0, -4.0, 0.0],
        [-4.0, 0.0, 0.0],
        [0.0, 4.0, 0.0],
    ];
    let expected_v = [
        [0.0, 0.0, 1.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 1.0, 0.0],
    ];
    for order in 1..=9 {
        assert_vector_eq(
            &torus.eval_derivative(0.0, 0.0, order, 0),
            &expected_u[(order - 1) % 4],
        );
        assert_vector_eq(
            &torus.eval_derivative(0.0, 0.0, 0, order),
            &expected_v[(order - 1) % 4],
        );
    }
}

#[test]
fn mixed_derivatives_are_correct() {
    let torus = reference_torus();

    assert_vector_eq(&torus.eval_derivative(0.0, 0.0, 1, 1), &[0.0, 0.0, 0.0]);
    assert_vector_eq(
        &torus.eval_derivative(0.0, FRAC_PI_2, 1, 1),
        &[-1.0, 0.0, 0.0],
    );
    assert_vector_eq(
        &torus.eval_derivative(0.0, FRAC_PI_2, 2, 1),
        &[0.0, 1.0, 0.0],
    );
    assert_vector_eq(&torus.eval_derivative(0.0, 0.0, 1, 2), &[-1.0, 0.0, 0.0]);
    assert_vector_eq(&torus.eval_derivative(0.0, 0.0, 2, 2), &[0.0, 1.0, 0.0]);
    assert_vector_eq(&torus.eval_derivative(0.0, 0.0, 4, 4), &[0.0, 1.0, 0.0]);
}

#[test]
fn derivatives_match_finite_differences_away_from_cardinal_points() {
    let torus = reference_torus();
    let (u, v, h) = (0.7, -0.4, 1.0e-5);
    for (m, n) in [(0, 0), (1, 0), (0, 1), (1, 1), (3, 2)] {
        let finite_u = (torus.eval_derivative(u + h, v, m, n)
            - torus.eval_derivative(u - h, v, m, n))
            / (2.0 * h);
        let finite_v = (torus.eval_derivative(u, v + h, m, n)
            - torus.eval_derivative(u, v - h, m, n))
            / (2.0 * h);
        let derivative_u = torus.eval_derivative(u, v, m + 1, n);
        let derivative_v = torus.eval_derivative(u, v, m, n + 1);
        for i in 0..3 {
            assert_abs_diff_eq!(finite_u[i], derivative_u[i], epsilon = 1.0e-9);
            assert_abs_diff_eq!(finite_v[i], derivative_v[i], epsilon = 1.0e-9);
        }
    }
}

#[test]
fn all_derivatives_use_column_major_parameter_order() {
    let torus = reference_torus();
    let mut derivatives = vec![vector([99.0, 99.0, 99.0]); 6];

    torus.eval_all_derivatives(0.0, FRAC_PI_2, 2, 1, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[1.0, 1.0, 4.0]);
    assert_vector_eq(&derivatives[1], &[3.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[0.0, -3.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[0.0, -1.0, 0.0]);
    assert_vector_eq(&derivatives[4], &[-1.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[5], &[0.0, 1.0, 0.0]);
}

#[test]
fn points_and_derivatives_are_periodic_in_both_parameters() {
    let torus = reference_torus();
    for m in 0..=4 {
        for n in 0..=4 {
            let expected = torus.eval_derivative(0.7, -0.4, m, n);
            assert_vector_eq(
                &torus.eval_derivative(0.7 + TAU, -0.4, m, n),
                expected.as_ref(),
            );
            assert_vector_eq(
                &torus.eval_derivative(0.7, -0.4 + TAU, m, n),
                expected.as_ref(),
            );
        }
    }
}

#[test]
fn normals_follow_parameter_orientation() {
    let torus = reference_torus();

    assert_vector_eq(&torus.eval_normal(0.0, 0.0, false), &[0.0, -4.0, 0.0]);
    assert_vector_eq(&torus.eval_normal(0.0, 0.0, true), &[0.0, -1.0, 0.0]);
    assert_vector_eq(&torus.eval_normal(0.0, FRAC_PI_2, false), &[0.0, 0.0, -3.0]);
    assert_vector_eq(&torus.eval_normal(0.0, PI, true), &[0.0, 1.0, 0.0]);
}

#[test]
fn rotated_axis_preserves_torus_geometry() {
    let center = vector([1.0, -2.0, 3.0]);
    let axis = vector([1.0, 2.0, 3.0]).normalize();
    let torus = Torus::new(center, axis, 3.0, 1.0);

    for (u, v) in [(0.0, 0.0), (0.7, -0.4), (1.2, 2.3), (PI, PI)] {
        let offset = torus.eval_point(u, v) - center;
        let height = offset.dot(&axis);
        let radial = offset - height * axis;
        assert_abs_diff_eq!(height, v.sin(), epsilon = TOLERANCE);
        assert_abs_diff_eq!(radial.norm(), 3.0 + v.cos(), epsilon = TOLERANCE);
        assert_abs_diff_eq!(
            (radial.norm() - 3.0).powi(2) + height.powi(2),
            1.0,
            epsilon = TOLERANCE
        );
    }
}

#[test]
#[should_panic]
fn rejects_zero_major_radius() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), 0.0, 1.0);
}

#[test]
#[should_panic]
fn rejects_negative_major_radius() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), -3.0, 1.0);
}

#[test]
#[should_panic]
fn rejects_zero_minor_radius() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), 3.0, 0.0);
}

#[test]
#[should_panic]
fn rejects_negative_minor_radius() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), 3.0, -1.0);
}

#[test]
#[should_panic]
fn rejects_equal_radii() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), 1.0, 1.0);
}

#[test]
#[should_panic]
fn rejects_major_radius_smaller_than_minor_radius() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 1.0]), 1.0, 2.0);
}

#[test]
#[should_panic]
fn rejects_non_unit_axis() {
    Torus::new(Vec3::zeros(), vector([0.0, 0.0, 2.0]), 3.0, 1.0);
}

#[test]
#[should_panic]
fn rejects_zero_axis() {
    Torus::new(Vec3::zeros(), Vec3::zeros(), 3.0, 1.0);
}
