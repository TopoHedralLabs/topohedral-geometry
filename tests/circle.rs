//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Circle, Curve};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::assert_abs_diff_eq;
use topohedral_linalg::{DVector, VecType::Col, VectorOps};
//}}}
//--------------------------------------------------------------------------------------------------

const TOLERANCE: f64 = 1.0e-12;

fn vector(values: &[f64]) -> DVector<f64> {
    DVector::from_slice_vec(values, values.len(), Col)
}

fn assert_vector_eq(
    actual: &DVector<f64>,
    expected: &[f64],
) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = TOLERANCE);
    }
}

fn reference_circle() -> Circle {
    Circle::new(vector(&[1.0, -2.0, 3.0]), vector(&[0.0, 0.0, 1.0]), 2.0)
}

#[test]
fn dimension_and_cardinal_points() {
    let circle = reference_circle();

    assert_eq!(circle.dimension(), 3);
    assert_vector_eq(&circle.eval_point(0.0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(
        &circle.eval_point(std::f64::consts::FRAC_PI_2),
        &[3.0, -2.0, 3.0],
    );
    assert_vector_eq(&circle.eval_point(std::f64::consts::PI), &[1.0, -4.0, 3.0]);
    assert_vector_eq(
        &circle.eval_point(3.0 * std::f64::consts::FRAC_PI_2),
        &[-1.0, -2.0, 3.0],
    );
}

#[test]
fn derivatives_repeat_every_four_orders() {
    let circle = reference_circle();

    assert_vector_eq(&circle.eval_derivative(0.0, 0), &[1.0, 0.0, 3.0]);
    assert_vector_eq(&circle.eval_derivative(0.0, 1), &[2.0, 0.0, 0.0]);
    assert_vector_eq(&circle.eval_derivative(0.0, 2), &[0.0, -2.0, 0.0]);
    assert_vector_eq(&circle.eval_derivative(0.0, 3), &[-2.0, 0.0, 0.0]);
    assert_vector_eq(&circle.eval_derivative(0.0, 4), &[0.0, 2.0, 0.0]);
    assert_vector_eq(&circle.eval_derivative(0.0, 5), &[2.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_are_written_to_the_output_slice() {
    let circle = reference_circle();
    let mut derivatives = vec![vector(&[99.0, 99.0, 99.0]); 5];

    circle.eval_all_derivatives(0.0, 4, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[1.0, 0.0, 3.0]);
    assert_vector_eq(&derivatives[1], &[2.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[2], &[0.0, -2.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[-2.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[4], &[0.0, 2.0, 0.0]);
}

#[test]
fn tangent_and_frame_follow_the_circle_orientation() {
    let circle = reference_circle();

    assert_vector_eq(&circle.eval_tangent(0.0, true), &[1.0, 0.0, 0.0]);
    assert_vector_eq(&circle.eval_normal(0.0, true), &[0.0, -1.0, 0.0]);
    assert_vector_eq(&circle.eval_binormal(0.0, true), &[0.0, 0.0, -1.0]);
}
