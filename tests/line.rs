//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Curve, Line};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::assert_abs_diff_eq;
use topohedral_linalg::{DVector, VecType::Col, VectorOps};
//}}}
//--------------------------------------------------------------------------------------------------

fn vector(values: &[f64]) -> DVector<f64> {
    DVector::from_slice_vec(values, values.len(), Col)
}

fn assert_vector_eq(
    actual: &DVector<f64>,
    expected: &[f64],
) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1.0e-12);
    }
}

fn reference_line() -> Line {
    Line::new(vector(&[1.0, -2.0, 3.0]), vector(&[2.0, 4.0, -1.0]))
}

#[test]
fn dimension_and_point_evaluation() {
    let line = reference_line();

    assert_eq!(line.dimension(), 3);
    assert_vector_eq(&line.eval_point(0.0), &[1.0, -2.0, 3.0]);
    assert_vector_eq(&line.eval_point(0.5), &[2.0, 0.0, 2.5]);
    assert_vector_eq(&line.eval_point(-2.0), &[-3.0, -10.0, 5.0]);
}

#[test]
fn derivatives_are_constant_or_zero() {
    let line = reference_line();

    assert_vector_eq(&line.eval_derivative(0.25, 0), &[1.5, -1.0, 2.75]);
    assert_vector_eq(&line.eval_derivative(0.25, 1), &[2.0, 4.0, -1.0]);
    assert_vector_eq(&line.eval_derivative(0.25, 2), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_derivative(0.25, 4), &[0.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_are_written_to_the_output_slice() {
    let line = reference_line();
    let mut derivatives = vec![vector(&[99.0, 99.0, 99.0]); 4];

    line.eval_all_derivatives(0.5, 3, &mut derivatives);

    assert_vector_eq(&derivatives[0], &[2.0, 0.0, 2.5]);
    assert_vector_eq(&derivatives[1], &[2.0, 4.0, -1.0]);
    assert_vector_eq(&derivatives[2], &[0.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[0.0, 0.0, 0.0]);
}

#[test]
fn tangent_can_be_normalized() {
    let line = Line::new(vector(&[0.0, 0.0]), vector(&[3.0, 4.0]));

    assert_vector_eq(&line.eval_tangent(0.75, false), &[3.0, 4.0]);
    assert_vector_eq(&line.eval_tangent(0.75, true), &[0.6, 0.8]);
}

#[test]
fn straight_line_has_zero_three_dimensional_frame() {
    let line = reference_line();

    assert_vector_eq(&line.eval_normal(0.5, false), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_normal(0.5, true), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_binormal(0.5, false), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_binormal(0.5, true), &[0.0, 0.0, 0.0]);
}
