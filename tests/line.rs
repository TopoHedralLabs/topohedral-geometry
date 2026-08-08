//! Short Description of module
//!
//! Longer description of module
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
use topohedral_geometry::{Curve, Line, Vec3};
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::assert_abs_diff_eq;
use topohedral_linalg::VectorOps;
//}}}
//--------------------------------------------------------------------------------------------------

fn vector(values: [f64; 3]) -> Vec3 {
    Vec3::from_col_slice(&values)
}

fn assert_vector_eq(
    actual: &Vec3,
    expected: &[f64],
) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_abs_diff_eq!(actual, expected, epsilon = 1.0e-12);
    }
}

fn reference_line() -> Line {
    Line::new(
        vector([1.0, -2.0, 3.0]),
        vector([
            2.0 / 21.0_f64.sqrt(),
            4.0 / 21.0_f64.sqrt(),
            -1.0 / 21.0_f64.sqrt(),
        ]),
    )
}

#[test]
fn dimension_and_point_evaluation() {
    let line = reference_line();

    assert_eq!(line.dimension(), 3);
    assert_vector_eq(&line.eval_point(0.0), &[1.0, -2.0, 3.0]);
    assert_vector_eq(
        &line.eval_point(0.5),
        &[
            1.0 + 1.0 / 21.0_f64.sqrt(),
            -2.0 + 2.0 / 21.0_f64.sqrt(),
            3.0 - 0.5 / 21.0_f64.sqrt(),
        ],
    );
    assert_vector_eq(
        &line.eval_point(-2.0),
        &[
            1.0 - 4.0 / 21.0_f64.sqrt(),
            -2.0 - 8.0 / 21.0_f64.sqrt(),
            3.0 + 2.0 / 21.0_f64.sqrt(),
        ],
    );
}

#[test]
fn derivatives_are_constant_or_zero() {
    let line = reference_line();

    assert_vector_eq(
        &line.eval_derivative(0.25, 0),
        &[
            1.0 + 0.5 / 21.0_f64.sqrt(),
            -2.0 + 1.0 / 21.0_f64.sqrt(),
            3.0 - 0.25 / 21.0_f64.sqrt(),
        ],
    );
    assert_vector_eq(
        &line.eval_derivative(0.25, 1),
        &[
            2.0 / 21.0_f64.sqrt(),
            4.0 / 21.0_f64.sqrt(),
            -1.0 / 21.0_f64.sqrt(),
        ],
    );
    assert_vector_eq(&line.eval_derivative(0.25, 2), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_derivative(0.25, 4), &[0.0, 0.0, 0.0]);
}

#[test]
fn all_derivatives_are_written_to_the_output_slice() {
    let line = reference_line();
    let mut derivatives = vec![vector([99.0, 99.0, 99.0]); 4];

    line.eval_all_derivatives(0.5, 3, &mut derivatives);

    assert_vector_eq(
        &derivatives[0],
        &[
            1.0 + 1.0 / 21.0_f64.sqrt(),
            -2.0 + 2.0 / 21.0_f64.sqrt(),
            3.0 - 0.5 / 21.0_f64.sqrt(),
        ],
    );
    assert_vector_eq(
        &derivatives[1],
        &[
            2.0 / 21.0_f64.sqrt(),
            4.0 / 21.0_f64.sqrt(),
            -1.0 / 21.0_f64.sqrt(),
        ],
    );
    assert_vector_eq(&derivatives[2], &[0.0, 0.0, 0.0]);
    assert_vector_eq(&derivatives[3], &[0.0, 0.0, 0.0]);
}

#[test]
fn tangent_is_already_normalized() {
    let line = Line::new(vector([0.0, 0.0, 0.0]), vector([0.6, 0.8, 0.0]));

    assert_vector_eq(&line.eval_tangent(0.75, false), &[0.6, 0.8, 0.0]);
    assert_vector_eq(&line.eval_tangent(0.75, true), &[0.6, 0.8, 0.0]);
}

#[test]
#[should_panic]
fn rejects_non_unit_direction() {
    Line::new(vector([0.0, 0.0, 0.0]), vector([3.0, 4.0, 0.0]));
}

#[test]
fn straight_line_has_zero_three_dimensional_frame() {
    let line = reference_line();

    assert_vector_eq(&line.eval_normal(0.5, false), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_normal(0.5, true), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_binormal(0.5, false), &[0.0, 0.0, 0.0]);
    assert_vector_eq(&line.eval_binormal(0.5, true), &[0.0, 0.0, 0.0]);
}
