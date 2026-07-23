//! Spline evaluation
//!
//--------------------------------------------------------------------------------------------------

//{{{ crate imports
//}}}
//{{{ std imports
//}}}
//{{{ dep imports
use approx::ulps_eq;
//}}}
//--------------------------------------------------------------------------------------------------

//{{{ collection: constants
/// This is the maximum allowable order of a bspline basis. It is an arbitrary number.
pub const PMAX: usize = 8;
/// This is the tolerance with which two knots are considered equal
pub const KNOT_ULPS: u32 = 32;
//}}}
//{{{ fun: knot_lt
/// Tolerant less-thant for knots
fn knot_lt(
    u1: f64,
    u2: f64,
) -> bool {
    u1 < u2 || ulps_eq!(u1, u2, max_ulps = KNOT_ULPS)
}
//}}}
//{{{ fun: knot_gt
/// Tolerant greater-thant for knots
fn knot_gt(
    u1: f64,
    u2: f64,
) -> bool {
    u1 > u2 || ulps_eq!(u1, u2, max_ulps = KNOT_ULPS)
}
//}}}
//{{{ fun: knot_eq
/// Checks if two floating-point values are considered equal within a tolerance.
pub fn knot_eq(
    u1: f64,
    u2: f64,
) -> bool {
    ulps_eq!(u1, u2, max_ulps = KNOT_ULPS)
}
//}}}
//{{{ fun: upper_bound
/// Finds the upper bound index for a given value in a sorted array of floating-point numbers.
///
/// This function uses a binary search algorithm to find the insertion point for `value` in `arr`.
/// It takes into account the custom comparison function `knot_lt` for floating-point values.
///
/// # Parameters
///
/// - `arr`: A slice of f64 values, assumed to be sorted in ascending order.
/// - `value`: The f64 value to find the upper bound for.
///
/// # Returns
///
/// The index where `value` should be inserted to maintain sorted order.
fn upper_bound(
    arr: &[f64],
    value: f64,
) -> usize {
    let comp = |val: &f64| {
        if knot_lt(*val, value) {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    };
    arr.binary_search_by(comp).unwrap_or_else(|index| index)
}
//}}}
//{{{ fun: is_member
/// Checks if a given parameter value is within the range of the knot vector.
///
/// This function determines whether a parameter value `u` is within the range defined by the
/// knot vector, using tolerant comparisons for floating-point values.
///
/// # Parameters
///
/// - `knots`: A slice of f64 values representing the knot vector.
/// - `u`: The parameter value to check for membership.
///
/// # Returns
///
/// Returns `true` if `u` is within the range of the knot vector, `false` otherwise.
pub fn is_member(
    knots: &[f64],
    u: f64,
) -> bool {
    let umin = knots.first().unwrap();
    let umax = knots.last().unwrap();
    knot_gt(u, *umin) && knot_lt(u, *umax)
}
//}}}
//{{{ fun: find_span
/// Finds the index of the knot vector that contains the given parameter value `u`.
///
/// This function determines the index of the knot vector that contains the given parameter value `u`,
/// assuming the knot vector is sorted in ascending order. It uses a binary search algorithm to
/// efficiently locate the appropriate index.
///
/// # Parameters
///
/// - `knots`: A slice of `f64` values representing the knot vector.
/// - `u`: The parameter value to find the containing knot vector index for.
/// - `p`: The degree of the spline.
///
/// # Returns
///
/// The index of the knot vector that contains the given parameter value `u`.
pub fn find_span(
    knots: &[f64],
    u: f64,
    p: usize,
) -> usize {
    let n = knots.len() - p - 1;
    if knot_eq(u, knots[n]) {
        n - 1
    } else {
        let low = p;
        let idx = upper_bound(&knots[low..], u) + low;
        idx - 1
    }
}
//}}}
//{{{ fun: non_zero_basis
/// Finds the indices of the non-zero basis functions for the given parameter value `u`.
///
/// This function determines the indices of the non-zero basis functions for the given parameter value `u`,
/// assuming the knot vector is sorted in ascending order. It uses the `find_span` function to efficiently
/// locate the appropriate index range.
///
/// # Parameters
///
/// - `knots`: A slice of `f64` values representing the knot vector.
/// - `u`: The parameter value to find the non-zero basis function indices for.
/// - `p`: The degree of the spline.
///
/// # Returns
///
/// A tuple containing the start index, end index, and the number of non-zero basis functions.
pub fn non_zero_basis(
    knots: &[f64],
    u: f64,
    p: usize,
) -> (usize, usize, usize) {
    debug_assert!(is_member(knots, u));

    let mi = knots.len() as i32;
    let span_i = find_span(knots, u, p) as i32;
    let pi = p as i32;

    let start = std::cmp::max(0, span_i - pi) as usize;
    let end = (std::cmp::min(span_i, mi - 2 - pi) + 1) as usize;
    let num_basis = end - start;
    (start, end, num_basis)
}
//}}}
//{{{ fun: eval
/// Evaluates the B-spline basis functions for the given parameter value `u`.
///
/// This function computes the values of the B-spline basis functions for the given parameter value `u`,
/// assuming the knot vector is sorted in ascending order. It uses the `find_span` function to efficiently
/// locate the appropriate index range.
///
/// # Parameters
///
/// - `knots`: A slice of `f64` values representing the knot vector.
/// - `u`: The parameter value to evaluate the basis functions for.
/// - `p`: The degree of the spline.
/// - `shape_funs`: A mutable slice of `f64` values to store the computed basis function values.
///
/// # Panics
///
/// - If `shape_funs.len() < p + 1`, indicating the buffer is too small to hold the results.
/// - If `u` is not a member of the parameter range defined by the knot vector.
pub fn eval(
    knots: &[f64],
    u: f64,
    p: usize,
    shape_funs: &mut [f64],
) {
    debug_assert!(shape_funs.len() > p, "Buffer too small to hold results");
    debug_assert!(is_member(knots, u), "u is outside of ");

    shape_funs.fill(0.0);
    shape_funs[0] = 1.0;

    let mut left = [0.0; PMAX];
    let mut right = [0.0; PMAX];

    let i = find_span(knots, u, p);

    for j in 1..p + 1 {
        left[j - 1] = u - knots[i - p + j];
        right[j - 1] = knots[i + j] - u;
    }

    for j in 1..p + 1 {
        let mut saved = 0.0;
        for r in 0..j {
            let ri = right[r];
            let le = left[p - j + r];
            let temp = shape_funs[r] / (ri + le);
            shape_funs[r] = saved + (ri * temp);
            saved = le * temp;
        }
        shape_funs[j] = saved;
    }
}
//}}}
//{{{ fun: eval_diff
/// Evaluates the derivatives of the B-spline basis functions up to the `k`-th order.
///
/// # Parameters
///
/// - `knots`: The knot vector.
/// - `u`: The parameter value at which to evaluate the derivatives.
/// - `p`: The degree of the spline.
/// - `k`: The order of the derivative to compute.
/// - `shape_ders`: A mutable slice to store the computed basis function derivatives.
///
/// # Panics
///
/// - If `u` is not a member of the parameter range defined by the knot vector.
/// - If `p` exceeds the maximum allowable order (`PMAX`).
/// - If the length of `shape_ders` is less than `(p + 1)`, indicating the buffer is too small to hold the results.
pub fn eval_diff(
    knots: &[f64],
    u: f64,
    p: usize,
    k: usize,
    shape_ders: &mut [f64],
) {
    debug_assert!(is_member(knots, u), "u is not member of parameter range");
    debug_assert!(p <= PMAX, "order exceeds max allowable");
    // debug_assert!(k <= p, "order of derivative exceeds order");
    debug_assert!(shape_ders.len() > p, "Derivative buffer too small");

    let i = find_span(knots, u, p);
    let p_k = p - k;
    shape_ders.fill(0.0);
    if k <= p {
        /* doc: compute shape functions
        This code-block initialises the first p_k+1 elements of
        shape_ders to the shape functions:

            N_{i-p_k, p_k}^{(0)} ... N_{i-p_k, p_k}^{(0)}

        This is the first phase of the cox-deboor recursion.
        */
        shape_ders[0] = 1.0;
        let mut left = [0.0; PMAX];
        let mut right = [0.0; PMAX];
        for j in 1..p_k + 1 {
            left[j - 1] = u - knots[i - p_k + j];
            right[j - 1] = knots[i + j] - u;
        }
        for j in 1..p_k + 1 {
            let mut saved = 0.0;
            for r in 0..j {
                let ri = right[r];
                let le = left[p_k - j + r];
                let temp = shape_ders[r] / (ri + le);
                shape_ders[r] = saved + (ri * temp);
                saved = le * temp;
            }
            shape_ders[j] = saved;
        }
    }

    for k2 in 1..k + 1
    // loop over derivative
    {
        let p_k2 = p_k + k2;
        let p_k2_f = p_k2 as f64;
        let mut saved = 0.0;
        for r in 0..p_k2 {
            let denom = knots[i + 1 + r] - knots[i + 1 + r - p_k2];
            let temp = (p_k2_f / denom) * shape_ders[r];
            shape_ders[r] = saved - temp;
            saved = temp;
        }
        shape_ders[p_k2] = saved;
    }
}
//}}}
//{{{ fun: eval_diff_all
pub fn eval_diff_all(
    knots: &[f64],
    u: f64,
    p: usize,
    k_in: usize,
    shape_ders: &mut [f64],
) {
    debug_assert!(is_member(knots, u), "u is not member of parameter range");
    debug_assert!(p < PMAX, "order exceeds max allowable");
    // debug_assert!(k <= p, "order of derivative exceeds order");
    debug_assert!(
        shape_ders.len() >= (k_in + 1) * (p + 1),
        "Derivative buffer too small"
    );

    let k = p.min(k_in);
    shape_ders.fill(0.0);
    let row_size = p + 1;
    eval(knots, u, p, &mut shape_ders[..row_size]);
    for k2 in 1..k + 1 {
        let start = k2 * row_size;
        eval_diff(knots, u, p, k2, &mut shape_ders[start..start + row_size]);
    }
}
//}}}
//{{{ fun: multiplicites
pub fn multiplicites(knots: &[f64]) -> Vec<(f64, usize)> {
    let mut out = Vec::with_capacity(knots.len());

    let mut cur_knot = knots[0];
    let mut cur_mult = 0;

    for &knot in knots.iter().skip(1) {
        cur_mult += 1;
        if !knot_eq(knot, cur_knot) {
            out.push((cur_knot, cur_mult));
            cur_knot = knot;
            cur_mult = 0;
        }
    }

    out.push((cur_knot, cur_mult + 1));
    out
}
//}}}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const U0: f64 = 0.6060606060606061;
    const U1: f64 = 2.4242424242424243;
    const U2: f64 = 4.545454545454546;
    const TOLERANCE: f64 = 1.0e-14;

    const KNOTS_P0: &[f64] = &[0.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0];
    const KNOTS_P1: &[f64] = &[0.0, 0.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0, 6.0];
    const KNOTS_P2: &[f64] = &[0.0, 0.0, 0.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0, 6.0, 6.0];
    const KNOTS_P3: &[f64] = &[
        0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0, 6.0, 6.0, 6.0,
    ];
    const KNOTS_P4: &[f64] = &[
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 4.0, 4.0, 5.0, 6.0, 6.0, 6.0, 6.0, 6.0,
    ];

    fn knots(p: usize) -> &'static [f64] {
        match p {
            0 => KNOTS_P0,
            1 => KNOTS_P1,
            2 => KNOTS_P2,
            3 => KNOTS_P3,
            4 => KNOTS_P4,
            _ => unreachable!(),
        }
    }

    fn assert_slice_eq(
        actual: &[f64],
        expected: &[f64],
    ) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            assert_abs_diff_eq!(actual, expected, epsilon = TOLERANCE);
        }
    }

    #[test]
    fn knot_comparison_and_membership_match_reference_range() {
        assert!(knot_eq(KNOTS_P4[0], KNOTS_P4[4]));
        assert!(!knot_eq(KNOTS_P4[4], KNOTS_P4[5]));

        assert!(is_member(KNOTS_P4, 0.0));
        assert!(is_member(KNOTS_P4, U1));
        assert!(is_member(KNOTS_P4, 6.0));
        assert!(!is_member(KNOTS_P4, -0.1));
        assert!(!is_member(KNOTS_P4, 6.1));
    }

    #[test]
    fn spans_and_non_zero_ranges_match_reference_data() {
        let reference_spans = [[0, 2, 4], [1, 3, 5], [2, 4, 6], [3, 5, 7], [4, 6, 8]];

        for (p, spans) in reference_spans.into_iter().enumerate() {
            for (u, expected_span) in [U0, U1, U2].into_iter().zip(spans) {
                assert_eq!(find_span(knots(p), u, p), expected_span);
                assert_eq!(
                    non_zero_basis(knots(p), u, p),
                    (expected_span - p, expected_span + 1, p + 1),
                );
            }
        }
    }

    #[test]
    fn basis_functions_match_reference_data() {
        let cases: &[(usize, f64, &[f64])] = &[
            (0, U0, &[1.0]),
            (0, U1, &[1.0]),
            (0, U2, &[1.0]),
            (1, U0, &[0.3939393939393939, 0.6060606060606061]),
            (1, U1, &[0.7878787878787878, 0.21212121212121215]),
            (1, U2, &[0.45454545454545414, 0.5454545454545459]),
            (
                2,
                U0,
                &[0.155188246097337, 0.6611570247933884, 0.18365472910927458],
            ),
            (
                2,
                U1,
                &[0.41383532292623193, 0.5411692684419958, 0.04499540863177228],
            ),
            (
                2,
                U2,
                &[0.20661157024793353, 0.6446280991735538, 0.14876033057851262],
            ),
            (
                3,
                U0,
                &[
                    0.06113476361410245,
                    0.5548598937028689,
                    0.3561788685755628,
                    0.027826474107465846,
                ],
            ),
            (
                3,
                U1,
                &[
                    0.16302603630427318,
                    0.5350598114601787,
                    0.2955511651563075,
                    0.006362987079240525,
                ],
            ),
            (
                3,
                U2,
                &[
                    0.03130478337089899,
                    0.644127222639619,
                    0.28399699474079665,
                    0.04057099924868529,
                ],
            ),
            (
                4,
                U0,
                &[
                    0.02408339172676763,
                    0.42377190386206165,
                    0.4703517350649833,
                    0.07757683690566236,
                    0.004216132440525129,
                ],
            ),
            (
                4,
                U1,
                &[
                    0.06422237793804701,
                    0.30958479621417534,
                    0.5145957117810639,
                    0.11092225180073362,
                    0.0006748622659800558,
                ],
            ),
            (
                4,
                U2,
                &[
                    0.0035573617466930635,
                    0.2619755025840673,
                    0.6164424105821551,
                    0.10695990711017025,
                    0.011064817976914179,
                ],
            ),
        ];

        for &(p, u, expected) in cases {
            let mut actual = [f64::NAN; PMAX + 1];
            eval(knots(p), u, p, &mut actual);
            assert_slice_eq(&actual[..p + 1], expected);
        }
    }

    #[test]
    fn basis_derivatives_match_reference_data() {
        // Each consecutive p + 1 values is one derivative order at U1.
        let cases: &[(usize, &[f64])] = &[
            (0, &[1.0]),
            (1, &[0.7878787878787878, 0.21212121212121215, -0.5, 0.5]),
            (
                2,
                &[
                    0.41383532292623193,
                    0.5411692684419958,
                    0.04499540863177228,
                    -0.5252525252525252,
                    0.31313131313131304,
                    0.21212121212121215,
                    0.3333333333333333,
                    -0.8333333333333333,
                    0.5,
                ],
            ),
            (
                3,
                &[
                    0.16302603630427318,
                    0.5350598114601787,
                    0.2955511651563075,
                    0.006362987079240525,
                    -0.31037649219467395,
                    -0.23079277624732178,
                    0.49617385981022344,
                    0.04499540863177228,
                    0.3939393939393939,
                    -0.7070707070707068,
                    0.10101010101010088,
                    0.21212121212121215,
                    -0.25,
                    1.083333333333333,
                    -1.3333333333333333,
                    0.5,
                ],
            ),
            (
                4,
                &[
                    0.06422237793804701,
                    0.30958479621417534,
                    0.5145957117810639,
                    0.11092225180073362,
                    0.0006748622659800558,
                    -0.16302603630427318,
                    -0.3720337751559055,
                    0.23950864630387114,
                    0.289188178077067,
                    0.006362987079240525,
                    0.31037649219467395,
                    -0.07958371594735217,
                    -0.7269666360575453,
                    0.4511784511784511,
                    0.04499540863177228,
                    -0.3939393939393939,
                    1.1010101010101008,
                    -0.8080808080808077,
                    -0.11111111111111127,
                    0.21212121212121215,
                    0.25,
                    -1.3333333333333333,
                    2.416666666666667,
                    -1.8333333333333335,
                    0.5,
                ],
            ),
        ];

        for &(p, expected) in cases {
            let row_size = p + 1;
            for derivative in 0..=p {
                let mut actual = [f64::NAN; PMAX + 1];
                eval_diff(knots(p), U1, p, derivative, &mut actual);
                let start = derivative * row_size;
                assert_slice_eq(&actual[..row_size], &expected[start..start + row_size]);
            }

            let mut actual_all = vec![f64::NAN; row_size * row_size];
            eval_diff_all(knots(p), U1, p, p, &mut actual_all);
            assert_slice_eq(&actual_all, expected);
        }
    }

    #[test]
    fn knot_multiplicities_match_reference_knot_vector() {
        assert_eq!(
            multiplicites(KNOTS_P4),
            vec![(0.0, 5), (1.0, 1), (2.0, 1), (4.0, 2), (5.0, 1), (6.0, 5)],
        );
    }
}
//}}}
