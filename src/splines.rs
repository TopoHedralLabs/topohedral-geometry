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
pub const PMAX_SQR: usize = 64;
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
//{{{ fun: knot_multiplicity
pub fn knot_multiplicity(
    knots: &[f64],
    knot: f64,
) -> usize {
    let mut knot_mult = 0;
    for knot_i in knots {
        if knot_eq(*knot_i, knot) {
            knot_mult += 1;
        }
    }
    knot_mult
}
//}}}

//-------------------------------------------------------------------------------------------------
//{{{ mod: tests
#[cfg(test)]
mod tests {
    use super::*;

    const U0: f64 = 0.6060606060606061;
    const U1: f64 = 2.4242424242424243;
    const U2: f64 = 4.545454545454546;

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
}
//}}}
