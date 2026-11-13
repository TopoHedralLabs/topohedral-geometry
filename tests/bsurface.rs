//! Tests for tensor-product B-spline surface evaluation.

use approx::assert_relative_eq;
use topohedral_geometry::{knot_multiplicities, BSurface, Curve, Surface, Vec3};

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

//{{{ fun: reference_surface
/// A rational bicubic-by-biquadratic surface with interior knots in both directions.
///
/// The control net is a deterministic ripple over a 7-by-5 grid, with non-uniform weights so that
/// the surface is genuinely rational and all four of the rational code paths are exercised.
fn reference_surface() -> BSurface {
    let order_u = 3;
    let order_v = 2;
    let knots_u = vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0];

    let num_points_u = knots_u.len() - order_u - 1;
    let num_points_v = knots_v.len() - order_v - 1;

    let mut points = Vec::with_capacity(num_points_u * num_points_v);
    let mut weights = Vec::with_capacity(num_points_u * num_points_v);
    for i in 0..num_points_u {
        for j in 0..num_points_v {
            let x = i as f64;
            let y = 1.5 * j as f64;
            let z = ((i * j) % 5) as f64 - 2.0;
            points.push(point([x, y, z]));
            weights.push(1.0 + 0.25 * ((i + 2 * j) % 4) as f64);
        }
    }

    BSurface::new(order_u, order_v, knots_u, knots_v, points, weights)
}
//}}}
//{{{ fun: assert_vectors_close
fn assert_vectors_close(
    actual: &Vec3,
    expected: &Vec3,
) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert_relative_eq!(*actual, expected, epsilon = 1.0e-10, max_relative = 1.0e-10);
    }
}
//}}}
//{{{ collection: iso-curve extraction
//{{{ test: extract_isocurve_at_u_matches_surface
#[test]
fn extract_isocurve_at_u_matches_surface() {
    let surface = reference_surface();

    for u in [0.0, 0.25, 0.42, 0.5, 1.0] {
        let curve = surface.extract_isocurve_at_u(u);

        assert_eq!(curve.order(), surface.order_v());
        assert_eq!(curve.knots(), surface.knots_v());

        for v in [0.0, 0.3, 0.47, 0.6, 1.0] {
            assert_vectors_close(&curve.eval_point(v), &surface.eval_point(u, v));
        }
    }
}
//}}}
//{{{ test: extract_isocurve_at_v_matches_surface
#[test]
fn extract_isocurve_at_v_matches_surface() {
    let surface = reference_surface();

    for v in [0.0, 0.3, 0.47, 0.6, 1.0] {
        let curve = surface.extract_isocurve_at_v(v);

        assert_eq!(curve.order(), surface.order_u());
        assert_eq!(curve.knots(), surface.knots_u());

        for u in [0.0, 0.25, 0.42, 0.5, 1.0] {
            assert_vectors_close(&curve.eval_point(u), &surface.eval_point(u, v));
        }
    }
}
//}}}
//}}}
//{{{ fun: assert_surfaces_agree
/// Asserts that two surfaces are the same map, by comparing points and low-order derivatives on a
/// grid covering the parameter domain.
fn assert_surfaces_agree(
    expected: &BSurface,
    actual: &BSurface,
) {
    for i in 0..=10 {
        for j in 0..=10 {
            let u = i as f64 / 10.0;
            let v = j as f64 / 10.0;

            assert_vectors_close(&actual.eval_point(u, v), &expected.eval_point(u, v));

            for (m, n) in [(1, 0), (0, 1), (1, 1), (2, 0), (0, 2)] {
                assert_vectors_close(
                    &actual.eval_derivative(u, v, m, n),
                    &expected.eval_derivative(u, v, m, n),
                );
            }
        }
    }
}
//}}}
//{{{ collection: knot insertion
//{{{ test: insert_single_knot_u
#[test]
fn insert_single_knot_u() {
    let surface = reference_surface();
    let mut refined = reference_surface();
    refined.insert_knot_u(0.4, 1);

    assert_eq!(refined.knots_u().len(), surface.knots_u().len() + 1);
    assert_eq!(refined.num_points_u(), surface.num_points_u() + 1);
    // The v-direction is untouched.
    assert_eq!(refined.order_v(), surface.order_v());
    assert_eq!(refined.knots_v(), surface.knots_v());
    assert_eq!(refined.num_points_v(), surface.num_points_v());

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_two_knots_u
#[test]
fn insert_two_knots_u() {
    let surface = reference_surface();
    let mut refined = reference_surface();
    refined.insert_knot_u(0.4, 2);

    assert_eq!(refined.knots_u().len(), surface.knots_u().len() + 2);
    assert_eq!(refined.num_points_u(), surface.num_points_u() + 2);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_existing_knot_u
#[test]
fn insert_existing_knot_u() {
    let surface = reference_surface();
    // 0.5 is already a knot of multiplicity one, taking it up to the order of the basis.
    let refined = surface.with_inserted_knot_u(0.5, 2);

    assert_eq!(refined.knots_u().len(), surface.knots_u().len() + 2);
    assert_eq!(refined.num_points_u(), surface.num_points_u() + 2);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_single_knot_v
#[test]
fn insert_single_knot_v() {
    let surface = reference_surface();
    let mut refined = reference_surface();
    refined.insert_knot_v(0.45, 1);

    assert_eq!(refined.knots_v().len(), surface.knots_v().len() + 1);
    assert_eq!(refined.num_points_v(), surface.num_points_v() + 1);
    // The u-direction is untouched.
    assert_eq!(refined.order_u(), surface.order_u());
    assert_eq!(refined.knots_u(), surface.knots_u());
    assert_eq!(refined.num_points_u(), surface.num_points_u());

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_two_knots_v
#[test]
fn insert_two_knots_v() {
    let surface = reference_surface();
    let mut refined = reference_surface();
    refined.insert_knot_v(0.45, 2);

    assert_eq!(refined.knots_v().len(), surface.knots_v().len() + 2);
    assert_eq!(refined.num_points_v(), surface.num_points_v() + 2);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_existing_knot_v
#[test]
fn insert_existing_knot_v() {
    let surface = reference_surface();
    // 0.3 is already a knot of multiplicity one, taking it up to the order of the basis.
    let refined = surface.with_inserted_knot_v(0.3, 1);

    assert_eq!(refined.knots_v().len(), surface.knots_v().len() + 1);
    assert_eq!(refined.num_points_v(), surface.num_points_v() + 1);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_knots_in_both_directions
#[test]
fn insert_knots_in_both_directions() {
    let surface = reference_surface();
    let mut refined = reference_surface();
    refined.insert_knot_u(0.4, 2);
    refined.insert_knot_v(0.45, 1);
    refined.insert_knot_u(0.1, 1);

    assert_eq!(refined.knots_u().len(), surface.knots_u().len() + 3);
    assert_eq!(refined.num_points_u(), surface.num_points_u() + 3);
    assert_eq!(refined.knots_v().len(), surface.knots_v().len() + 1);
    assert_eq!(refined.num_points_v(), surface.num_points_v() + 1);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//{{{ test: insert_knot_u_rejects_knot_outside_range
#[test]
#[should_panic(expected = "assertion failed: is_member(knots_before, u)")]
fn insert_knot_u_rejects_knot_outside_range() {
    reference_surface().insert_knot_u(1.5, 1);
}
//}}}
//{{{ test: insert_knot_v_rejects_multiplicity_greater_than_order
#[test]
#[should_panic(
    expected = "assertion failed: num_insertions >= 1 && num_insertions <= (ord - knot_mult)"
)]
fn insert_knot_v_rejects_multiplicity_greater_than_order() {
    reference_surface().insert_knot_v(0.3, 2);
}
//}}}
//{{{ test: with_inserted_knot_u_leaves_original_untouched
#[test]
fn with_inserted_knot_u_leaves_original_untouched() {
    let surface = reference_surface();
    let refined = surface.with_inserted_knot_u(0.4, 1);

    assert_eq!(surface.knots_u(), reference_surface().knots_u());
    assert_eq!(surface.num_points_u(), reference_surface().num_points_u());
    assert_eq!(refined.num_points_u(), surface.num_points_u() + 1);
}
//}}}
//}}}
//{{{ collection: knot refinement
//{{{ test: refine_knots_u_matches_repeated_insertion
#[test]
fn refine_knots_u_matches_repeated_insertion() {
    let new_knots = [0.1, 0.25, 0.4, 0.4, 0.6, 0.95];

    let surface = reference_surface();
    let refined = surface.with_refined_knots_u(&new_knots);

    // Inserting the same knots one at a time must give the same surface.
    let mut inserted = reference_surface();
    for u in new_knots {
        inserted.insert_knot_u(u, 1);
    }

    assert_eq!(
        refined.knots_u().len(),
        surface.knots_u().len() + new_knots.len()
    );
    assert_eq!(
        refined.num_points_u(),
        surface.num_points_u() + new_knots.len()
    );
    assert_eq!(refined.knots_u(), inserted.knots_u());

    assert_surfaces_agree(&surface, &refined);
    assert_surfaces_agree(&inserted, &refined);
}
//}}}
//{{{ test: refine_knots_v_matches_repeated_insertion
#[test]
fn refine_knots_v_matches_repeated_insertion() {
    let new_knots = [0.15, 0.3, 0.45, 0.8];

    let surface = reference_surface();
    let refined = surface.with_refined_knots_v(&new_knots);

    let mut inserted = reference_surface();
    for v in new_knots {
        inserted.insert_knot_v(v, 1);
    }

    assert_eq!(
        refined.knots_v().len(),
        surface.knots_v().len() + new_knots.len()
    );
    assert_eq!(
        refined.num_points_v(),
        surface.num_points_v() + new_knots.len()
    );
    assert_eq!(refined.knots_v(), inserted.knots_v());
    // The u-direction is untouched.
    assert_eq!(refined.knots_u(), surface.knots_u());
    assert_eq!(refined.num_points_u(), surface.num_points_u());

    assert_surfaces_agree(&surface, &refined);
    assert_surfaces_agree(&inserted, &refined);
}
//}}}
//{{{ test: refine_knots_in_both_directions
#[test]
fn refine_knots_in_both_directions() {
    let surface = reference_surface();
    let refined = surface
        .with_refined_knots_u(&[0.1, 0.4, 0.95])
        .with_refined_knots_v(&[0.15, 0.8]);

    assert_eq!(refined.num_points_u(), surface.num_points_u() + 3);
    assert_eq!(refined.num_points_v(), surface.num_points_v() + 2);

    assert_surfaces_agree(&surface, &refined);
}
//}}}
//}}}
//{{{ collection: degree elevation
//{{{ test: elevate_degree_u_by_one
#[test]
fn elevate_degree_u_by_one() {
    let surface = reference_surface();
    let num_unique_knots = knot_multiplicities(surface.knots_u()).len();

    let mut elevated = reference_surface();
    elevated.elevate_degree_u(1);

    assert_eq!(elevated.order_u(), surface.order_u() + 1);
    assert_eq!(
        elevated.knots_u().len(),
        surface.knots_u().len() + num_unique_knots
    );
    assert_eq!(
        elevated.num_points_u(),
        surface.num_points_u() + num_unique_knots - 1
    );
    // The v-direction is untouched.
    assert_eq!(elevated.order_v(), surface.order_v());
    assert_eq!(elevated.knots_v(), surface.knots_v());
    assert_eq!(elevated.num_points_v(), surface.num_points_v());

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_u_by_two
#[test]
fn elevate_degree_u_by_two() {
    let surface = reference_surface();
    let num_unique_knots = knot_multiplicities(surface.knots_u()).len();

    let elevated = surface.with_elevated_degree_u(2);

    assert_eq!(elevated.order_u(), surface.order_u() + 2);
    assert_eq!(
        elevated.knots_u().len(),
        surface.knots_u().len() + 2 * num_unique_knots
    );
    assert_eq!(
        elevated.num_points_u(),
        surface.num_points_u() + 2 * (num_unique_knots - 1)
    );

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_v_by_one
#[test]
fn elevate_degree_v_by_one() {
    let surface = reference_surface();
    let num_unique_knots = knot_multiplicities(surface.knots_v()).len();

    let mut elevated = reference_surface();
    elevated.elevate_degree_v(1);

    assert_eq!(elevated.order_v(), surface.order_v() + 1);
    assert_eq!(
        elevated.knots_v().len(),
        surface.knots_v().len() + num_unique_knots
    );
    assert_eq!(
        elevated.num_points_v(),
        surface.num_points_v() + num_unique_knots - 1
    );
    // The u-direction is untouched.
    assert_eq!(elevated.order_u(), surface.order_u());
    assert_eq!(elevated.knots_u(), surface.knots_u());
    assert_eq!(elevated.num_points_u(), surface.num_points_u());

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_v_by_two
#[test]
fn elevate_degree_v_by_two() {
    let surface = reference_surface();
    let num_unique_knots = knot_multiplicities(surface.knots_v()).len();

    let elevated = surface.with_elevated_degree_v(2);

    assert_eq!(elevated.order_v(), surface.order_v() + 2);
    assert_eq!(
        elevated.knots_v().len(),
        surface.knots_v().len() + 2 * num_unique_knots
    );
    assert_eq!(
        elevated.num_points_v(),
        surface.num_points_v() + 2 * (num_unique_knots - 1)
    );

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_by_zero_is_a_no_op
#[test]
fn elevate_degree_by_zero_is_a_no_op() {
    let surface = reference_surface();
    let elevated = surface.with_elevated_degree_u(0).with_elevated_degree_v(0);

    assert_eq!(elevated.order_u(), surface.order_u());
    assert_eq!(elevated.order_v(), surface.order_v());
    assert_eq!(elevated.knots_u(), surface.knots_u());
    assert_eq!(elevated.knots_v(), surface.knots_v());

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_in_both_directions
#[test]
fn elevate_degree_in_both_directions() {
    let surface = reference_surface();
    let elevated = surface.with_elevated_degree_u(1).with_elevated_degree_v(1);

    assert_eq!(elevated.order_u(), surface.order_u() + 1);
    assert_eq!(elevated.order_v(), surface.order_v() + 1);

    assert_surfaces_agree(&surface, &elevated);
}
//}}}
//{{{ test: elevate_degree_after_refinement
#[test]
fn elevate_degree_after_refinement() {
    let surface = reference_surface();
    let modified = surface
        .with_refined_knots_u(&[0.1, 0.4])
        .with_elevated_degree_u(1)
        .with_inserted_knot_v(0.45, 1)
        .with_elevated_degree_v(1);

    assert_eq!(modified.order_u(), surface.order_u() + 1);
    assert_eq!(modified.order_v(), surface.order_v() + 1);

    assert_surfaces_agree(&surface, &modified);
}
//}}}
//}}}
