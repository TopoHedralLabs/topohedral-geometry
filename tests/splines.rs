use approx::assert_abs_diff_eq;
use topohedral_geometry::{
    knot_multiplicities, knot_multiplicity, spline_eval, spline_eval_diff, spline_eval_diff_all,
};

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
        let mut actual = vec![f64::NAN; p + 1];
        spline_eval(knots(p), u, p, &mut actual);
        assert_slice_eq(&actual, expected);
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
            let mut actual = vec![f64::NAN; row_size];
            spline_eval_diff(knots(p), U1, p, derivative, &mut actual);
            let start = derivative * row_size;
            assert_slice_eq(&actual, &expected[start..start + row_size]);
        }

        let mut actual_all = vec![f64::NAN; row_size * row_size];
        spline_eval_diff_all(knots(p), U1, p, p, &mut actual_all);
        assert_slice_eq(&actual_all, expected);
    }
}

#[test]
fn knot_multiplicities_match_reference_knot_vector() {
    assert_eq!(
        knot_multiplicities(KNOTS_P4),
        vec![(0.0, 5), (1.0, 1), (2.0, 1), (4.0, 2), (5.0, 1), (6.0, 5)],
    );
}

#[test]
fn knot_multiplicity_matches_reference_knot_vector() {
    let cases = [(0.0, 5), (1.0, 1), (2.0, 1), (4.0, 2), (5.0, 1), (6.0, 5)];

    for (knot, expected) in cases {
        assert_eq!(knot_multiplicity(KNOTS_P4, knot), expected);
    }

    assert_eq!(knot_multiplicity(KNOTS_P4, 3.0), 0);
    assert_eq!(knot_multiplicity(&[], 3.0), 0);
}

#[test]
fn knot_multiplicity_uses_tolerant_knot_comparison() {
    let one_ulp_above_four = f64::from_bits(4.0_f64.to_bits() + 1);

    assert_eq!(knot_multiplicity(KNOTS_P4, one_ulp_above_four), 2);
}
