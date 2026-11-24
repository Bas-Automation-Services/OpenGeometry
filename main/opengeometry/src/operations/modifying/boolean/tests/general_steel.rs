use super::super::general::general_boolean;
use super::super::types::BooleanOp;
use crate::brep::{Accuracy, BrepEnvelope, Frame3};
use crate::exchange::export_step;
use crate::math::Point3;
use crate::primitives::{self, ProfileEdge};
use crate::tessellation::measure;
use std::f64::consts::{FRAC_PI_2, PI};

pub(super) const ACCURACY: Accuracy = Accuracy::STANDARD;

const H: f64 = 0.15;
const B: f64 = 0.075;
const TW: f64 = 0.00355;
const TF: f64 = 0.0107;
const R: f64 = 0.015;

fn line(from: [f64; 2], to: [f64; 2]) -> ProfileEdge {
    ProfileEdge::Line { from, to }
}
fn arc(center: [f64; 2], start_angle: f64, sweep_angle: f64) -> ProfileEdge {
    ProfileEdge::Arc {
        center,
        radius: R,
        start_angle,
        sweep_angle,
    }
}

pub(super) fn ipe_outline() -> Vec<ProfileEdge> {
    let (y0, y1) = (-H + TF, H - TF);
    vec![
        line([-B, -H], [B, -H]),
        line([B, -H], [B, y0]),
        line([B, y0], [TW + R, y0]),
        arc([TW + R, y0 + R], -FRAC_PI_2, -FRAC_PI_2),
        line([TW, y0 + R], [TW, y1 - R]),
        arc([TW + R, y1 - R], PI, -FRAC_PI_2),
        line([TW + R, y1], [B, y1]),
        line([B, y1], [B, H]),
        line([B, H], [-B, H]),
        line([-B, H], [-B, y1]),
        line([-B, y1], [-TW - R, y1]),
        arc([-TW - R, y1 - R], FRAC_PI_2, -FRAC_PI_2),
        line([-TW, y1 - R], [-TW, y0 + R]),
        arc([-TW - R, y0 + R], 0.0, -FRAC_PI_2),
        line([-TW - R, y0], [-B, y0]),
        line([-B, y0], [-B, -H]),
    ]
}

pub(super) fn ipe_area() -> f64 {
    2.0 * (2.0 * B * TF) + (2.0 * H - 2.0 * TF) * 2.0 * TW + 4.0 * R * R * (1.0 - PI / 4.0)
}

pub(super) fn ipe(id: &str, length: f64) -> BrepEnvelope {
    primitives::arc_edged_extrusion(id.into(), Frame3::IDENTITY, ipe_outline(), length, ACCURACY)
        .unwrap()
}

pub(super) fn half_space(id: &str, origin: Point3, normal: Point3, reach: f64) -> BrepEnvelope {
    let helper = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let frame = Frame3::from_axis(origin, normal, helper).unwrap();
    let corner = frame.point([-reach, -reach, 0.0]);
    let frame = Frame3 {
        origin: corner,
        ..frame
    };
    primitives::cuboid(
        id.into(),
        frame,
        [2.0 * reach, 2.0 * reach, reach],
        ACCURACY,
    )
    .unwrap()
}

fn assert_volume(body: &BrepEnvelope, exact: f64) {
    body.validate().unwrap();
    let m = measure(body, 1e-4).unwrap();
    assert!(
        (m.volume - exact).abs() <= m.error_bound + 1e-12,
        "volume {} vs exact {} (bound {})",
        m.volume,
        exact,
        m.error_bound
    );
}

fn cut(member: &BrepEnvelope, cutter: &BrepEnvelope) -> BrepEnvelope {
    general_boolean(member, cutter, BooleanOp::Subtraction, "cut".into())
        .unwrap()
        .brep
}

#[test]
fn a_cuboid_less_a_corner_cuboid() {
    let a = primitives::cuboid("a".into(), Frame3::IDENTITY, [2.0, 2.0, 2.0], ACCURACY).unwrap();
    let frame = Frame3 {
        origin: [1.0, 1.0, 1.0],
        ..Frame3::IDENTITY
    };
    let b = primitives::cuboid("b".into(), frame, [2.0, 2.0, 2.0], ACCURACY).unwrap();
    assert_volume(&cut(&a, &b), 7.0);
}

#[test]
fn an_ipe_mitred_at_45_degrees() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let member = ipe("member", 1.0);
    let cutter = half_space("cutter", [0.0, 0.0, 0.8], [0.0, s, s], 2.0);
    let result = cut(&member, &cutter);
    assert_volume(&result, 0.8 * ipe_area());
    export_step(&result, "metre").unwrap();
}

#[test]
fn an_ipe_cut_at_a_skew_angle() {
    let n = [0.3, 0.4, (1.0f64 - 0.25).sqrt()];
    let member = ipe("member", 1.0);
    let cutter = half_space("cutter", [0.0, 0.0, 0.7], n, 2.0);
    assert_volume(&cut(&member, &cutter), 0.7 * ipe_area());
}

fn shs_loop(half: f64, radius: f64) -> Vec<ProfileEdge> {
    let k = half - radius;
    let corner = |c: [f64; 2], start: f64| ProfileEdge::Arc {
        center: c,
        radius,
        start_angle: start,
        sweep_angle: FRAC_PI_2,
    };
    vec![
        line([-k, -half], [k, -half]),
        corner([k, -k], -FRAC_PI_2),
        line([half, -k], [half, k]),
        corner([k, k], 0.0),
        line([k, half], [-k, half]),
        corner([-k, k], FRAC_PI_2),
        line([-half, k], [-half, -k]),
        corner([-k, -k], PI),
    ]
}

fn rounded_square_area(half: f64, radius: f64) -> f64 {
    (2.0 * half).powi(2) - (4.0 - PI) * radius * radius
}

#[test]
fn a_hollow_section_mitred() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let member = primitives::arc_edged_extrusion_with_holes(
        "shs".into(),
        Frame3::IDENTITY,
        shs_loop(0.05, 0.0075),
        vec![shs_loop(0.045, 0.0025)],
        1.0,
        ACCURACY,
    )
    .unwrap();
    let cutter = half_space("cutter", [0.0, 0.0, 0.8], [s, 0.0, s], 2.0);
    let area = rounded_square_area(0.05, 0.0075) - rounded_square_area(0.045, 0.0025);
    let result = cut(&member, &cutter);
    assert_volume(&result, 0.8 * area);
    export_step(&result, "metre").unwrap();
}

#[test]
fn a_round_bar_cut_obliquely_across_its_seam() {
    let member = primitives::cylinder("bar".into(), Frame3::IDENTITY, 0.05, 1.0, ACCURACY).unwrap();
    for n in [[0.0, 0.6, 0.8], [0.6, 0.0, 0.8], [-0.48, 0.36, 0.8]] {
        let cutter = half_space("cutter", [0.0, 0.0, 0.6], n, 2.0);
        assert_volume(&cut(&member, &cutter), PI * 0.05 * 0.05 * 0.6);
    }
}
