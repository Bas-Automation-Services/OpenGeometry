// A structural I-section (an arc-edged extrusion along z, with root radii) coped at one end:
// a notch through its top flange, the notch tool an extrusion across the member along x.
use super::support::volume_at_deflection;
use crate::brep::{coarse_accuracy, BrepEnvelope, Frame3, GeometryError};
use crate::operations::modifying::boolean::dispatch::boolean_brep;
use crate::operations::modifying::boolean::types::BooleanOp;
use crate::primitives::{self, ProfileEdge};
use std::f64::consts::{FRAC_PI_2, PI};

// IPE300 in metres: depth, flange width, web and flange thickness, root radius; member length.
const H: f64 = 0.3;
const B: f64 = 0.15;
const TW: f64 = 0.0071;
const TF: f64 = 0.0107;
const R: f64 = 0.015;
const L: f64 = 1.5;
// Notch: back to the flange tip, its floor below the root radius, corner radius in the web.
const NOTCH_LENGTH: f64 = 0.075;
const NOTCH_FLOOR: f64 = 0.11;
const CORNER: f64 = 0.01;

fn line(from: [f64; 2], to: [f64; 2]) -> ProfileEdge {
    ProfileEdge::Line { from, to }
}

fn arc(center: [f64; 2], radius: f64, start_angle: f64, sweep_angle: f64) -> ProfileEdge {
    ProfileEdge::Arc {
        center,
        radius,
        start_angle,
        sweep_angle,
    }
}

/// The I-section, counter-clockwise; `r` = 0 gives sharp root corners.
fn i_section(r: f64) -> Vec<ProfileEdge> {
    let (b, h, w, f) = (B / 2.0, H / 2.0, TW / 2.0, TF);
    let mut edges = vec![line([-b, -h], [b, -h]), line([b, -h], [b, -h + f])];
    let mut fillet = |edges: &mut Vec<ProfileEdge>, center: [f64; 2], start: f64| {
        if r > 0.0 {
            edges.push(arc(center, r, start, -FRAC_PI_2));
        }
    };
    edges.push(line([b, -h + f], [w + r, -h + f]));
    fillet(&mut edges, [w + r, -h + f + r], -FRAC_PI_2);
    edges.push(line([w, -h + f + r], [w, h - f - r]));
    fillet(&mut edges, [w + r, h - f - r], PI);
    edges.extend([
        line([w + r, h - f], [b, h - f]),
        line([b, h - f], [b, h]),
        line([b, h], [-b, h]),
        line([-b, h], [-b, h - f]),
        line([-b, h - f], [-w - r, h - f]),
    ]);
    fillet(&mut edges, [-w - r, h - f - r], FRAC_PI_2);
    edges.push(line([-w, h - f - r], [-w, -h + f + r]));
    fillet(&mut edges, [-w - r, -h + f + r], 0.0);
    edges.extend([
        line([-w - r, -h + f], [-b, -h + f]),
        line([-b, -h + f], [-b, -h]),
    ]);
    edges
}

fn member(r: f64) -> BrepEnvelope {
    primitives::arc_edged_extrusion(
        "member".into(),
        Frame3::IDENTITY,
        i_section(r),
        L,
        coarse_accuracy(1e-5),
    )
    .unwrap()
}

/// The notch tool in the member's (z, y) plane, run along −x through the whole width: open
/// towards the z = L end and the top, its inner corner rounded by `corner` (0: square).
fn notch(corner: f64) -> BrepEnvelope {
    let (z, y, reach, half_width) = (L - NOTCH_LENGTH, NOTCH_FLOOR, 0.5, 0.2);
    let frame =
        Frame3::from_axis([half_width, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]).unwrap();
    let mut edges = vec![
        line([z + corner, y], [z + reach, y]),
        line([z + reach, y], [z + reach, y + reach]),
        line([z + reach, y + reach], [z, y + reach]),
        line([z, y + reach], [z, y + corner]),
    ];
    if corner > 0.0 {
        edges.push(arc([z + corner, y + corner], corner, PI, FRAC_PI_2));
    }
    primitives::arc_edged_extrusion(
        "notch".into(),
        frame,
        edges,
        2.0 * half_width,
        coarse_accuracy(1e-5),
    )
    .unwrap()
}

fn section_area(r: f64) -> f64 {
    2.0 * B * TF + TW * (H - 2.0 * TF) + 4.0 * r * r * (1.0 - PI / 4.0)
}

/// Section area above the notch floor: the flange, its two root radii and the web below them.
fn area_above_floor(r: f64) -> f64 {
    B * TF + 2.0 * r * r * (1.0 - PI / 4.0) + TW * (H / 2.0 - TF - NOTCH_FLOOR)
}

fn assert_volume(brep: &BrepEnvelope, exact: f64) {
    let volume = volume_at_deflection(brep, 1e-7);
    assert!((volume / exact - 1.0).abs() < 1e-6, "{volume} vs {exact}");
}

#[test]
fn a_square_notch_copes_an_i_section_with_root_radii_exactly() {
    let result = boolean_brep(
        &member(R),
        &notch(0.0),
        BooleanOp::Subtraction,
        "coped".into(),
    )
    .unwrap()
    .brep;
    assert_volume(
        &result,
        section_area(R) * L - area_above_floor(R) * NOTCH_LENGTH,
    );
}

#[test]
fn a_square_notch_copes_a_sharp_i_section_exactly() {
    let result = boolean_brep(
        &member(0.0),
        &notch(0.0),
        BooleanOp::Subtraction,
        "coped".into(),
    )
    .unwrap()
    .brep;
    assert_volume(
        &result,
        section_area(0.0) * L - area_above_floor(0.0) * NOTCH_LENGTH,
    );
}

/// The rounded corner sits wholly in the web: it keeps a web-thick sliver of the corner square
/// outside the quarter circle.
fn rounded_corner_kept() -> f64 {
    TW * CORNER * CORNER * (1.0 - PI / 4.0)
}

#[test]
#[ignore = "CoverageGap: no handler yet for an arc-edged extrusion minus an arc-edged extrusion across its axis"]
fn a_rounded_notch_copes_an_i_section_exactly() {
    let result = boolean_brep(
        &member(R),
        &notch(CORNER),
        BooleanOp::Subtraction,
        "coped".into(),
    )
    .unwrap()
    .brep;
    assert_volume(
        &result,
        section_area(R) * L - area_above_floor(R) * NOTCH_LENGTH + rounded_corner_kept(),
    );
}

#[test]
fn a_rounded_notch_is_a_coverage_gap_on_an_i_section_today() {
    for r in [0.0, R] {
        let result = boolean_brep(
            &member(r),
            &notch(CORNER),
            BooleanOp::Subtraction,
            "coped".into(),
        );
        assert!(
            matches!(result, Err(GeometryError::CoverageGap { .. })),
            "r = {r}"
        );
    }
}
