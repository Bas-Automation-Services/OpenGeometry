use super::super::dispatch::boolean_brep_outcome_with_handlers;
use super::super::types::BooleanOp;
use crate::brep::{Accuracy, BrepEnvelope, Frame3};
use crate::primitives::{self, ProfileEdge};
use crate::tessellation::measure;

const ACCURACY: Accuracy = Accuracy {
    geometric: 6.3e-8,
    intersection: 1.575e-8,
    tessellation: 0.01,
    exchange: 1e-5,
};

fn polar(radius: f64, angle: f64) -> [f64; 2] {
    [radius * angle.cos(), radius * angle.sin()]
}

fn arc(radius: f64, start_angle: f64, sweep_angle: f64) -> ProfileEdge {
    ProfileEdge::Arc {
        center: [0.0, 0.0],
        radius,
        start_angle,
        sweep_angle,
    }
}

fn band() -> BrepEnvelope {
    let (outer, inner, start, sweep) = (3.0, 2.7, 0.8, -1.6);
    let profile = vec![
        arc(outer, start, sweep),
        ProfileEdge::Line {
            from: polar(outer, start + sweep),
            to: polar(inner, start + sweep),
        },
        arc(inner, start + sweep, -sweep),
        ProfileEdge::Line {
            from: polar(inner, start),
            to: polar(outer, start),
        },
    ];
    let ground = Frame3 {
        x: [1.0, 0.0, 0.0],
        y: [0.0, 0.0, -1.0],
        z: [0.0, 1.0, 0.0],
        ..Frame3::IDENTITY
    };
    primitives::arc_edged_extrusion("band".into(), ground, profile, 3.0, ACCURACY).unwrap()
}

fn oblique_opening() -> BrepEnvelope {
    let (sin, cos) = 0.5f64.sin_cos();
    let frame = Frame3 {
        origin: [2.2, 1.0, 0.9],
        x: [cos, 0.0, -sin],
        y: [0.0, 1.0, 0.0],
        z: [sin, 0.0, cos],
    };
    primitives::cuboid("oblique".into(), frame, [1.0, 0.5, 0.4], ACCURACY).unwrap()
}

fn round_opening() -> BrepEnvelope {
    let frame = Frame3 {
        origin: [2.65, 1.5, 0.0],
        x: [0.0, 0.0, -1.0],
        y: [0.0, 1.0, 0.0],
        z: [1.0, 0.0, 0.0],
    };
    primitives::cylinder("round".into(), frame, 0.5, 0.4, ACCURACY).unwrap()
}

fn boolean(a: &BrepEnvelope, b: &BrepEnvelope, operation: BooleanOp) -> BrepEnvelope {
    let (result, handlers) = boolean_brep_outcome_with_handlers(a, b, operation, "out".into());
    assert_eq!(handlers, ["general_boolean"], "{operation:?}");
    result.unwrap().brep
}

#[test]
fn a_curved_band_less_an_oblique_then_a_round_opening_meets_a_traced_curve_at_its_end() {
    let (host, _) = boolean_brep_outcome_with_handlers(
        &band(),
        &oblique_opening(),
        BooleanOp::Subtraction,
        "host".into(),
    );
    let host = host.unwrap().brep;
    let round = round_opening();
    let rest = boolean(&host, &round, BooleanOp::Subtraction);
    let common = boolean(&host, &round, BooleanOp::Intersection);
    rest.validate().unwrap();
    let [whole, rest, common] = [&host, &rest, &common].map(|body| measure(body, 1e-4).unwrap());
    let missing = whole.volume - rest.volume - common.volume;
    assert!(
        missing.abs() <= 0.01 * common.volume,
        "{} - {} - {}",
        whole.volume,
        rest.volume,
        common.volume
    );
}
