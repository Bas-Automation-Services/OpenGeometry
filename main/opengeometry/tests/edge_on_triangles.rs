//! A hem's 180-degree bend with a hole through its corner leaves a sliver of
//! cylinder; its trim samples along one straight edge (collinear in uv up to
//! rounding) earcut into a triangle seen edge-on from the surface normal.
use opengeometry::brep::BrepEnvelope;
use opengeometry::tessellation::measure;

const SLIVER: &str = include_str!("fixtures/cases/hem-corner-sliver.brep.json");
/// The body's exact volume (from its sheet metal model).
const VOLUME: f64 = 8322.23524780223;

#[test]
fn a_triangle_seen_edge_on_turns_like_its_face() {
    let body: BrepEnvelope = serde_json::from_str(SLIVER).unwrap();
    for deflection in [0.02, 0.002] {
        let m = measure(&body, deflection).unwrap();
        assert!(
            (m.volume - VOLUME).abs() < 1e-3 * VOLUME,
            "{deflection}: {}",
            m.volume
        );
    }
}
