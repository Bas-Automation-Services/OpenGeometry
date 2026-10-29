//! Builds analytic solids from outside the crate, using only the public
//! `analytic::primitives::Builder` API.

use opengeometry::analytic::exchange::export_step;
use opengeometry::analytic::primitives::{annular_cylinder, plane_boundary, Builder};
use opengeometry::analytic::query::validate_orientation;
use opengeometry::analytic::topology::{Accuracy, Orientation};
use opengeometry::analytic::{CurveGeometry, Frame3, GeometryError, Point3, SurfaceGeometry};
use opengeometry::math::interval::Interval;

fn accuracy() -> Accuracy {
    Accuracy {
        geometric: 1e-7,
        intersection: 1e-8,
        tessellation: 1e-3,
        exchange: 1e-6,
    }
}

fn sub(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// An axis-aligned box `size` at the origin, built face by face.
fn box_solid(size: [f64; 3]) -> Result<Builder, GeometryError> {
    let mut b = Builder::new("box".into(), accuracy())?;
    let points: Vec<Point3> = (0..8)
        .map(|i| {
            std::array::from_fn(|axis| {
                if i & (1 << axis) == 0 {
                    0.0
                } else {
                    size[axis]
                }
            })
        })
        .collect();
    for p in &points {
        b.vertex(*p);
    }
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut edges = [[u32::MAX; 3]; 8];
    for i in 0..8 {
        for axis in 0..3 {
            if i & (1 << axis) == 0 {
                edges[i][axis] = b.edge(
                    CurveGeometry::Line {
                        origin: points[i],
                        direction: axes[axis],
                    },
                    Interval::new(0.0, size[axis])?,
                    false,
                );
            }
        }
    }
    for (key, idx, normal) in [
        ("bottom", [0, 2, 3, 1], [0.0, 0.0, -1.0]),
        ("top", [4, 5, 7, 6], [0.0, 0.0, 1.0]),
        ("left", [0, 4, 6, 2], [-1.0, 0.0, 0.0]),
        ("right", [1, 3, 7, 5], [1.0, 0.0, 0.0]),
        ("front", [0, 1, 5, 4], [0.0, -1.0, 0.0]),
        ("back", [2, 6, 7, 3], [0.0, 1.0, 0.0]),
    ] {
        let frame = Frame3::from_axis(points[idx[0]], normal, sub(points[idx[1]], points[idx[0]]))?;
        let mut bounds = [[f64::INFINITY, f64::NEG_INFINITY]; 2];
        let mut uses = Vec::new();
        for i in 0..4 {
            let (from, to) = (idx[i], idx[(i + 1) % 4]);
            let axis = (from ^ to).trailing_zeros() as usize;
            let sense = if from < to {
                Orientation::Forward
            } else {
                Orientation::Reverse
            };
            uses.push(plane_boundary(
                &b,
                frame,
                edges[from.min(to)][axis],
                from as u32,
                to as u32,
                sense,
            )?);
            let uv = frame.local(points[from]);
            for j in 0..2 {
                bounds[j][0] = bounds[j][0].min(uv[j] - 1e-7);
                bounds[j][1] = bounds[j][1].max(uv[j] + 1e-7);
            }
        }
        b.face(key, SurfaceGeometry::Plane { frame }, bounds, uses)?;
    }
    Ok(b)
}

#[test]
fn a_box_built_through_the_public_builder_validates_and_exports() {
    let brep = box_solid([2.0, 3.0, 4.0]).unwrap().finish().unwrap();

    assert_eq!(brep.topology.faces.len(), 6);
    assert_eq!(brep.topology.edges.len(), 12);
    assert!(brep
        .topology
        .edges
        .iter()
        .all(|e| e.twin_halfedge.is_some()));
    let keys: Vec<_> = brep.topology.faces.iter().map(|f| f.key.as_str()).collect();
    assert_eq!(keys, ["bottom", "top", "left", "right", "front", "back"]);
    assert!(brep
        .topology
        .faces
        .iter()
        .all(|f| f.provenance.sources[0].key == f.key));

    let (step, report) = export_step(&brep, "millimetre").unwrap();
    assert_eq!(report.faces, 6);
    assert_eq!(step.matches("=ADVANCED_FACE(").count(), 6);
}

#[test]
fn faces_get_sequential_ids_and_their_sense_can_be_set() {
    let mut b = Builder::new("one-face".into(), accuracy()).unwrap();
    let v: Vec<u32> = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        .into_iter()
        .map(|p| b.vertex(p))
        .collect();
    assert_eq!(v, [0, 1, 2]);

    assert!(matches!(
        b.set_face_sense(0, Orientation::Reverse),
        Err(GeometryError::MissingReference { .. })
    ));

    let mut b = box_solid([1.0, 1.0, 1.0]).unwrap();
    b.set_face_sense(5, Orientation::Reverse).unwrap();
    assert_eq!(b.brep().topology.faces[5].sense, Orientation::Reverse);
    // Face 5 now has its material normal pointing into the box. Turning its
    // loop to match leaves it running the same way as its neighbours' loops.
    let err = b.finish().unwrap_err();
    assert!(
        matches!(&err, GeometryError::InvalidTopology(m) if m.contains("twins") && m.contains("'back'")),
        "{err:?}"
    );
}

#[test]
fn orientation_check_accepts_reversed_bores_and_thin_plates() {
    // The bore of a tube is a cylinder face with sense Reverse: valid.
    let tube =
        annular_cylinder("tube".into(), Frame3::IDENTITY, 0.5, 1.0, 2.0, accuracy()).unwrap();
    validate_orientation(&tube).unwrap();

    // A 1 mm sheet, 1000 x 500 mm: probes must not step through the plate.
    let plate = box_solid([1000.0, 500.0, 1.0]).unwrap().finish().unwrap();
    validate_orientation(&plate).unwrap();
}

#[test]
fn an_edge_used_by_a_third_loop_is_rejected() {
    let mut b = box_solid([1.0, 1.0, 1.0]).unwrap();
    let frame = Frame3::from_axis([0.0; 3], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]).unwrap();
    // Edge 0 (vertex 0 -> 1) already has both of its uses.
    let reuse = plane_boundary(&b, frame, 0, 0, 1, Orientation::Forward).unwrap();
    let result = b.face(
        "extra",
        SurfaceGeometry::Plane { frame },
        [[-1.0, 1.0], [-1.0, 1.0]],
        vec![reuse],
    );
    assert!(matches!(result, Err(GeometryError::InvalidTopology(_))));
}
