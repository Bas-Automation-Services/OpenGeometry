//! A bend sector with a slanted side, built through the public builder: the
//! slanted side is a wrapped plane bounded by two wrapped lines (helices).

use opengeometry::analytic::exchange::export_step_in;
use opengeometry::analytic::geometry::Wrap;
use opengeometry::analytic::primitives::{boundary, plane_boundary, uv_line, Builder};
use opengeometry::analytic::query::{classify_point, PointClassification};
use opengeometry::analytic::tessellation::tessellate;
use opengeometry::analytic::topology::{Accuracy, Orientation};
use opengeometry::analytic::{BrepEnvelope, CurveGeometry, Frame3, Point3, SurfaceGeometry};
use opengeometry::math::interval::Interval;
use opengeometry::units::LengthUnit;

const R: f64 = 2.0; // inner radius
const T: f64 = 2.0; // thickness
const RN: f64 = 2.8; // wrap radius (the neutral layer)
const L: f64 = 4.0; // flat length around the axis
const S: f64 = 1.0; // the slanted side rises this much over L
const W: f64 = 10.0; // flat width along the axis

fn accuracy() -> Accuracy {
    Accuracy {
        geometric: 1e-6,
        intersection: 1e-7,
        tessellation: 1e-3,
        exchange: 1e-5,
    }
}

fn wrap() -> Wrap {
    Wrap {
        frame: Frame3::from_axis([1.0, 2.0, 3.0], [0.2, 1.0, 0.1], [1.0, 0.0, 0.3]).unwrap(),
        radius: RN,
    }
}

fn wrapped(curve: CurveGeometry) -> CurveGeometry {
    CurveGeometry::Wrapped {
        wrap: wrap(),
        curve: Box::new(curve),
    }
}

/// The flat region (0,0) (L,S) (L,W) (0,W), z = R .. R + T, wrapped.
fn sector() -> BrepEnvelope {
    let w = wrap();
    let f = w.frame;
    let mut b = Builder::new("sector".into(), accuracy()).unwrap();
    let corners = [[0.0, 0.0], [L, S], [L, W], [0.0, W]];
    let z = [R, R + T];
    // v[corner][level]
    let v: Vec<[u32; 2]> = corners
        .iter()
        .map(|c| z.map(|z| b.vertex(w.point([c[0], c[1], z]))))
        .collect();
    let slant = (L * L + S * S).sqrt();
    let d = [L / slant, S / slant, 0.0];
    let radial = |x: f64| f.vector([(x / RN).cos(), (x / RN).sin(), 0.0]);
    let line = |b: &mut Builder, origin: Point3, direction: Point3, length: f64| {
        b.edge(
            CurveGeometry::Line { origin, direction },
            Interval::new(0.0, length).unwrap(),
            false,
        )
    };
    // Per level: AB (wrapped line, from A), BC and DA (axial lines from B and
    // A), CD (arc from D around to C).
    let mut ab = [0; 2];
    let mut bc = [0; 2];
    let mut cd = [0; 2];
    let mut da = [0; 2];
    for (level, &z) in z.iter().enumerate() {
        ab[level] = b.edge(
            wrapped(CurveGeometry::Line {
                origin: [0.0, 0.0, z],
                direction: d,
            }),
            Interval::new(0.0, slant).unwrap(),
            false,
        );
        bc[level] = line(&mut b, w.point([L, S, z]), f.z, W - S);
        da[level] = line(&mut b, w.point([0.0, 0.0, z]), f.z, W);
        cd[level] = b.edge(
            CurveGeometry::Circle {
                frame: Frame3 {
                    origin: f.point([0.0, 0.0, W]),
                    ..f
                },
                radius: z,
            },
            Interval::new(0.0, L / RN).unwrap(),
            false,
        );
    }
    // Radial edges at each corner, from the inner level out.
    let up: Vec<u32> = corners
        .iter()
        .map(|c| line(&mut b, w.point([c[0], c[1], R]), radial(c[0]), T))
        .collect();
    let fwd = Orientation::Forward;
    let rev = Orientation::Reverse;
    let pad = 1e-3;

    // Inner and outer cylinders: UV = (angle, axial).
    for (level, &z) in z.iter().enumerate() {
        let uses = vec![
            boundary(
                ab[level],
                v[0][level],
                v[1][level],
                fwd,
                uv_line([0.0, 0.0], [d[0] / RN, d[1]]),
            ),
            boundary(
                bc[level],
                v[1][level],
                v[2][level],
                fwd,
                uv_line([L / RN, S], [0.0, 1.0]),
            ),
            boundary(
                cd[level],
                v[2][level],
                v[3][level],
                rev,
                uv_line([0.0, W], [1.0, 0.0]),
            ),
            boundary(
                da[level],
                v[3][level],
                v[0][level],
                rev,
                uv_line([0.0, 0.0], [0.0, 1.0]),
            ),
        ];
        let face = b
            .face(
                if level == 0 { "inner" } else { "outer" },
                SurfaceGeometry::Cylinder {
                    frame: f,
                    radius: z,
                },
                [[-pad, L / RN + pad], [-pad, W + pad]],
                uses,
            )
            .unwrap();
        if level == 0 {
            b.set_face_sense(face, rev).unwrap();
        }
    }

    // Planar ends and the far side.
    let plane =
        |b: &mut Builder, key: &str, frame: Frame3, loop_: &[(u32, u32, u32, Orientation)]| {
            let uses = loop_
                .iter()
                .map(|&(e, from, to, sense)| plane_boundary(b, frame, e, from, to, sense).unwrap())
                .collect();
            let big = 4.0 * (R + T + W);
            b.face(
                key,
                SurfaceGeometry::Plane { frame },
                [[-big, big], [-big, big]],
                uses,
            )
            .unwrap();
        };
    plane(
        &mut b,
        "start",
        Frame3::from_axis(f.origin, f.vector([0.0, -1.0, 0.0]), f.x).unwrap(),
        &[
            (da[0], v[0][0], v[3][0], fwd),
            (up[3], v[3][0], v[3][1], fwd),
            (da[1], v[3][1], v[0][1], rev),
            (up[0], v[0][1], v[0][0], rev),
        ],
    );
    let end = L / RN;
    plane(
        &mut b,
        "end",
        Frame3::from_axis(f.origin, f.vector([-end.sin(), end.cos(), 0.0]), radial(L)).unwrap(),
        &[
            (bc[0], v[1][0], v[2][0], fwd),
            (up[2], v[2][0], v[2][1], fwd),
            (bc[1], v[2][1], v[1][1], rev),
            (up[1], v[1][1], v[1][0], rev),
        ],
    );
    plane(
        &mut b,
        "far",
        Frame3 {
            origin: f.point([0.0, 0.0, W]),
            ..f
        },
        &[
            (cd[0], v[2][0], v[3][0], rev),
            (up[3], v[3][0], v[3][1], fwd),
            (cd[1], v[3][1], v[2][1], fwd),
            (up[2], v[2][1], v[2][0], rev),
        ],
    );

    // The slanted side: a flat vertical plane through A and B, normal
    // pointing out of the region; UV = (along AB, z).
    let n = [S / slant, -L / slant, 0.0];
    let flat = Frame3 {
        origin: [0.0; 3],
        x: d,
        y: [0.0, 0.0, 1.0],
        z: n,
    };
    let uses = vec![
        boundary(ab[0], v[0][0], v[1][0], fwd, uv_line([0.0, R], [1.0, 0.0])),
        boundary(
            up[1],
            v[1][0],
            v[1][1],
            fwd,
            uv_line([slant, R], [0.0, 1.0]),
        ),
        boundary(
            ab[1],
            v[1][1],
            v[0][1],
            rev,
            uv_line([0.0, R + T], [1.0, 0.0]),
        ),
        boundary(up[0], v[0][1], v[0][0], rev, uv_line([0.0, R], [0.0, 1.0])),
    ];
    b.face(
        "slant",
        SurfaceGeometry::Wrapped {
            wrap: w,
            surface: Box::new(SurfaceGeometry::Plane { frame: flat }),
        },
        [[-pad, slant + pad], [R - pad, R + T + pad]],
        uses,
    )
    .unwrap();
    b.finish().unwrap()
}

/// Flat area times the radial integral of the Jacobian z / RN.
fn exact_volume() -> f64 {
    let area = L * W - L * S / 2.0;
    area * ((R + T).powi(2) - R * R) / (2.0 * RN)
}

#[test]
fn a_wrapped_side_validates_and_meshes_to_the_exact_volume() {
    let brep = sector();
    let mesh = tessellate(&brep, 1e-4, 10_000_000).unwrap();
    let p = |i: u32| -> Point3 { std::array::from_fn(|j| mesh.positions[i as usize * 3 + j]) };
    let volume: f64 = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum();
    let exact = exact_volume();
    assert!(
        ((volume - exact) / exact).abs() < 1e-5,
        "{volume} vs {exact}"
    );
}

#[test]
fn points_on_either_side_of_the_wrapped_wall_classify_correctly() {
    let brep = sector();
    let w = wrap();
    // Just inside and just outside the slanted wall, mid-thickness.
    let on = |x: f64| [x, x * S / L, R + T / 2.0];
    for x in [0.5, 2.0, 3.5] {
        let p = on(x);
        let inside = w.point([p[0], p[1] + 0.01, p[2]]);
        let outside = w.point([p[0], p[1] - 0.01, p[2]]);
        assert_eq!(
            classify_point(&brep, inside).unwrap(),
            PointClassification::Inside
        );
        assert_eq!(
            classify_point(&brep, outside).unwrap(),
            PointClassification::Outside
        );
    }
    let middle = w.point([L / 2.0, W / 2.0, R + T / 2.0]);
    assert_eq!(
        classify_point(&brep, middle).unwrap(),
        PointClassification::Inside
    );
}

#[test]
fn a_wrapped_side_exports_to_step_as_b_splines_within_budget() {
    let brep = sector();
    let (step, report) =
        export_step_in(&brep, LengthUnit::Millimetre, LengthUnit::Millimetre).unwrap();
    assert_eq!(report.faces, 6);
    assert_eq!(step.matches("B_SPLINE_SURFACE_WITH_KNOTS").count(), 1);
    assert_eq!(step.matches("B_SPLINE_CURVE_WITH_KNOTS").count(), 2);
    assert!(report.exchange_error_bound <= accuracy().exchange);
}
