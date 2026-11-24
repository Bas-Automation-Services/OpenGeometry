use super::arrangement::Use;
use super::curve::{meeting, CurveRef};
use super::flat::{FaceMap, Flat};
use crate::brep::{
    boundary, Accuracy, BrepEnvelope, Builder, CurveGeometry, FaceProvenance, GeometryError,
    IntersectionSide, Orientation, PcurveGeometry, Shell, SolidRegion, SurfaceGeometry, UV,
};
use crate::math::{add, dot, norm, scale, sub, Interval, Point3};
use crate::operations::modifying::boolean::assembly::{assign_shell_faces, enclosing_solid};
use crate::query::{classify_point_in_shell, PointClassification};
use std::collections::HashMap;
use std::f64::consts::TAU;

pub(super) struct Kept<'k, 'a> {
    pub(super) map: &'k FaceMap,
    pub(super) loops: Vec<Vec<Use<'a>>>,
    pub(super) sense: Orientation,
    pub(super) key: String,
    pub(super) provenance: FaceProvenance,
    pub(super) side: usize,
    pub(super) inside: Flat,
}

struct EdgeSlot {
    edge: u32,
    curve: u32,
    range: Interval,
    definition: Option<u32>,
}

type UseKey = (usize, usize, usize);
type Welded<'a> = Vec<(Point3, Vec<(CurveRef<'a>, f64)>)>;
type Group = (usize, usize, Point3, Vec<UseKey>);
type Slots = HashMap<UseKey, (usize, Orientation)>;

fn unresolved(message: impl Into<String>) -> GeometryError {
    GeometryError::UnresolvedIntersection(format!("general boolean: {}", message.into()))
}

fn same_curve(a: CurveRef<'_>, b: CurveRef<'_>) -> bool {
    std::ptr::eq(a.store, b.store) && a.id == b.id
}

fn each_use<'k, 'a>(kept: &'k [Kept<'_, 'a>]) -> impl Iterator<Item = (UseKey, &'k Use<'a>)> {
    kept.iter().enumerate().flat_map(|(k, region)| {
        region.loops.iter().enumerate().flat_map(move |(l, ring)| {
            ring.iter()
                .enumerate()
                .map(move |(u, piece)| ((k, l, u), piece))
        })
    })
}

fn use_at<'k, 'a>(kept: &'k [Kept<'_, 'a>], (k, l, u): UseKey) -> &'k Use<'a> {
    &kept[k].loops[l][u]
}

pub(super) fn sew(
    kept: Vec<Kept<'_, '_>>,
    id: String,
    accuracy: Accuracy,
    snap: f64,
) -> Result<BrepEnvelope, GeometryError> {
    let (points, ends) = weld(&kept, snap)?;
    let mut builder = Builder::new(id, accuracy)?;
    let mut vertices = Vec::with_capacity(points.len());
    for (p, curves) in &points {
        vertices.push(builder.vertex(meeting(curves, *p)?));
    }
    let groups = group_edges(&kept, &ends, snap)?;
    let (slots, edges) = make_edges(&mut builder, &kept, &ends, &vertices, &groups)?;
    let wanted = make_faces(&mut builder, &kept, &ends, &vertices, &slots, &edges, snap)?;
    lift_intersection_uses(&mut builder, &wanted)?;
    finish(builder, &kept, snap)
}

fn weld<'a>(
    kept: &[Kept<'_, 'a>],
    snap: f64,
) -> Result<(Welded<'a>, HashMap<UseKey, [usize; 2]>), GeometryError> {
    let mut points: Welded<'a> = Vec::new();
    let mut ends = HashMap::new();
    for (key, piece) in each_use(kept) {
        let mut pair = [0; 2];
        for (slot, t) in [piece.from, piece.to].into_iter().enumerate() {
            let p = piece.curve.point(t)?;
            let index = match points.iter().position(|(q, _)| norm(sub(*q, p)) <= snap) {
                Some(i) => i,
                None => {
                    points.push((p, Vec::new()));
                    points.len() - 1
                }
            };
            let curves = &mut points[index].1;
            if !curves.iter().any(|(c, _)| same_curve(*c, piece.curve)) {
                curves.push((piece.curve, t));
            }
            pair[slot] = index;
        }
        ends.insert(key, pair);
    }
    Ok((points, ends))
}

fn group_edges(
    kept: &[Kept<'_, '_>],
    ends: &HashMap<UseKey, [usize; 2]>,
    snap: f64,
) -> Result<Vec<Group>, GeometryError> {
    let mut groups: Vec<Group> = Vec::new();
    for (key, piece) in each_use(kept) {
        let [a, b] = ends[&key];
        let mid = piece.curve.point(0.5 * (piece.from + piece.to))?;
        let found = groups.iter().position(|(x, y, m, _)| {
            ((*x, *y) == (a, b) || (*x, *y) == (b, a)) && norm(sub(*m, mid)) <= 16.0 * snap
        });
        match found {
            Some(g) => groups[g].3.push(key),
            None => groups.push((a, b, mid, vec![key])),
        }
    }
    Ok(groups)
}

fn make_edges(
    builder: &mut Builder,
    kept: &[Kept<'_, '_>],
    ends: &HashMap<UseKey, [usize; 2]>,
    vertices: &[u32],
    groups: &[Group],
) -> Result<(Slots, Vec<EdgeSlot>), GeometryError> {
    let positions: Vec<Point3> = vertices
        .iter()
        .map(|&v| builder.brep().topology.vertices[v as usize].position)
        .collect();
    let mut slots = HashMap::new();
    let mut edges = Vec::new();
    let mut definitions = HashMap::new();
    for (a, b, _, members) in groups {
        if members.len() != 2 {
            return Err(unresolved(format!(
                "an edge of the result has {} faces",
                members.len()
            )));
        }
        let (piece, other) = (use_at(kept, members[0]), use_at(kept, members[1]));
        let closed = a == b;
        let t0 = piece.curve.parameter(positions[*a], piece.from)?;
        let mut t1 = piece.curve.parameter(positions[*b], piece.to)?;
        if closed {
            t1 = t0 + if piece.to > piece.from { TAU } else { -TAU };
        }
        let (range, sense, start) = if t0 < t1 {
            (Interval::new(t0, t1)?, Orientation::Forward, *a)
        } else {
            (Interval::new(t1, t0)?, Orientation::Reverse, *b)
        };
        let (geometry, definition) = edge_geometry(builder, piece.curve, &mut definitions);
        let curve = builder.brep().geometry.curves.len() as u32;
        let edge = builder.edge(geometry, range, members[0].0 == members[1].0);
        let other_sense = if closed {
            closed_sense(piece, other, sense)?
        } else if ends[&members[1]][0] == start {
            Orientation::Forward
        } else {
            Orientation::Reverse
        };
        slots.insert(members[0], (edges.len(), sense));
        slots.insert(members[1], (edges.len(), other_sense));
        edges.push(EdgeSlot {
            edge,
            curve,
            range,
            definition,
        });
    }
    Ok((slots, edges))
}

fn edge_geometry(
    builder: &mut Builder,
    curve: CurveRef<'_>,
    definitions: &mut HashMap<(usize, u32), u32>,
) -> (CurveGeometry, Option<u32>) {
    let geometry = curve.geometry().clone();
    let CurveGeometry::Intersection { definition } = geometry else {
        return (geometry, None);
    };
    let key = (curve.store as *const _ as usize, definition);
    let index = *definitions.entry(key).or_insert_with(|| {
        let intersections = &mut builder.brep.geometry.intersections;
        intersections.push(curve.store.intersections[definition as usize].clone());
        intersections.len() as u32 - 1
    });
    (
        CurveGeometry::Intersection { definition: index },
        Some(index),
    )
}

fn closed_sense(
    piece: &Use<'_>,
    other: &Use<'_>,
    sense: Orientation,
) -> Result<Orientation, GeometryError> {
    let along = dot(
        piece.curve.tangent(piece.from)?,
        other.curve.tangent(other.from)?,
    ) * (piece.to - piece.from).signum()
        * (other.to - other.from).signum();
    Ok(if along > 0.0 { sense } else { sense.reversed() })
}

fn make_faces<'k>(
    builder: &mut Builder,
    kept: &'k [Kept<'k, '_>],
    ends: &HashMap<UseKey, [usize; 2]>,
    vertices: &[u32],
    slots: &Slots,
    edges: &[EdgeSlot],
    snap: f64,
) -> Result<Vec<(f64, UV, &'k FaceMap)>, GeometryError> {
    let mut wanted = Vec::new();
    for (k, region) in kept.iter().enumerate() {
        let surface_id = builder.brep().geometry.surfaces.len() as u32;
        let mut rings = Vec::new();
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for (l, ring) in region.loops.iter().enumerate() {
            let mut uses = Vec::new();
            for (u, piece) in ring.iter().enumerate() {
                for flat in &piece.flat {
                    let uv = region.map.uv(*flat);
                    for axis in 0..2 {
                        lo[axis] = lo[axis].min(uv[axis]);
                        hi[axis] = hi[axis].max(uv[axis]);
                    }
                }
                let (slot, sense) = slots[&(k, l, u)];
                let e = &edges[slot];
                let [a, b] = ends[&(k, l, u)];
                let pcurve = pcurve(region, piece, e, sense, builder, surface_id)?;
                if let Some(d) = e.definition {
                    builder.brep.geometry.intersections[d as usize].surfaces[region.side] =
                        surface_id;
                }
                wanted.push((piece.from, region.map.uv(piece.flat[0]), region.map));
                uses.push(boundary(e.edge, vertices[a], vertices[b], sense, pcurve));
            }
            rings.push(uses);
        }
        let pad = 4.0 * snap;
        let bounds = [[lo[0] - pad, hi[0] + pad], [lo[1] - pad, hi[1] + pad]];
        let outer = rings.remove(0);
        let face = builder.face_with_holes(
            &region.key,
            region.map.surface.clone(),
            bounds,
            outer,
            rings,
        )?;
        builder.set_face_sense(face, region.sense)?;
        builder.brep.topology.faces[face as usize].provenance = region.provenance.clone();
    }
    Ok(wanted)
}

fn pcurve(
    region: &Kept<'_, '_>,
    piece: &Use<'_>,
    e: &EdgeSlot,
    sense: Orientation,
    builder: &Builder,
    surface: u32,
) -> Result<PcurveGeometry, GeometryError> {
    let geometry = &builder.brep().geometry.curves[e.curve as usize];
    if let Some(definition) = e.definition {
        let side = if region.side == 0 {
            IntersectionSide::A
        } else {
            IntersectionSide::B
        };
        return Ok(PcurveGeometry::IntersectionSide { definition, side });
    }
    if let SurfaceGeometry::Plane { frame } = region.map.surface {
        let local = |p: Point3| [dot(p, frame.x), dot(p, frame.y)];
        match geometry {
            CurveGeometry::Line { origin, direction } => {
                return Ok(PcurveGeometry::Line2 {
                    origin: local(sub(*origin, frame.origin)),
                    direction: local(*direction),
                })
            }
            CurveGeometry::Circle { frame: c, radius } => {
                return Ok(PcurveGeometry::Conic2 {
                    origin: local(sub(c.origin, frame.origin)),
                    axis_a: local(scale(c.x, *radius)),
                    axis_b: local(scale(c.y, *radius)),
                })
            }
            CurveGeometry::Ellipse {
                frame: c,
                major_radius,
                minor_radius,
            } => {
                return Ok(PcurveGeometry::Conic2 {
                    origin: local(sub(c.origin, frame.origin)),
                    axis_a: local(scale(c.x, *major_radius)),
                    axis_b: local(scale(c.y, *minor_radius)),
                })
            }
            _ => {}
        }
    }
    let (first, last) = (piece.flat[0], piece.flat[piece.flat.len() - 1]);
    let (at_lo, at_hi) = if sense == Orientation::Forward {
        (first, last)
    } else {
        (last, first)
    };
    let (uv_lo, uv_hi) = (region.map.uv(at_lo), region.map.uv(at_hi));
    let width = e.range.width();
    Ok(PcurveGeometry::ProjectedCurve {
        curve: e.curve,
        surface,
        chart: 0,
        uv_hint: uv_lo,
        uv_rate: [(uv_hi[0] - uv_lo[0]) / width, (uv_hi[1] - uv_lo[1]) / width],
        parameter_origin: e.range.lo,
    })
}

fn lift_intersection_uses(
    builder: &mut Builder,
    wanted: &[(f64, UV, &FaceMap)],
) -> Result<(), GeometryError> {
    let brep = &mut builder.brep;
    for (h, &(t, uv, map)) in wanted.iter().enumerate() {
        let Some(pcurve) = brep.topology.halfedges[h].geometry_use.pcurve else {
            continue;
        };
        if !matches!(
            brep.geometry.pcurves[pcurve as usize],
            PcurveGeometry::IntersectionSide { .. }
        ) {
            continue;
        }
        let have = brep.geometry.pcurve_at(pcurve, t)?;
        let mut lift = [0; 2];
        for axis in 0..2 {
            if let Some(period) = map.uv_period(axis) {
                lift[axis] = ((uv[axis] - have[axis]) / period).round() as i32;
            }
        }
        brep.topology.halfedges[h].geometry_use.periodic_lift = lift;
    }
    Ok(())
}

fn finish(
    mut builder: Builder,
    kept: &[Kept<'_, '_>],
    snap: f64,
) -> Result<BrepEnvelope, GeometryError> {
    let count = builder.brep.topology.faces.len();
    builder
        .brep
        .orient_loops()
        .map_err(|e| unresolved(format!("loops: {e}")))?;
    let mut assigned = vec![false; count];
    let mut outers: Vec<(u32, f64)> = Vec::new();
    let mut cavities = Vec::new();
    for start in 0..count {
        if assigned[start] {
            continue;
        }
        let shell = builder.brep.topology.shells.len() as u32;
        let faces = assign_shell_faces(
            &mut builder.brep,
            &mut assigned,
            start,
            shell,
            || unresolved("the result has an open boundary"),
            || unresolved("the result has an open boundary"),
        )?;
        let face = faces[0] as usize;
        let map = kept[face].map;
        let p = map.point(kept[face].inside)?;
        let n = map.normal(kept[face].inside)?;
        let n = scale(n, kept[face].sense.multiplier() * map.sense.multiplier());
        let probe = add(p, scale(n, 4.0 * snap / norm(n)));
        builder.brep.topology.shells.push(Shell {
            id: shell,
            faces: faces.clone(),
            is_closed: true,
        });
        let single = faces.len() == count;
        match if single {
            PointClassification::Outside
        } else {
            classify_point_in_shell(&builder.brep, &faces, probe)?
        } {
            PointClassification::Inside => cavities.push(shell),
            PointClassification::Outside => outers.push((shell, shell_box(&builder.brep, &faces))),
            _ => return Err(unresolved("cannot tell a shell's inside from its outside")),
        }
    }
    for &(outer_shell, _) in &outers {
        builder.brep.solids.push(SolidRegion {
            outer_shell,
            cavity_shells: Vec::new(),
        });
    }
    for cavity in cavities {
        let index = enclosing_solid(
            &builder.brep,
            &outers,
            cavity,
            |brep, faces, point| {
                Ok(classify_point_in_shell(brep, faces, point)? == PointClassification::Inside)
            },
            || unresolved("a void lies in no solid"),
        )?;
        builder.brep.solids[index].cavity_shells.push(cavity);
    }
    builder
        .brep
        .orient_loops()
        .and_then(|()| builder.brep.validate())
        .map_err(|error| match error {
            GeometryError::InvalidTopology(m) | GeometryError::InvalidGeometry(m) => {
                unresolved(format!("the sewn result is not valid: {m}"))
            }
            other => other,
        })?;
    Ok(builder.brep)
}

fn shell_box(brep: &BrepEnvelope, faces: &[u32]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for h in &brep.topology.halfedges {
        if h.face.is_some_and(|f| faces.contains(&f)) {
            let p = brep.topology.vertices[h.from as usize].position;
            for i in 0..3 {
                lo[i] = lo[i].min(p[i]);
                hi[i] = hi[i].max(p[i]);
            }
        }
    }
    (0..3).map(|i| hi[i] - lo[i]).product()
}
