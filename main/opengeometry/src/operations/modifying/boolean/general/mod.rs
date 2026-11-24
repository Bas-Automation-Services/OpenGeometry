mod arrangement;
mod curve;
mod flat;
mod imprint;
mod pieces;
mod sew;

use super::assembly::analytic_face_mappings;
use super::types::{BooleanOp, BooleanReport, BooleanResult};
use crate::brep::{
    Accuracy, BrepEnvelope, FaceProvenance, FaceRole, FaceSource, GeometryError, GeometryQuality,
    PatchBounds,
};
use crate::intersection::{intersect_breps, BodyIntersectionGraph};
use crate::math::dot;
use crate::query::{Classifier, PointClassification};
use arrangement::{reversed_loop, ring};
use curve::CurveRef;
use flat::{inside_point, signed_area, FaceMap, Flat};
use imprint::imprint;
use pieces::{Cut, Operand, Region};
use sew::{sew, Kept};
use std::collections::HashMap;

#[derive(Clone, Copy)]
enum Located {
    Out,
    In,
    OnSame,
    OnOpposite,
}

fn gap(message: &str) -> GeometryError {
    GeometryError::CoverageGap {
        families: ["general boolean".into(), message.into()],
    }
}

fn unresolved(message: &str) -> GeometryError {
    GeometryError::UnresolvedIntersection(format!("general boolean: {message}"))
}

pub(crate) fn general_boolean(
    a: &BrepEnvelope,
    b: &BrepEnvelope,
    operation: BooleanOp,
    id: String,
) -> Result<BooleanResult, GeometryError> {
    let accuracy = Accuracy::combined(a.accuracy, b.accuracy);
    let snap = 100.0 * accuracy.geometric;
    let graph = intersect_breps(a, b)?;
    let mut operands = [Operand::new(a, snap)?, Operand::new(b, snap)?];
    cut(&mut operands, &graph)?;
    let kept = keep(&operands, operation, snap)?;
    let mut brep = if kept.is_empty() {
        BrepEnvelope::new(id, accuracy)?
    } else {
        sew(kept, id, accuracy, snap)?
    };
    brep.revision = a.revision.max(b.revision).saturating_add(1);
    let face_mappings = analytic_face_mappings(&brep, [a, b]);
    Ok(BooleanResult {
        brep,
        report: BooleanReport {
            operation,
            quality: GeometryQuality::Analytic,
            contacts: graph
                .pairs
                .iter()
                .flat_map(|pair| pair.graph.contacts.iter().copied())
                .collect(),
            coincident: false,
            face_mappings,
        },
    })
}

fn cut<'a>(
    operands: &mut [Operand<'a>; 2],
    graph: &'a BodyIntersectionGraph,
) -> Result<(), GeometryError> {
    for pair in &graph.pairs {
        if pair.graph.coincident {
            imprint(operands, pair.faces)?;
        }
        for branch in &pair.graph.branches {
            let curve = CurveRef {
                store: &pair.graph.geometry,
                id: branch.curve,
            };
            let (lo, hi) = (branch.range.lo, branch.range.hi);
            let mid = curve.point(0.5 * (lo + hi))?;
            for side in 0..2 {
                let face = pair.faces[side];
                if !operands[side].on_boundary(face, mid)? {
                    operands[side].add_cut(face, Cut { curve, lo, hi });
                }
                for p in branch.endpoints {
                    operands[side].split_at(face, p)?;
                }
            }
        }
    }
    let mut crossings = HashMap::new();
    for operand in operands.iter_mut() {
        for face in 0..operand.brep.topology.faces.len() as u32 {
            operand.seam_crossings(face, &mut crossings)?;
        }
    }
    for operand in operands.iter_mut() {
        for face in 0..operand.brep.topology.faces.len() as u32 {
            operand.trace_cuts(face, &crossings)?;
        }
    }
    Ok(())
}

fn keep<'k, 'a>(
    operands: &'k [Operand<'a>; 2],
    operation: BooleanOp,
    snap: f64,
) -> Result<Vec<Kept<'k, 'a>>, GeometryError> {
    let classifiers = [
        Classifier::new(operands[0].brep),
        Classifier::new(operands[1].brep),
    ];
    let boxes = [
        operands[0].brep.bounds_unchecked()?,
        operands[1].brep.bounds_unchecked()?,
    ];
    let mut kept = Vec::new();
    let mut keys = HashMap::new();
    for (side, operand) in operands.iter().enumerate() {
        for face in 0..operand.brep.topology.faces.len() as u32 {
            for region in operand.regions(face)? {
                let map = &operand.maps[face as usize];
                let rings: Vec<Vec<Flat>> = region.loops.iter().map(|l| ring(l)).collect();
                if rings.iter().any(|r| sliver(r, snap)) {
                    return Err(unresolved(
                        "a face piece is a sliver below geometric resolution",
                    ));
                }
                let inside =
                    inside_point(&rings).ok_or_else(|| gap("no point found inside a face"))?;
                let place = if region.whole && apart(&rings, map, boxes[1 - side], snap)? {
                    Located::Out
                } else {
                    locate(&operands[1 - side], &classifiers[1 - side], map, inside)?
                };
                let Some(flip) = flipped(operation, side, place) else {
                    continue;
                };
                kept.push(kept_region(operand, side, region, flip, inside, &mut keys));
            }
        }
    }
    Ok(kept)
}

fn locate(
    other: &Operand<'_>,
    classifier: &Classifier<'_>,
    map: &FaceMap,
    inside: Flat,
) -> Result<Located, GeometryError> {
    let p = map.point(inside)?;
    Ok(match classifier.classify(p)? {
        PointClassification::Boundary => {
            let normal = other
                .normal_at(p)?
                .ok_or_else(|| gap("a face on the other solid's boundary"))?;
            if dot(normal, map.normal(inside)?) > 0.0 {
                Located::OnSame
            } else {
                Located::OnOpposite
            }
        }
        PointClassification::Inside => Located::In,
        PointClassification::Outside => Located::Out,
        PointClassification::Unknown => {
            return Err(unresolved(
                "a face cannot be located against the other solid",
            ))
        }
    })
}

fn flipped(operation: BooleanOp, side: usize, place: Located) -> Option<bool> {
    match (operation, side, place) {
        (BooleanOp::Union, 0, Located::Out | Located::OnSame)
        | (BooleanOp::Union, 1, Located::Out)
        | (BooleanOp::Intersection, 0, Located::In | Located::OnSame)
        | (BooleanOp::Intersection, 1, Located::In)
        | (BooleanOp::Subtraction, 0, Located::Out | Located::OnOpposite) => Some(false),
        (BooleanOp::Subtraction, 1, Located::In) => Some(true),
        _ => None,
    }
}

fn kept_region<'k, 'a>(
    operand: &'k Operand<'a>,
    side: usize,
    region: Region<'a>,
    flip: bool,
    inside: Flat,
    keys: &mut HashMap<String, usize>,
) -> Kept<'k, 'a> {
    let face = region.face;
    let map = &operand.maps[face as usize];
    let source = &operand.brep.topology.faces[face as usize];
    let role = if flip {
        FaceRole::Cut
    } else if region.whole {
        FaceRole::Preserved
    } else {
        FaceRole::Split
    };
    let base = format!("{}:{}", operand.brep.id, source.key);
    let count = keys.entry(base.clone()).or_insert(0);
    *count += 1;
    let key = if *count == 1 {
        base
    } else {
        format!("{base}-{count}")
    };
    Kept {
        map,
        loops: if flip {
            region.loops.iter().map(|l| reversed_loop(l)).collect()
        } else {
            region.loops
        },
        sense: if flip {
            map.sense.reversed()
        } else {
            map.sense
        },
        key,
        provenance: FaceProvenance {
            sources: vec![FaceSource {
                entity: operand.brep.id.clone(),
                body: operand.brep.id.clone(),
                key: source.key.clone(),
                face,
            }],
            role,
            reversed: flip,
        },
        side,
        inside,
    }
}

fn apart(
    rings: &[Vec<Flat>],
    map: &FaceMap,
    bounds: Option<PatchBounds>,
    snap: f64,
) -> Result<bool, GeometryError> {
    let Some(bounds) = bounds else {
        return Ok(true);
    };
    for p in rings.iter().flatten() {
        let p = map.point(*p)?;
        if (0..3).any(|i| p[i] < bounds.axes[i].lo - snap || p[i] > bounds.axes[i].hi + snap) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn sliver(ring: &[Flat], snap: f64) -> bool {
    let length: f64 = (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            (b[0] - a[0]).hypot(b[1] - a[1])
        })
        .sum();
    signed_area(ring).abs() <= snap * length
}
