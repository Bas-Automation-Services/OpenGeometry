//! One boolean for analytic solids, by the classic pipeline rather than by operand family.
//!
//! The stages, and what each takes from the kernel:
//!
//! 1. Intersect every face of one solid with every face of the other whose boxes meet
//!    (`intersect_breps`): closed-form curves where the surfaces have them (lines, circles,
//!    ellipses), traced intersection curves where not, each clipped to where it lies on both faces.
//! 2. Each such piece cuts the face it runs across, and its ends split the edges they lie on
//!    (`pieces`). A cut that runs over a cylinder's seam is split there, and so is the seam.
//! 3. Each face is cut along its cuts in its own flat (u, v) (`arrangement`): loops traced so the
//!    face stays on their left, at every node the first use clockwise from the way back; loops
//!    sorted into outer boundaries and the holes inside them.
//! 4. Each region is located against the other solid at a point well inside it (`Classifier`),
//!    and the regions the operation keeps (those of the second solid reversed in a subtraction)
//!    are sewn into one B-rep by their edges' geometry (`sew`): vertices where the curves through
//!    them meet, shells by adjacency, voids into the solid round them.
//!
//! Faces that lie on one surface (coincident) are refused for now.

mod arrangement;
mod curve;
mod flat;
mod pieces;
mod sew;

use super::assembly::analytic_face_mappings;
use super::types::{BooleanOp, BooleanReport, BooleanResult};
use crate::brep::{
    Accuracy, BrepEnvelope, FaceProvenance, FaceRole, FaceSource, GeometryError, GeometryQuality,
};
use crate::intersection::intersect_breps;
use crate::query::{Classifier, PointClassification};
use arrangement::{reversed_loop, ring};
use curve::CurveRef;
use flat::{inside_point, Flat};
use pieces::{Cut, Operand};
use sew::{sew, Kept};

fn gap(message: &str) -> GeometryError {
    GeometryError::CoverageGap {
        families: ["general boolean".into(), message.into()],
    }
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
    let mut contacts = Vec::new();
    for pair in &graph.pairs {
        if pair.graph.coincident {
            return Err(gap("faces on one surface"));
        }
        contacts.extend(pair.graph.contacts.iter().copied());
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
    let mut crossings = std::collections::HashMap::new();
    for operand in &mut operands {
        for face in 0..operand.brep.topology.faces.len() as u32 {
            operand.seam_crossings(face, &mut crossings)?;
        }
    }
    for operand in &mut operands {
        for face in 0..operand.brep.topology.faces.len() as u32 {
            operand.trace_cuts(face, &crossings)?;
        }
    }
    let classifiers = [Classifier::new(a), Classifier::new(b)];
    let mut kept = Vec::new();
    let mut keys: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for side in 0..2 {
        let operand = &operands[side];
        for face in 0..operand.brep.topology.faces.len() as u32 {
            for region in operand.regions(face)? {
                let map = &operand.maps[face as usize];
                let rings: Vec<Vec<Flat>> = region.loops.iter().map(|l| ring(l)).collect();
                let inside =
                    inside_point(&rings).ok_or_else(|| gap("no point found inside a face"))?;
                let place = classifiers[1 - side].classify(map.point(inside)?)?;
                let flip = match (operation, side, place) {
                    (_, _, PointClassification::Boundary) => {
                        return Err(gap("faces on one surface"))
                    }
                    (_, _, PointClassification::Unknown) => {
                        return Err(GeometryError::UnresolvedIntersection(
                            "general boolean: a face cannot be located against the other solid"
                                .into(),
                        ))
                    }
                    (BooleanOp::Union, _, PointClassification::Outside)
                    | (BooleanOp::Intersection, _, PointClassification::Inside)
                    | (BooleanOp::Subtraction, 0, PointClassification::Outside) => false,
                    (BooleanOp::Subtraction, 1, PointClassification::Inside) => true,
                    _ => continue,
                };
                let source = operand.brep.topology.faces[face as usize].clone();
                let role = if flip {
                    FaceRole::Cut
                } else if region.whole {
                    FaceRole::Preserved
                } else {
                    FaceRole::Split
                };
                kept.push(Kept {
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
                    key: {
                        let base = format!("{}:{}", operand.brep.id, source.key);
                        let count = keys.entry(base.clone()).or_insert(0);
                        *count += 1;
                        if *count == 1 {
                            base
                        } else {
                            format!("{base}-{count}")
                        }
                    },
                    provenance: FaceProvenance {
                        sources: vec![FaceSource {
                            entity: operand.brep.id.clone(),
                            body: operand.brep.id.clone(),
                            key: source.key,
                            face,
                        }],
                        role,
                        reversed: flip,
                    },
                    side,
                    inside,
                });
            }
        }
    }
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
            contacts,
            coincident: false,
            face_mappings,
        },
    })
}
