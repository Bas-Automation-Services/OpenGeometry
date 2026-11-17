use super::arrangement::{faces, reversed_loop, Use};
use super::curve::CurveRef;
use super::flat::{distance_to_polyline, halfedge_flat, winding, FaceMap, Flat};
use super::imprint::Place;
use crate::brep::{BrepEnvelope, EdgeGeometry, GeometryError, Orientation, Surface};
use crate::math::{norm, sub, Point3};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// A curve told apart from others by the store it lives in and its index there.
pub(super) type CurveKey = (usize, u32);

pub(super) fn curve_key(curve: CurveRef<'_>) -> CurveKey {
    (curve.store as *const _ as usize, curve.id)
}

/// A piece of a face-pair intersection curve that runs across a face.
#[derive(Clone, Copy)]
pub(super) struct Cut<'a> {
    pub(super) curve: CurveRef<'a>,
    pub(super) lo: f64,
    pub(super) hi: f64,
}

/// One solid of a boolean: its faces unrolled, where its edges split and the cuts on its faces.
pub(super) struct Operand<'a> {
    pub(super) brep: &'a BrepEnvelope,
    pub(super) maps: Vec<FaceMap>,
    rings: Vec<Vec<Vec<Flat>>>,
    splits: Vec<Vec<f64>>,
    cuts: Vec<Vec<Cut<'a>>>,
    traced: Vec<Vec<Use<'a>>>,
    snap: f64,
}

/// A region of a face after the cuts: loops of uses with the face on their left seen from
/// outside the solid, and whether the face came through whole.
pub(super) struct Region<'a> {
    pub(super) face: u32,
    pub(super) loops: Vec<Vec<Use<'a>>>,
    pub(super) whole: bool,
}

fn gap(message: &str) -> GeometryError {
    GeometryError::CoverageGap {
        families: ["general boolean".into(), message.into()],
    }
}

impl<'a> Operand<'a> {
    pub(super) fn new(brep: &'a BrepEnvelope, snap: f64) -> Result<Self, GeometryError> {
        let mut operand = Self {
            brep,
            maps: Vec::new(),
            rings: Vec::new(),
            splits: vec![Vec::new(); brep.topology.edges.len()],
            cuts: vec![Vec::new(); brep.topology.faces.len()],
            traced: vec![Vec::new(); brep.topology.faces.len()],
            snap,
        };
        for face in &brep.topology.faces {
            let mut map = FaceMap::new(brep, face)?;
            let mut rings = Vec::new();
            for halfedges in operand.face_loops(face.id) {
                let mut ring = Vec::new();
                for h in halfedges {
                    let range = brep.topology.edges
                        [brep.topology.halfedges[h as usize].edge as usize]
                        .geometry
                        .range();
                    let curve = operand.edge_curve(brep.topology.halfedges[h as usize].edge)?;
                    ring.extend(halfedge_flat(
                        brep,
                        &map,
                        h,
                        &curve.stations(range.lo, range.hi),
                    )?);
                }
                rings.push(ring);
            }
            map.fit_window(&rings);
            operand.maps.push(map);
            operand.rings.push(rings);
        }
        Ok(operand)
    }

    /// The halfedges of each loop of a face, the outer loop first.
    fn face_loops(&self, face: u32) -> Vec<Vec<u32>> {
        let f = &self.brep.topology.faces[face as usize];
        std::iter::once(f.trim.outer)
            .chain(f.trim.holes.iter().copied())
            .map(|l| {
                let start = self.brep.topology.loops[l as usize].start_halfedge;
                let mut result = vec![start];
                let mut h = self.brep.topology.halfedges[start as usize].next;
                while let Some(next) = h.filter(|&n| n != start) {
                    result.push(next);
                    h = self.brep.topology.halfedges[next as usize].next;
                }
                result
            })
            .collect()
    }

    pub(super) fn edge_curve(&self, edge: u32) -> Result<CurveRef<'a>, GeometryError> {
        match self.brep.topology.edges[edge as usize].geometry {
            EdgeGeometry::Curve { curve, .. } => Ok(CurveRef {
                store: &self.brep.geometry,
                id: curve,
            }),
            EdgeGeometry::Collapsed { .. } => Err(gap("collapsed edge")),
        }
    }

    /// Whether a point on a face's surface lies on the face's boundary.
    pub(super) fn on_boundary(&self, face: u32, p: Point3) -> Result<bool, GeometryError> {
        let flat = self.flat_of(face, p)?;
        Ok(self.rings[face as usize].iter().any(|ring| {
            let mut closed = ring.clone();
            closed.push(ring[0]);
            distance_to_polyline(flat, &closed) <= self.snap
        }))
    }

    fn flat_of(&self, face: u32, p: Point3) -> Result<Flat, GeometryError> {
        let map = &self.maps[face as usize];
        let uv = map.surface.project(p, Some(map.window))?;
        let mut flat = map.flat(uv);
        for axis in 0..2 {
            if let Some(period) = map.period(axis) {
                let start = map.flat(map.window)[axis] - self.snap;
                flat[axis] = start + (flat[axis] - start).rem_euclid(period);
            }
        }
        Ok(flat)
    }

    pub(super) fn add_cut(&mut self, face: u32, cut: Cut<'a>) {
        self.cuts[face as usize].push(cut);
    }

    /// Splits whichever edge of the face the point lies on, unless at one of its ends.
    pub(super) fn split_at(&mut self, face: u32, p: Point3) -> Result<(), GeometryError> {
        for halfedges in self.face_loops(face) {
            for h in halfedges {
                self.split_edge_at(self.brep.topology.halfedges[h as usize].edge, p)?;
            }
        }
        Ok(())
    }

    /// Splits the edge at the point when it lies on it, away from its ends.
    pub(super) fn split_edge_at(&mut self, edge: u32, p: Point3) -> Result<(), GeometryError> {
        let range = self.brep.topology.edges[edge as usize].geometry.range();
        let curve = self.edge_curve(edge)?;
        let t = curve.parameter(p, 0.5 * (range.lo + range.hi))?;
        let t = if curve.periodic() {
            range.lo + (t - range.lo).rem_euclid(TAU)
        } else {
            t
        };
        if t < range.lo || t > range.hi || norm(sub(curve.point(t)?, p)) > self.snap {
            return Ok(());
        }
        let ends = [curve.point(range.lo)?, curve.point(range.hi)?];
        if ends.iter().any(|e| norm(sub(*e, p)) <= self.snap) {
            return Ok(());
        }
        let mut fresh = true;
        for s in &self.splits[edge as usize] {
            fresh &= norm(sub(curve.point(*s)?, p)) > self.snap;
        }
        if fresh {
            self.splits[edge as usize].push(t);
        }
        Ok(())
    }

    pub(super) fn snap(&self) -> f64 {
        self.snap
    }

    /// The face's edges, each once, with their curves and ranges.
    pub(super) fn face_edges(
        &self,
        face: u32,
    ) -> Result<Vec<(u32, CurveRef<'a>, f64, f64)>, GeometryError> {
        let mut result: Vec<(u32, CurveRef<'a>, f64, f64)> = Vec::new();
        for halfedges in self.face_loops(face) {
            for h in halfedges {
                let edge = self.brep.topology.halfedges[h as usize].edge;
                if result.iter().any(|(e, ..)| *e == edge) {
                    continue;
                }
                let range = self.brep.topology.edges[edge as usize].geometry.range();
                result.push((edge, self.edge_curve(edge)?, range.lo, range.hi));
            }
        }
        Ok(result)
    }

    /// An edge of the face in its flat, in the edge's own direction, with the sample parameters.
    pub(super) fn edge_flat(
        &self,
        face: u32,
        edge: u32,
    ) -> Result<(Vec<f64>, Vec<Flat>), GeometryError> {
        let h = self
            .face_loops(face)
            .into_iter()
            .flatten()
            .find(|&h| self.brep.topology.halfedges[h as usize].edge == edge)
            .ok_or_else(|| gap("edge not on face"))?;
        let range = self.brep.topology.edges[edge as usize].geometry.range();
        let stations = self.edge_curve(edge)?.stations(range.lo, range.hi);
        let flat = halfedge_flat(self.brep, &self.maps[face as usize], h, &stations)?;
        Ok((stations, flat))
    }

    /// Where a point on the face's surface lies against the face.
    pub(super) fn locate(&self, face: u32, p: Point3) -> Result<Place, GeometryError> {
        let flat = self.flat_of(face, p)?;
        let rings = &self.rings[face as usize];
        if rings.iter().any(|ring| {
            let mut closed = ring.clone();
            closed.push(ring[0]);
            distance_to_polyline(flat, &closed) <= self.snap
        }) {
            return Ok(Place::Boundary);
        }
        Ok(if winding(flat, rings) != 0 {
            Place::Inside
        } else {
            Place::Outside
        })
    }

    /// The outward normal of the face of this solid a point lies on, when it lies on one.
    pub(super) fn normal_at(&self, p: Point3) -> Result<Option<Point3>, GeometryError> {
        for face in 0..self.maps.len() {
            let map = &self.maps[face];
            let Ok(uv) = map.surface.project(p, Some(map.window)) else {
                continue;
            };
            if norm(sub(map.surface.point_at(uv)?, p)) > self.snap {
                continue;
            }
            if self.locate(face as u32, p)? != Place::Outside {
                return Ok(Some(map.normal(map.flat(uv))?));
            }
        }
        Ok(None)
    }

    /// Where the face's cuts run across its seams: the seam edges split there, and the cut curves'
    /// parameters recorded, so the cut splits there on every face it runs across.
    pub(super) fn seam_crossings(
        &mut self,
        face: u32,
        crossings: &mut HashMap<CurveKey, Vec<f64>>,
    ) -> Result<(), GeometryError> {
        let mut seams = Vec::new();
        for cut in self.cuts[face as usize].clone() {
            let mut start = cut.lo;
            loop {
                let map = &self.maps[face as usize];
                let (stations, flat) = self.windowed_trace(face, cut.curve, start, cut.hi)?;
                let Some((i, axis, side)) = map.exit(&flat) else {
                    break;
                };
                let t = crossing(map, cut.curve, &stations, &flat, i, axis, side)?;
                crossings.entry(curve_key(cut.curve)).or_default().push(t);
                seams.push(cut.curve.point(t)?);
                start = t;
            }
        }
        for p in seams {
            self.split_at(face, p)?;
        }
        Ok(())
    }

    /// A traced piece moved, when it starts on a seam, to the side it runs into.
    fn windowed_trace(
        &self,
        face: u32,
        curve: CurveRef<'_>,
        from: f64,
        to: f64,
    ) -> Result<(Vec<f64>, Vec<Flat>), GeometryError> {
        let map = &self.maps[face as usize];
        let (stations, mut flat) = map.trace(curve, from, to)?;
        if let Some((0, axis, side)) = map.exit(&flat) {
            if (flat[0][axis] - side).abs() <= self.snap {
                let period = map.period(axis).unwrap_or(0.0);
                let shift = if flat[1][axis] < side {
                    period
                } else {
                    -period
                };
                flat.iter_mut().for_each(|p| p[axis] += shift);
            }
        }
        Ok((stations, flat))
    }

    /// Traces every cut of the face into flat pieces, split where any face's seam splits it.
    pub(super) fn trace_cuts(
        &mut self,
        face: u32,
        crossings: &HashMap<CurveKey, Vec<f64>>,
    ) -> Result<(), GeometryError> {
        let mut pieces = Vec::new();
        for cut in &self.cuts[face as usize] {
            let mut splits: Vec<f64> = crossings
                .get(&curve_key(cut.curve))
                .map(|ts| {
                    ts.iter()
                        .map(|&t| {
                            if cut.curve.periodic() {
                                cut.lo + (t - cut.lo).rem_euclid(TAU)
                            } else {
                                t
                            }
                        })
                        .filter(|&t| t > cut.lo + 1e-12 && t < cut.hi - 1e-12)
                        .collect()
                })
                .unwrap_or_default();
            splits.sort_by(f64::total_cmp);
            splits.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
            let closed = cut.curve.periodic() && (cut.hi - cut.lo - TAU).abs() <= 1e-9;
            let stations: Vec<f64> = if closed && !splits.is_empty() {
                let mut s = splits.clone();
                s.push(splits[0] + TAU);
                s
            } else {
                std::iter::once(cut.lo)
                    .chain(splits)
                    .chain(std::iter::once(cut.hi))
                    .collect()
            };
            for w in stations.windows(2) {
                let (_, flat) = self.windowed_trace(face, cut.curve, w[0], w[1])?;
                if self.maps[face as usize].exit(&flat).is_some() {
                    return Err(gap("a cut runs across a seam it was not split at"));
                }
                pieces.push(Use {
                    curve: cut.curve,
                    from: w[0],
                    to: w[1],
                    flat,
                });
            }
        }
        self.traced[face as usize] = pieces;
        Ok(())
    }

    /// The face cut into regions, or as it is when nothing cut or split it.
    pub(super) fn regions(&self, face: u32) -> Result<Vec<Region<'a>>, GeometryError> {
        let topology = &self.brep.topology;
        let mut uses = Vec::new();
        let mut loops = Vec::new();
        let mut touched = !self.traced[face as usize].is_empty();
        for halfedges in self.face_loops(face) {
            let mut ring = Vec::new();
            for h in halfedges {
                let halfedge = &topology.halfedges[h as usize];
                let edge = halfedge.edge;
                let range = topology.edges[edge as usize].geometry.range();
                let curve = self.edge_curve(edge)?;
                let mut params = vec![range.lo];
                let mut splits = self.splits[edge as usize].clone();
                splits.sort_by(f64::total_cmp);
                touched |= !splits.is_empty();
                params.extend(splits);
                params.push(range.hi);
                let mut pieces: Vec<(f64, f64)> = params.windows(2).map(|w| (w[0], w[1])).collect();
                if halfedge.geometry_use.sense == Orientation::Reverse {
                    pieces = pieces.into_iter().rev().map(|(a, b)| (b, a)).collect();
                }
                for (from, to) in pieces {
                    let flat = halfedge_flat(
                        self.brep,
                        &self.maps[face as usize],
                        h,
                        &curve.stations(from, to),
                    )?;
                    ring.push(Use {
                        curve,
                        from,
                        to,
                        flat,
                    });
                }
            }
            uses.extend(ring.iter().cloned());
            loops.push(ring);
        }
        if !touched {
            return Ok(vec![Region {
                face,
                loops,
                whole: true,
            }]);
        }
        for piece in &self.traced[face as usize] {
            uses.push(piece.clone());
            uses.push(piece.reversed());
        }
        let reverse = self.maps[face as usize].sense == Orientation::Reverse;
        if reverse {
            uses = reversed_loop(&uses);
        }
        Ok(faces(uses, self.snap)?
            .into_iter()
            .map(|loops| Region {
                face,
                loops: if reverse {
                    loops.iter().map(|l| reversed_loop(l)).collect()
                } else {
                    loops
                },
                whole: false,
            })
            .collect())
    }
}

/// The parameter where a traced piece crosses the window's edge between samples `i` and `i + 1`.
fn crossing(
    map: &FaceMap,
    curve: CurveRef<'_>,
    stations: &[f64],
    flat: &[Flat],
    i: usize,
    axis: usize,
    side: f64,
) -> Result<f64, GeometryError> {
    let (mut lo, mut hi) = (stations[i], stations[i + 1]);
    let below = flat[i][axis] < side;
    let hint = map.uv(flat[i]);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let uv = map.surface.project(curve.point(mid)?, Some(hint))?;
        if (map.flat(uv)[axis] < side) == below {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}
