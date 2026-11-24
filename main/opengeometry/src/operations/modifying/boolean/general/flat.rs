use super::curve::CurveRef;
use crate::brep::{BrepEnvelope, Face, GeometryError, Orientation, Surface, SurfaceGeometry, UV};
use crate::math::Point3;

pub(super) type Flat = [f64; 2];

pub(super) struct FaceMap {
    pub(super) surface: SurfaceGeometry,
    pub(super) sense: Orientation,
    scale: [f64; 2],
    periods: [Option<f64>; 2],
    pub(super) window: [f64; 2],
}

impl FaceMap {
    pub(super) fn new(brep: &BrepEnvelope, face: &Face) -> Result<Self, GeometryError> {
        let surface = brep.geometry.surface(face.surface)?.clone();
        let scale = match surface {
            SurfaceGeometry::Cylinder { radius, .. } => [radius, 1.0],
            _ => [1.0, 1.0],
        };
        let periods = surface.charts()[face.trim.chart as usize].periods;
        let bounds = face.trim.uv_bounds;
        Ok(Self {
            surface,
            sense: face.sense,
            scale,
            periods,
            window: [bounds[0].lo, bounds[1].lo],
        })
    }

    pub(super) fn fit_window(&mut self, loops: &[Vec<Flat>]) {
        for axis in 0..2 {
            let low = loops
                .iter()
                .flatten()
                .map(|p| p[axis])
                .fold(f64::INFINITY, f64::min);
            if low.is_finite() {
                self.window[axis] = low / self.scale[axis];
            }
        }
    }

    pub(super) fn uv_period(&self, axis: usize) -> Option<f64> {
        self.periods[axis]
    }

    pub(super) fn period(&self, axis: usize) -> Option<f64> {
        self.periods[axis].map(|p| p * self.scale[axis])
    }

    pub(super) fn flat(&self, uv: UV) -> Flat {
        [uv[0] * self.scale[0], uv[1] * self.scale[1]]
    }
    pub(super) fn uv(&self, flat: Flat) -> UV {
        [flat[0] / self.scale[0], flat[1] / self.scale[1]]
    }
    pub(super) fn point(&self, flat: Flat) -> Result<Point3, GeometryError> {
        self.surface.point_at(self.uv(flat))
    }
    pub(super) fn normal(&self, flat: Flat) -> Result<Point3, GeometryError> {
        let n = self.surface.normal_at(self.uv(flat))?;
        Ok(n.map(|x| x * self.sense.multiplier()))
    }

    fn windowed(&self, mut uv: UV, slack: f64) -> UV {
        for axis in 0..2 {
            if let Some(period) = self.periods[axis] {
                let start = self.window[axis] - slack;
                uv[axis] = start + (uv[axis] - start).rem_euclid(period);
            }
        }
        uv
    }

    pub(super) fn trace(
        &self,
        curve: CurveRef<'_>,
        from: f64,
        to: f64,
    ) -> Result<(Vec<f64>, Vec<Flat>), GeometryError> {
        let stations = curve.stations(from, to);
        let mut uvs: Vec<UV> = Vec::with_capacity(stations.len());
        for &t in &stations {
            let hint = uvs.last().copied();
            uvs.push(self.surface.project(curve.point(t)?, hint)?);
        }
        let first = uvs[0];
        let moved = self.windowed(first, 1e-9);
        let shift = [moved[0] - first[0], moved[1] - first[1]];
        Ok((
            stations,
            uvs.iter()
                .map(|uv| self.flat([uv[0] + shift[0], uv[1] + shift[1]]))
                .collect(),
        ))
    }

    pub(super) fn exit(&self, flat: &[Flat]) -> Option<(usize, usize, f64)> {
        for axis in 0..2 {
            let Some(period) = self.periods[axis] else {
                continue;
            };
            let lo = self.window[axis] * self.scale[axis];
            let hi = (self.window[axis] + period) * self.scale[axis];
            let slack = 1e-9 * self.scale[axis];
            for i in 1..flat.len() {
                let x = flat[i][axis];
                if x < lo - slack {
                    return Some((i - 1, axis, lo));
                }
                if x > hi + slack {
                    return Some((i - 1, axis, hi));
                }
            }
        }
        None
    }
}

pub(super) fn halfedge_flat(
    brep: &BrepEnvelope,
    map: &FaceMap,
    halfedge: u32,
    stations: &[f64],
) -> Result<Vec<Flat>, GeometryError> {
    let h = &brep.topology.halfedges[halfedge as usize];
    let pcurve = h
        .geometry_use
        .pcurve
        .ok_or_else(|| GeometryError::InvalidTopology("face halfedge has no pcurve".into()))?;
    stations
        .iter()
        .map(|&t| {
            let mut uv = brep.geometry.pcurve_at(pcurve, t)?;
            for axis in 0..2 {
                if let Some(period) = map.periods[axis] {
                    uv[axis] += f64::from(h.geometry_use.periodic_lift[axis]) * period;
                }
            }
            Ok(map.flat(uv))
        })
        .collect()
}

pub(super) fn distance_to_polyline(p: Flat, line: &[Flat]) -> f64 {
    line.windows(2)
        .map(|w| segment_distance(p, w[0], w[1]))
        .fold(f64::INFINITY, f64::min)
}

fn segment_distance(p: Flat, a: Flat, b: Flat) -> f64 {
    let d = [b[0] - a[0], b[1] - a[1]];
    let len2 = d[0] * d[0] + d[1] * d[1];
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / len2).clamp(0.0, 1.0)
    };
    let q = [a[0] + t * d[0] - p[0], a[1] + t * d[1] - p[1]];
    q[0].hypot(q[1])
}

pub(super) fn winding(p: Flat, loops: &[Vec<Flat>]) -> i32 {
    let mut count = 0;
    for ring in loops {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if a[1] <= p[1] {
                if b[1] > p[1] && cross(a, b, p) > 0.0 {
                    count += 1;
                }
            } else if b[1] <= p[1] && cross(a, b, p) < 0.0 {
                count -= 1;
            }
        }
    }
    count
}

fn cross(a: Flat, b: Flat, p: Flat) -> f64 {
    (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])
}

pub(super) fn signed_area(ring: &[Flat]) -> f64 {
    (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        / 2.0
}

pub(super) fn inside_point(rings: &[Vec<Flat>]) -> Option<Flat> {
    let ys = rings.iter().flatten().map(|p| p[1]);
    let (low, high) = ys.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), y| {
        (lo.min(y), hi.max(y))
    });
    let mut best: Option<(Flat, f64)> = None;
    for f in [
        0.5137, 0.3291, 0.7243, 0.1517, 0.8669, 0.4211, 0.6083, 0.0711, 0.9377,
    ] {
        let y = low + (high - low) * f;
        let mut xs: Vec<f64> = Vec::new();
        for ring in rings {
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if (a[1] - y) * (b[1] - y) < 0.0 {
                    xs.push(a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]));
                }
            }
        }
        if xs.len() % 2 != 0 {
            continue;
        }
        xs.sort_by(f64::total_cmp);
        for pair in xs.chunks_exact(2) {
            let p = [(pair[0] + pair[1]) / 2.0, y];
            let clearance = rings
                .iter()
                .map(|ring| {
                    let mut closed = ring.clone();
                    closed.push(ring[0]);
                    distance_to_polyline(p, &closed)
                })
                .fold(f64::INFINITY, f64::min);
            if best.is_none_or(|(_, c)| clearance > c) {
                best = Some((p, clearance));
            }
        }
    }
    best.map(|(p, _)| p)
}
