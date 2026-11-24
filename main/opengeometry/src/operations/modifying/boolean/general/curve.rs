use crate::brep::{Curve, CurveGeometry, GeometryError, GeometryStore};
use crate::math::{add, dot, norm, scale, sub, Point3};
use std::f64::consts::TAU;

#[derive(Clone, Copy)]
pub(super) struct CurveRef<'a> {
    pub(super) store: &'a GeometryStore,
    pub(super) id: u32,
}

impl CurveRef<'_> {
    pub(super) fn geometry(&self) -> &CurveGeometry {
        &self.store.curves[self.id as usize]
    }
    pub(super) fn point(&self, t: f64) -> Result<Point3, GeometryError> {
        self.store.curve(self.id)?.point_at(t)
    }
    pub(super) fn tangent(&self, t: f64) -> Result<Point3, GeometryError> {
        self.store.curve(self.id)?.tangent_at(t)
    }
    pub(super) fn periodic(&self) -> bool {
        matches!(
            self.geometry(),
            CurveGeometry::Circle { .. } | CurveGeometry::Ellipse { .. }
        )
    }

    pub(super) fn stations(&self, lo: f64, hi: f64) -> Vec<f64> {
        let count = match self.geometry() {
            CurveGeometry::Line { .. } => 1,
            CurveGeometry::Circle { .. } | CurveGeometry::Ellipse { .. } => {
                ((hi - lo).abs() / TAU * 256.0).ceil().max(16.0) as usize
            }
            _ => 64,
        };
        (0..=count)
            .map(|i| lo + (hi - lo) * i as f64 / count as f64)
            .collect()
    }

    pub(super) fn parameter(&self, p: Point3, near: f64) -> Result<f64, GeometryError> {
        let mut t = match *self.geometry() {
            CurveGeometry::Line { origin, direction } => return Ok(dot(sub(p, origin), direction)),
            CurveGeometry::Circle { frame, .. } => {
                let local = frame.local(p);
                local[1].atan2(local[0])
            }
            CurveGeometry::Ellipse {
                frame,
                major_radius,
                minor_radius,
            } => {
                let local = frame.local(p);
                (local[1] / minor_radius).atan2(local[0] / major_radius)
            }
            _ => near,
        };
        if self.periodic() {
            t += ((near - t) / TAU).round() * TAU;
        }
        let ends = self.ends();
        for _ in 0..32 {
            let tangent = self.tangent(t)?;
            let speed = dot(tangent, tangent);
            if speed == 0.0 {
                break;
            }
            let next = t + dot(sub(p, self.point(t)?), tangent) / speed;
            let next = ends.map_or(next, |[lo, hi]| next.clamp(lo, hi));
            let step = next - t;
            t = next;
            if step.abs() <= 1e-15 * t.abs().max(1.0) {
                break;
            }
        }
        Ok(t)
    }

    fn ends(&self) -> Option<[f64; 2]> {
        let CurveGeometry::Intersection { definition } = *self.geometry() else {
            return None;
        };
        let anchors = &self.store.intersections[definition as usize].anchors;
        Some([anchors.first()?.parameter, anchors.last()?.parameter])
    }
}

pub(super) fn meeting(curves: &[(CurveRef<'_>, f64)], p: Point3) -> Result<Point3, GeometryError> {
    if curves.len() < 2 {
        return Ok(p);
    }
    let mut point = p;
    let mut params: Vec<f64> = curves.iter().map(|(_, t)| *t).collect();
    for _ in 0..64 {
        let mut sum = [0.0; 3];
        for ((curve, _), t) in curves.iter().zip(params.iter_mut()) {
            *t = curve.parameter(point, *t)?;
            sum = add(sum, curve.point(*t)?);
        }
        let next = scale(sum, 1.0 / curves.len() as f64);
        let moved = norm(sub(next, point));
        point = next;
        if moved <= 1e-16 {
            break;
        }
    }
    Ok(point)
}
