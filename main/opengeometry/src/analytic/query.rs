use super::{
    geometry::{add, cross, dot, norm, scale, sub, unit, PatchBounds, Surface},
    topology::{BrepEnvelope, Face, PcurveGeometry},
    GeometryError, Point3, SurfaceGeometry, UV,
};
use crate::math::{
    interval::Interval,
    roots::{isolate_candidates, quadratic, QuadraticRoots},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PointClassification {
    Inside,
    Outside,
    Boundary,
    Unknown,
}

fn lifted_uv(brep: &BrepEnvelope, halfedge: u32, parameter: f64) -> Result<UV, GeometryError> {
    let use_ = &brep.topology.halfedges[halfedge as usize];
    let face = &brep.topology.faces[use_
        .face
        .ok_or_else(|| GeometryError::InvalidTopology("trim query halfedge has no face".into()))?
        as usize];
    let pcurve = use_.geometry_use.pcurve.ok_or_else(|| {
        GeometryError::InvalidTopology("trim query halfedge has no pcurve".into())
    })?;
    let mut uv = brep.geometry.pcurve_at(pcurve, parameter)?;
    let periods = brep.geometry.surface(face.surface)?.charts()[face.trim.chart as usize].periods;
    for axis in 0..2 {
        uv[axis] += periods[axis].map_or(0.0, |period| {
            period * f64::from(use_.geometry_use.periodic_lift[axis])
        });
    }
    Ok(uv)
}

fn trim_samples(brep: &BrepEnvelope, face: &Face, loop_id: u32) -> Result<Vec<UV>, GeometryError> {
    let loop_ = &brep.topology.loops[loop_id as usize];
    let start = loop_.start_halfedge;
    let mut current = start;
    let mut points = Vec::new();
    for _ in 0..=brep.topology.halfedges.len() {
        let use_ = &brep.topology.halfedges[current as usize];
        if use_.face != Some(face.id) {
            return Err(GeometryError::InvalidTopology(
                "trim query crossed into another face".into(),
            ));
        }
        let edge = &brep.topology.edges[use_.edge as usize];
        let range = edge.geometry.range();
        let pcurve = use_.geometry_use.pcurve.ok_or_else(|| {
            GeometryError::InvalidTopology("trim query halfedge has no pcurve".into())
        })?;
        let subdivisions = match brep.geometry.pcurves[pcurve as usize] {
            PcurveGeometry::Line2 { direction, .. } => {
                let periods =
                    brep.geometry.surface(face.surface)?.charts()[face.trim.chart as usize].periods;
                if periods.iter().enumerate().any(|(axis, period)| {
                    period.is_some_and(|period| {
                        ((direction[axis] * range.width()).abs() - period).abs()
                            <= 4.0 * brep.accuracy.geometric
                    })
                }) {
                    64
                } else {
                    1
                }
            }
            PcurveGeometry::Conic2 { .. } => 64,
            PcurveGeometry::ProjectedCurve { .. } | PcurveGeometry::IntersectionSide { .. } => 128,
        };
        for index in 0..subdivisions {
            let fraction = if use_.geometry_use.sense == super::topology::Orientation::Forward {
                index as f64 / subdivisions as f64
            } else {
                1.0 - index as f64 / subdivisions as f64
            };
            points.push(lifted_uv(
                brep,
                current,
                range.lo + range.width() * fraction,
            )?);
        }
        current = use_
            .next
            .ok_or_else(|| GeometryError::InvalidTopology("open trim loop in query".into()))?;
        if current == start {
            break;
        }
    }
    if current != start || points.len() < 3 {
        return Err(GeometryError::InvalidTopology(
            "trim query requires a closed nondegenerate loop".into(),
        ));
    }
    Ok(points)
}

fn segment_distance_2d(point: UV, a: UV, b: UV) -> f64 {
    let direction = [b[0] - a[0], b[1] - a[1]];
    let length_squared = dot(
        [direction[0], direction[1], 0.0],
        [direction[0], direction[1], 0.0],
    );
    if length_squared == 0.0 {
        return (point[0] - a[0]).hypot(point[1] - a[1]);
    }
    let parameter = (((point[0] - a[0]) * direction[0] + (point[1] - a[1]) * direction[1])
        / length_squared)
        .clamp(0.0, 1.0);
    (point[0] - a[0] - parameter * direction[0]).hypot(point[1] - a[1] - parameter * direction[1])
}

fn in_loop(point: UV, polygon: &[UV], tolerance: f64) -> Option<bool> {
    let mut inside = false;
    for index in 0..polygon.len() {
        let a = polygon[index];
        let b = polygon[(index + 1) % polygon.len()];
        if segment_distance_2d(point, a, b) <= tolerance {
            return None;
        }
        if (a[1] > point[1]) != (b[1] > point[1]) {
            let x = a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if x > point[0] {
                inside = !inside;
            }
        }
    }
    Some(inside)
}

fn lift_to_face(surface: &SurfaceGeometry, face: &Face, mut uv: UV) -> UV {
    let periods = surface.charts()[face.trim.chart as usize].periods;
    for axis in 0..2 {
        if let Some(period) = periods[axis] {
            let center = face.trim.uv_bounds[axis].midpoint();
            uv[axis] += ((center - uv[axis]) / period).round() * period;
        }
    }
    uv
}

fn periodic_cylinder_band_limit(
    brep: &BrepEnvelope,
    face: &Face,
    loop_id: u32,
) -> Result<Option<f64>, GeometryError> {
    let loop_ = &brep.topology.loops[loop_id as usize];
    let use_ = &brep.topology.halfedges[loop_.start_halfedge as usize];
    if use_.next != Some(loop_.start_halfedge) || use_.from != use_.to {
        return Ok(None);
    }
    let edge = &brep.topology.edges[use_.edge as usize];
    let range = edge.geometry.range();
    let Some(pcurve) = use_.geometry_use.pcurve else {
        return Ok(None);
    };
    let PcurveGeometry::Line2 { direction, .. } = brep.geometry.pcurves[pcurve as usize] else {
        return Ok(None);
    };
    let period = brep.geometry.surface(face.surface)?.charts()[face.trim.chart as usize].periods[0];
    let Some(period) = period else {
        return Ok(None);
    };
    if direction[1].abs() * range.width() > brep.accuracy.geometric
        || ((direction[0] * range.width()).abs() - period).abs() > 4.0 * brep.accuracy.geometric
    {
        return Ok(None);
    }
    Ok(Some(lifted_uv(brep, loop_.start_halfedge, range.lo)?[1]))
}

pub(crate) fn face_contains_uv(
    brep: &BrepEnvelope,
    face: &Face,
    uv: UV,
) -> Result<Option<bool>, GeometryError> {
    let surface = brep.geometry.surface(face.surface)?;
    let uv = lift_to_face(surface, face, uv);
    let jet = surface.derivatives(uv)?;
    let metric = norm(jet.du).min(norm(jet.dv));
    let uv_tolerance = if metric > 0.0 {
        brep.accuracy.geometric / metric
    } else {
        return Ok(None);
    };
    if !face.trim.uv_bounds[0].contains(uv[0]) || !face.trim.uv_bounds[1].contains(uv[1]) {
        return Ok(Some(false));
    }
    if matches!(surface, SurfaceGeometry::Cylinder { .. }) && face.trim.holes.len() == 1 {
        if let (Some(outer), Some(inner)) = (
            periodic_cylinder_band_limit(brep, face, face.trim.outer)?,
            periodic_cylinder_band_limit(brep, face, face.trim.holes[0])?,
        ) {
            let lo = outer.min(inner);
            let hi = outer.max(inner);
            if (uv[1] - lo).abs() <= uv_tolerance || (uv[1] - hi).abs() <= uv_tolerance {
                return Ok(None);
            }
            return Ok(Some(uv[1] > lo && uv[1] < hi));
        }
    }
    let outer = trim_samples(brep, face, face.trim.outer)?;
    let Some(mut inside) = in_loop(uv, &outer, uv_tolerance) else {
        return Ok(None);
    };
    if !inside {
        return Ok(Some(false));
    }
    for hole in &face.trim.holes {
        let polygon = trim_samples(brep, face, *hole)?;
        match in_loop(uv, &polygon, uv_tolerance) {
            None => return Ok(None),
            Some(true) => inside = false,
            Some(false) => {}
        }
    }
    Ok(Some(inside))
}

fn polynomial_value(coefficients: &[f64], value: f64) -> f64 {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |sum, coefficient| sum.mul_add(value, *coefficient))
}

fn refine_polynomial_root(coefficients: &[f64], interval: Interval) -> f64 {
    let mut t = interval.midpoint();
    let derivative = coefficients
        .iter()
        .enumerate()
        .skip(1)
        .map(|(degree, coefficient)| *coefficient * degree as f64)
        .collect::<Vec<_>>();
    for _ in 0..12 {
        let value = polynomial_value(coefficients, t);
        let slope = polynomial_value(&derivative, t);
        if slope == 0.0 || !slope.is_finite() {
            break;
        }
        let candidate = t - value / slope;
        if !candidate.is_finite() || candidate < interval.lo || candidate > interval.hi {
            break;
        }
        t = candidate;
    }
    t
}

fn support_roots(
    surface: &SurfaceGeometry,
    origin: Point3,
    direction: Point3,
    domain: Interval,
    tolerance: f64,
) -> Result<Option<Vec<f64>>, GeometryError> {
    let o = surface.frame().local(origin);
    let d = [
        dot(direction, surface.frame().x),
        dot(direction, surface.frame().y),
        dot(direction, surface.frame().z),
    ];
    let quadratic_roots = |a, b, c| -> Result<Option<Vec<f64>>, GeometryError> {
        Ok(Some(match quadratic(a, b, c)? {
            QuadraticRoots::Empty => Vec::new(),
            QuadraticRoots::One(root) => vec![root],
            QuadraticRoots::Two(roots) => roots.to_vec(),
            QuadraticRoots::IdenticallyZero => return Ok(None),
        }))
    };
    let result = match surface {
        SurfaceGeometry::Plane { .. } => {
            if d[2] == 0.0 {
                if o[2].abs() <= tolerance {
                    None
                } else {
                    Some(Vec::new())
                }
            } else {
                Some(vec![-o[2] / d[2]])
            }
        }
        SurfaceGeometry::Sphere { radius, .. } => {
            quadratic_roots(dot(d, d), 2.0 * dot(o, d), dot(o, o) - radius * radius)?
        }
        SurfaceGeometry::Cylinder { radius, .. } => quadratic_roots(
            d[0].mul_add(d[0], d[1] * d[1]),
            2.0 * o[0].mul_add(d[0], o[1] * d[1]),
            o[0].mul_add(o[0], o[1] * o[1]) - radius * radius,
        )?,
        SurfaceGeometry::Cone { semi_angle, .. } => {
            let k2 = semi_angle.tan().powi(2);
            quadratic_roots(
                d[0].mul_add(d[0], d[1] * d[1]) - k2 * d[2] * d[2],
                2.0 * (o[0].mul_add(d[0], o[1] * d[1]) - k2 * o[2] * d[2]),
                o[0].mul_add(o[0], o[1] * o[1]) - k2 * o[2] * o[2],
            )?
            .map(|roots| {
                roots
                    .into_iter()
                    .filter(|root| o[2] + root * d[2] >= -tolerance)
                    .collect()
            })
        }
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            let radial = [
                o[0].mul_add(o[0], o[1] * o[1]),
                2.0 * o[0].mul_add(d[0], o[1] * d[1]),
                d[0].mul_add(d[0], d[1] * d[1]),
            ];
            let sum = [
                dot(o, o) + major_radius * major_radius - minor_radius * minor_radius,
                2.0 * dot(o, d),
                dot(d, d),
            ];
            let r4 = 4.0 * major_radius * major_radius;
            let coefficients = [
                sum[0] * sum[0] - r4 * radial[0],
                2.0 * sum[0] * sum[1] - r4 * radial[1],
                sum[1] * sum[1] + 2.0 * sum[0] * sum[2] - r4 * radial[2],
                2.0 * sum[1] * sum[2],
                sum[2] * sum[2],
            ];
            let candidates = isolate_candidates(&coefficients, domain, tolerance, 100_000)
                .map_err(|error| {
                    GeometryError::UnresolvedIntersection(format!(
                        "torus ray root isolation: {error:?}"
                    ))
                })?;
            Some(
                candidates
                    .intervals
                    .into_iter()
                    .map(|candidate| refine_polynomial_root(&coefficients, candidate))
                    .collect(),
            )
        }
        // No implicit form: see wrapped_ray_roots.
        SurfaceGeometry::Wrapped { .. } => None,
    };
    Ok(result.map(|roots| {
        roots
            .into_iter()
            .filter(|root| domain.contains(*root))
            .collect()
    }))
}

/// The ray parameters where o + t d lies in a box, within `domain`.
fn ray_in_box(bounds: PatchBounds, o: Point3, d: Point3, domain: Interval) -> Option<Interval> {
    let (mut lo, mut hi) = (domain.lo, domain.hi);
    for axis in 0..3 {
        let (a, b) = (bounds.axes[axis].lo, bounds.axes[axis].hi);
        if d[axis] == 0.0 {
            if o[axis] < a || o[axis] > b {
                return None;
            }
            continue;
        }
        let (mut t0, mut t1) = ((a - o[axis]) / d[axis], (b - o[axis]) / d[axis]);
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }
        lo = lo.max(t0);
        hi = hi.min(t1);
    }
    Interval::new(lo, hi).ok().filter(|_| lo <= hi)
}

/// Ray roots on a wrapped surface within a UV box: subdivide while the
/// surface's conservative box meets the ray, then Newton from each small box.
/// `None` when the search doesn't settle.
fn wrapped_ray_roots(
    surface: &SurfaceGeometry,
    uv: [Interval; 2],
    o: Point3,
    d: Point3,
    domain: Interval,
    tolerance: f64,
) -> Result<Option<Vec<f64>>, GeometryError> {
    let diagonal = |b: PatchBounds| norm(std::array::from_fn(|i| b.axes[i].width()));
    let leaf = diagonal(surface.enclose(uv)?) / 64.0;
    let mut roots: Vec<f64> = Vec::new();
    let mut stack = vec![(uv, 0usize)];
    let mut visits = 0usize;
    while let Some((b, depth)) = stack.pop() {
        visits += 1;
        if visits > 50_000 {
            return Ok(None);
        }
        let bounds = surface.enclose(b)?;
        let Some(span) = ray_in_box(bounds, o, d, domain) else {
            continue;
        };
        if diagonal(bounds) <= leaf || depth >= 24 {
            let start = [b[0].midpoint(), b[1].midpoint(), span.midpoint()];
            if let Some([u, v, t]) = ray_newton(surface, start, o, d, tolerance)? {
                if uv[0].contains(u)
                    && uv[1].contains(v)
                    && domain.contains(t)
                    && roots.iter().all(|r| (r - t).abs() > tolerance)
                {
                    roots.push(t);
                }
            }
            continue;
        }
        let jet = surface.derivatives([b[0].midpoint(), b[1].midpoint()])?;
        let axis = usize::from(norm(jet.dv) * b[1].width() > norm(jet.du) * b[0].width());
        let m = b[axis].midpoint();
        for half in [Interval::new(b[axis].lo, m)?, Interval::new(m, b[axis].hi)?] {
            let mut child = b;
            child[axis] = half;
            stack.push((child, depth + 1));
        }
    }
    Ok(Some(roots))
}

/// Solves S(u, v) = o + t d by Newton's method from `start` = (u, v, t).
fn ray_newton(
    surface: &SurfaceGeometry,
    start: [f64; 3],
    o: Point3,
    d: Point3,
    tolerance: f64,
) -> Result<Option<[f64; 3]>, GeometryError> {
    let mut x = start;
    for _ in 0..32 {
        let jet = match surface.derivatives([x[0], x[1]]) {
            Ok(jet) => jet,
            Err(GeometryError::SingularParameterization) => return Ok(None),
            Err(error) => return Err(error),
        };
        let r = sub(jet.point, add(o, scale(d, x[2])));
        if norm(r) <= tolerance {
            return Ok(Some(x));
        }
        // [du dv -d] dx = -r, by Cramer's rule.
        let m = [jet.du, jet.dv, scale(d, -1.0)];
        let det = dot(m[0], cross(m[1], m[2]));
        if det.abs() <= f64::EPSILON * norm(m[0]) * norm(m[1]) {
            return Ok(None);
        }
        let rhs = scale(r, -1.0);
        let step = [
            dot(rhs, cross(m[1], m[2])) / det,
            dot(m[0], cross(rhs, m[2])) / det,
            dot(m[0], cross(m[1], rhs)) / det,
        ];
        for i in 0..3 {
            x[i] += step[i];
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Ok(None);
        }
    }
    Ok(None)
}

fn ray_box_domain(
    brep: &BrepEnvelope,
    origin: Point3,
    direction: Point3,
) -> Result<Option<Interval>, GeometryError> {
    let Some(bounds) = brep.bounds_unchecked()? else {
        return Ok(None);
    };
    let mut lo: f64 = 0.0;
    let mut hi = f64::MAX;
    for axis in 0..3 {
        if direction[axis] == 0.0 {
            if !bounds.axes[axis].contains(origin[axis]) {
                return Ok(None);
            }
            continue;
        }
        let mut a = (bounds.axes[axis].lo - origin[axis]) / direction[axis];
        let mut b = (bounds.axes[axis].hi - origin[axis]) / direction[axis];
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        lo = lo.max(a);
        hi = hi.min(b);
    }
    if hi < lo.max(0.0) {
        Ok(None)
    } else {
        Ok(Some(Interval::new(lo.max(0.0), hi)?))
    }
}

fn classify_with_ray<'a>(
    brep: &BrepEnvelope,
    point: Point3,
    direction: Point3,
    faces: impl Iterator<Item = &'a Face>,
) -> Result<PointClassification, GeometryError> {
    let Some(domain) = ray_box_domain(brep, point, direction)? else {
        return Ok(PointClassification::Outside);
    };
    let mut hits: Vec<(f64, Point3)> = Vec::new();
    for face in faces {
        let surface = brep.geometry.surface(face.surface)?;
        let roots = if matches!(surface, SurfaceGeometry::Wrapped { .. }) {
            wrapped_ray_roots(
                surface,
                face.trim.uv_bounds,
                point,
                direction,
                domain,
                brep.accuracy.intersection,
            )?
        } else {
            support_roots(
                surface,
                point,
                direction,
                domain,
                brep.accuracy.intersection,
            )?
        };
        let Some(roots) = roots else {
            return Ok(PointClassification::Unknown);
        };
        for parameter in roots {
            let position = [
                point[0] + parameter * direction[0],
                point[1] + parameter * direction[1],
                point[2] + parameter * direction[2],
            ];
            let uv = match surface.project(position, None) {
                Ok(uv) => uv,
                Err(super::GeometryError::SingularParameterization) => {
                    return Ok(PointClassification::Unknown)
                }
                Err(error) => return Err(error),
            };
            if norm(sub(surface.point_at(uv)?, position)) > brep.accuracy.intersection {
                continue;
            }
            match face_contains_uv(brep, face, uv)? {
                Some(true) => {}
                Some(false) => continue,
                None => return Ok(PointClassification::Unknown),
            }
            let normal = scale(surface.normal_at(uv)?, face.sense.multiplier());
            if dot(normal, direction).abs() <= 64.0 * f64::EPSILON {
                return Ok(PointClassification::Unknown);
            }
            if parameter <= brep.accuracy.geometric {
                return Ok(PointClassification::Boundary);
            }
            if hits
                .iter()
                .all(|(existing, _)| (existing - parameter).abs() > brep.accuracy.geometric)
            {
                hits.push((parameter, position));
            }
        }
    }
    Ok(if hits.len() % 2 == 1 {
        PointClassification::Inside
    } else {
        PointClassification::Outside
    })
}

pub fn classify_point(
    brep: &BrepEnvelope,
    point: Point3,
) -> Result<PointClassification, GeometryError> {
    brep.validate()?;
    classify_point_validated(brep, point)
}

pub(crate) fn classify_point_validated(
    brep: &BrepEnvelope,
    point: Point3,
) -> Result<PointClassification, GeometryError> {
    if point.into_iter().any(|value| !value.is_finite()) {
        return Err(GeometryError::InvalidGeometry(
            "point classification requires finite coordinates".into(),
        ));
    }
    let directions = [
        unit([1.0, 0.371, 0.127])?,
        unit([0.193, 1.0, 0.419])?,
        unit([0.311, 0.233, 1.0])?,
    ];
    for direction in directions {
        let result = classify_with_ray(brep, point, direction, brep.topology.faces.iter())?;
        if result != PointClassification::Unknown {
            return Ok(result);
        }
    }
    Ok(PointClassification::Unknown)
}

pub(crate) fn classify_point_in_shell(
    brep: &BrepEnvelope,
    faces: &[u32],
    point: Point3,
) -> Result<PointClassification, GeometryError> {
    for direction in [
        unit([1.0, 0.371, 0.127])?,
        unit([0.193, 1.0, 0.419])?,
        unit([0.311, 0.233, 1.0])?,
    ] {
        let result = classify_with_ray(
            brep,
            point,
            direction,
            faces
                .iter()
                .map(|face| &brep.topology.faces[*face as usize]),
        )?;
        if result != PointClassification::Unknown {
            return Ok(result);
        }
    }
    Ok(PointClassification::Unknown)
}

/// Checks that every face's material normal (surface normal times `sense`)
/// points out of the solid. Loop winding in this kernel doesn't encode the
/// material side, so `validate()` can't see a face whose `sense` is flipped;
/// this probes it geometrically instead. Per face, a point just in front of
/// the face must not be inside while the point just behind it is outside.
pub fn validate_orientation(brep: &BrepEnvelope) -> Result<(), GeometryError> {
    brep.validate()?;
    validate_orientation_validated(brep)
}

pub(crate) fn validate_orientation_validated(brep: &BrepEnvelope) -> Result<(), GeometryError> {
    if brep.solids.is_empty() {
        return Ok(());
    }
    let Some(bounds) = brep.bounds_unchecked()? else {
        return Ok(());
    };
    let diagonal = norm(std::array::from_fn(|i| bounds.axes[i].width()));
    // Below ~1e-12 of the model size the tolerance is at the edge of f64
    // precision and point classification stops being reliable.
    if brep.accuracy.geometric < 1e-12 * diagonal {
        return Ok(());
    }
    let deflection = (0.02 * diagonal).max(4.0 * brep.accuracy.geometric);
    // A body the tessellator can't mesh (extreme scales, precision limits)
    // can't be probed either; that is reported by tessellation itself.
    let mesh = match super::tessellation::tessellate(brep, deflection, 1 << 22) {
        Ok(mesh) => mesh,
        Err(GeometryError::LimitExceeded(_) | GeometryError::UnresolvedTessellation(_)) => {
            return Ok(())
        }
        Err(error) => return Err(error),
    };

    // One probe point per face: the centroid of its largest triangle.
    let mut probes: Vec<Option<(f64, Point3)>> = vec![None; brep.topology.faces.len()];
    for (triangle, ids) in mesh.indices.chunks_exact(3).enumerate() {
        let p: [Point3; 3] = std::array::from_fn(|k| {
            std::array::from_fn(|j| mesh.positions[ids[k] as usize * 3 + j])
        });
        let area = norm(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        let probe = &mut probes[mesh.triangle_face_ids[triangle] as usize];
        if probe.is_none_or(|(best, _)| area > best) {
            *probe = Some((
                area,
                std::array::from_fn(|j| (p[0][j] + p[1][j] + p[2][j]) / 3.0),
            ));
        }
    }

    let offset = 100.0 * brep.accuracy.geometric;
    for (face, probe) in brep.topology.faces.iter().zip(probes) {
        let Some((_, centroid)) = probe else {
            continue;
        };
        let surface = brep.geometry.surface(face.surface)?;
        let uv = surface.project(centroid, None)?;
        let Ok(normal) = surface.normal_at(uv) else {
            continue;
        };
        let step = scale(unit(normal)?, face.sense.multiplier() * offset);
        let on_face = surface.point_at(uv)?;
        let front = classify_point_validated(brep, add(on_face, step))?;
        let back = classify_point_validated(brep, sub(on_face, step))?;
        if front == PointClassification::Inside && back == PointClassification::Outside {
            return Err(GeometryError::InvalidTopology(format!(
                "face {} ('{}') has its material normal pointing into the solid",
                face.id, face.key
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analytic::{
        booleans::{boolean_brep, BooleanOp},
        primitives,
        topology::Accuracy,
        Frame3,
    };

    fn accuracy() -> Accuracy {
        Accuracy {
            geometric: 1e-8,
            intersection: 1e-9,
            tessellation: 0.01,
            exchange: 1e-5,
        }
    }

    #[test]
    fn analytic_ray_classification_covers_all_surface_families() {
        let bodies = [
            (
                primitives::cuboid("box".into(), Frame3::IDENTITY, [2.0; 3], accuracy()).unwrap(),
                [1.0, 1.0, 1.0],
                [3.0, 3.0, 3.0],
            ),
            (
                primitives::sphere("sphere".into(), Frame3::IDENTITY, 1.0, accuracy()).unwrap(),
                [0.0; 3],
                [2.0, 0.0, 0.0],
            ),
            (
                primitives::cylinder("cylinder".into(), Frame3::IDENTITY, 1.0, 2.0, accuracy())
                    .unwrap(),
                [0.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
            ),
            (
                primitives::cone("cone".into(), Frame3::IDENTITY, 1.0, 2.0, accuracy()).unwrap(),
                [0.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
            ),
            (
                primitives::torus("torus".into(), Frame3::IDENTITY, 3.0, 1.0, accuracy()).unwrap(),
                [3.0, 0.0, 0.0],
                [0.0; 3],
            ),
        ];
        for (body, inside, outside) in bodies {
            assert_eq!(
                classify_point(&body, inside).unwrap(),
                PointClassification::Inside,
                "{}",
                body.id
            );
            assert_eq!(
                classify_point(&body, outside).unwrap(),
                PointClassification::Outside,
                "{}",
                body.id
            );
        }
    }

    #[test]
    fn cavity_shell_uses_regularized_material_occupancy() {
        let host =
            primitives::torus("host".into(), Frame3::IDENTITY, 3.0, 1.0, accuracy()).unwrap();
        let cutter =
            primitives::torus("cutter".into(), Frame3::IDENTITY, 3.0, 0.4, accuracy()).unwrap();
        let result = boolean_brep(&host, &cutter, BooleanOp::Subtraction, "shell".into()).unwrap();
        assert_eq!(
            classify_point(&result.brep, [3.0, 0.0, 0.0]).unwrap(),
            PointClassification::Outside
        );
        assert_eq!(
            classify_point(&result.brep, [3.7, 0.0, 0.0]).unwrap(),
            PointClassification::Inside
        );
        assert_eq!(
            classify_point(&result.brep, [0.0; 3]).unwrap(),
            PointClassification::Outside
        );
    }
}
