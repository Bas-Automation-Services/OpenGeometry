use super::curve::{meeting, CurveRef};
use super::flat::Flat;
use super::pieces::{Cut, Operand, Place};
use crate::brep::GeometryError;
use crate::math::{norm, sub, Point3};

pub(super) fn imprint(
    operands: &mut [Operand<'_>; 2],
    faces: [u32; 2],
) -> Result<(), GeometryError> {
    for from in 0..2 {
        let onto = 1 - from;
        for (edge, curve, lo, hi) in operands[from].face_edges(faces[from])? {
            let mut ts = vec![lo, hi];
            for p in crossings(operands, from, faces, curve, lo, hi)? {
                let t = curve.parameter(p, 0.5 * (lo + hi))?;
                if t > lo && t < hi {
                    ts.push(t);
                }
                operands[onto].split_at(faces[onto], p)?;
                operands[from].split_edge_at(edge, p)?;
            }
            ts.sort_by(f64::total_cmp);
            ts.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
            for w in ts.windows(2) {
                let mid = curve.point(0.5 * (w[0] + w[1]))?;
                match operands[onto].locate(faces[onto], mid)? {
                    Place::Outside => continue,
                    Place::Inside => operands[onto].add_cut(
                        faces[onto],
                        Cut {
                            curve,
                            lo: w[0],
                            hi: w[1],
                        },
                    ),
                    Place::Boundary => {}
                }
                for t in [w[0], w[1]] {
                    let p = curve.point(t)?;
                    operands[onto].split_at(faces[onto], p)?;
                    operands[from].split_edge_at(edge, p)?;
                }
            }
        }
    }
    Ok(())
}

fn crossings(
    operands: &[Operand<'_>; 2],
    from: usize,
    faces: [u32; 2],
    curve: CurveRef<'_>,
    lo: f64,
    hi: f64,
) -> Result<Vec<Point3>, GeometryError> {
    let onto = 1 - from;
    let snap = operands[onto].snap();
    let map = &operands[onto].maps[faces[onto] as usize];
    let (stations, flat) = map.trace(curve, lo, hi)?;
    let mut points = Vec::new();
    for (edge, other, olo, ohi) in operands[onto].face_edges(faces[onto])? {
        let (ostations, oflat) = operands[onto].edge_flat(faces[onto], edge)?;
        for i in 0..flat.len() - 1 {
            for j in 0..oflat.len() - 1 {
                if let Some((s, r)) = segment_crossing(flat[i], flat[i + 1], oflat[j], oflat[j + 1])
                {
                    let t = stations[i] + s * (stations[i + 1] - stations[i]);
                    let u = ostations[j] + r * (ostations[j + 1] - ostations[j]);
                    points.push(meeting(&[(curve, t), (other, u)], curve.point(t)?)?);
                }
            }
        }
        for (c, a, b, d) in [(other, olo, ohi, curve), (curve, lo, hi, other)] {
            for t in [a, b] {
                let p = c.point(t)?;
                let s = d.parameter(p, t)?;
                if norm(sub(d.point(s)?, p)) <= snap {
                    points.push(p);
                }
            }
        }
    }
    let mut unique: Vec<Point3> = Vec::new();
    for p in points {
        if unique.iter().all(|q| norm(sub(*q, p)) > snap) {
            unique.push(p);
        }
    }
    Ok(unique)
}

fn segment_crossing(a: Flat, b: Flat, c: Flat, d: Flat) -> Option<(f64, f64)> {
    let r = [b[0] - a[0], b[1] - a[1]];
    let s = [d[0] - c[0], d[1] - c[1]];
    let denominator = r[0] * s[1] - r[1] * s[0];
    if denominator.abs() <= 1e-300 {
        return None;
    }
    let q = [c[0] - a[0], c[1] - a[1]];
    let t = (q[0] * s[1] - q[1] * s[0]) / denominator;
    let u = (q[0] * r[1] - q[1] * r[0]) / denominator;
    ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then_some((t, u))
}
