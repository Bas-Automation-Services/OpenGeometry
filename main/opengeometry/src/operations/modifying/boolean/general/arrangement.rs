use super::curve::CurveRef;
use super::flat::{signed_area, winding, Flat};
use crate::brep::GeometryError;
use std::f64::consts::TAU;

#[derive(Clone)]
pub(super) struct Use<'a> {
    pub(super) curve: CurveRef<'a>,
    pub(super) from: f64,
    pub(super) to: f64,
    pub(super) flat: Vec<Flat>,
}

impl Use<'_> {
    pub(super) fn reversed(&self) -> Self {
        let mut flat = self.flat.clone();
        flat.reverse();
        Self {
            curve: self.curve,
            from: self.to,
            to: self.from,
            flat,
        }
    }
    pub(super) fn start(&self) -> Flat {
        self.flat[0]
    }
    pub(super) fn end(&self) -> Flat {
        self.flat[self.flat.len() - 1]
    }
}

pub(super) fn reversed_loop<'a>(uses: &[Use<'a>]) -> Vec<Use<'a>> {
    uses.iter().rev().map(Use::reversed).collect()
}

pub(super) fn ring(uses: &[Use<'_>]) -> Vec<Flat> {
    uses.iter()
        .flat_map(|u| u.flat[..u.flat.len() - 1].iter().copied())
        .collect()
}

fn unresolved(message: &str) -> GeometryError {
    GeometryError::UnresolvedIntersection(format!("general boolean: {message}"))
}

pub(super) fn faces<'a>(
    uses: Vec<Use<'a>>,
    snap: f64,
) -> Result<Vec<Vec<Vec<Use<'a>>>>, GeometryError> {
    let mut nodes: Vec<Flat> = Vec::new();
    let mut node = |p: Flat| -> usize {
        if let Some(i) = nodes
            .iter()
            .position(|q| (q[0] - p[0]).hypot(q[1] - p[1]) <= snap)
        {
            return i;
        }
        nodes.push(p);
        nodes.len() - 1
    };
    let from: Vec<usize> = uses.iter().map(|u| node(u.start())).collect();
    let to: Vec<usize> = uses.iter().map(|u| node(u.end())).collect();
    let next = |h: usize| -> Option<usize> {
        let back = direction(&uses[h].flat, true);
        let mut best: Option<(usize, f64)> = None;
        for c in (0..uses.len()).filter(|&c| from[c] == to[h]) {
            let angle = direction(&uses[c].flat, false);
            let mut turn = (back - angle).rem_euclid(TAU);
            if turn < 1e-12 {
                turn = TAU;
            }
            if best.is_none_or(|(_, t)| turn < t) {
                best = Some((c, turn));
            }
        }
        best.map(|(c, _)| c)
    };
    let mut used = vec![false; uses.len()];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for first in 0..uses.len() {
        if used[first] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut h = first;
        while !used[h] {
            used[h] = true;
            cycle.push(h);
            h = next(h).ok_or_else(|| unresolved("a face's edges do not close into loops"))?;
        }
        if h != first {
            return Err(unresolved("a face's edges do not close into loops"));
        }
        cycles.push(cycle);
    }
    let rings: Vec<Vec<Flat>> = cycles
        .iter()
        .map(|c| ring(&c.iter().map(|&i| uses[i].clone()).collect::<Vec<_>>()))
        .collect();
    let areas: Vec<f64> = rings.iter().map(|r| signed_area(r)).collect();
    let mut outers: Vec<usize> = (0..cycles.len()).filter(|&i| areas[i] > 0.0).collect();
    outers.sort_by(|&a, &b| areas[a].total_cmp(&areas[b]));
    let mut faces: Vec<Vec<usize>> = outers.iter().map(|&o| vec![o]).collect();
    for hole in (0..cycles.len()).filter(|&i| areas[i] <= 0.0) {
        let probe = left_of(&uses[cycles[hole][0]].flat, snap * 4.0);
        let owner = outers
            .iter()
            .position(|&o| winding(probe, std::slice::from_ref(&rings[o])) != 0)
            .ok_or_else(|| unresolved("a hole of a face lies outside every outer loop"))?;
        faces[owner].push(hole);
    }
    Ok(faces
        .into_iter()
        .map(|loops| {
            loops
                .into_iter()
                .map(|c| cycles[c].iter().map(|&i| uses[i].clone()).collect())
                .collect()
        })
        .collect())
}

fn direction(flat: &[Flat], at_end: bool) -> f64 {
    let (a, b) = if at_end {
        (flat[flat.len() - 1], flat[flat.len() - 2])
    } else {
        (flat[0], flat[1])
    };
    (b[1] - a[1]).atan2(b[0] - a[0])
}

fn left_of(flat: &[Flat], distance: f64) -> Flat {
    let i = (flat.len() - 1) / 2;
    let (a, b) = (flat[i], flat[i + 1]);
    let d = [b[0] - a[0], b[1] - a[1]];
    let len = d[0].hypot(d[1]);
    let m = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
    [m[0] - d[1] / len * distance, m[1] + d[0] / len * distance]
}
