use crate::accuracy::accuracy;
use opengeometry::tessellation::Measure;
use opengeometry::world_graph::{CreateOptions, EditScope, ModifyingOperation, WorldGraph};

pub fn graph() -> WorldGraph {
    WorldGraph::new(accuracy(1e-6)).unwrap()
}

pub fn named(id: &str) -> CreateOptions {
    CreateOptions {
        og_id: Some(id.into()),
        ..CreateOptions::default()
    }
}

fn trial(world: &mut WorldGraph, target: &str, operation: ModifyingOperation, tool: &str) -> Measure {
    let mark = world.mark().unwrap();
    world
        .operate(target, operation, &[tool.to_string()], EditScope::Node)
        .unwrap();
    let measured = world.measure(target, Some(1e-4)).unwrap();
    world.rollback(mark).unwrap();
    world.release(mark).unwrap();
    measured
}

pub fn assert_cut_keeps_volume(world: &mut WorldGraph, target: &str, tool: &str) {
    let whole = world.measure(target, Some(1e-4)).unwrap();
    let rest = trial(world, target, ModifyingOperation::Subtract, tool);
    let common = trial(world, target, ModifyingOperation::Intersect, tool);
    let missing = whole.volume - rest.volume - common.volume;
    let smaller = rest.volume.min(common.volume);
    assert!(
        missing.abs() <= 0.01 * smaller,
        "{target} {} - {tool} cut {} + common {}: {missing} over 1% of {smaller}",
        whole.volume,
        rest.volume,
        common.volume
    );
}
