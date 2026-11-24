use opengeometry::world_graph::{EditScope, ModifyingOperation, WorldGraph};
use opengeometry_test_support::scenes::acceptance::recut_scene;
use opengeometry_test_support::scenes::prismatic::box_and_leg;
use opengeometry_test_support::world_graph::assert_cut_keeps_volume;

fn assert_cut_by_general(world: &mut WorldGraph, target: &str, tool: &str) {
    assert_cut_keeps_volume(world, target, tool);
    world
        .operate(
            target,
            ModifyingOperation::Subtract,
            &[tool.to_string()],
            EditScope::Node,
        )
        .unwrap();
    assert_eq!(
        world.report(target).unwrap().unwrap()["handlers"],
        serde_json::json!(["general_boolean"])
    );
}

#[test]
fn rail_minus_wall_cuts_through_the_general_boolean_and_keeps_volume() {
    let mut scene = recut_scene();
    let (rail, wall) = (scene.rail.clone(), scene.wall.clone());
    assert_cut_by_general(&mut scene.graph, &rail, &wall);
}

#[test]
fn cuboid_minus_oblique_cylinder_cuts_through_the_general_boolean_and_keeps_volume() {
    let mut world = box_and_leg(10.0, [0.5, -0.1, 0.0]);
    assert_cut_by_general(&mut world, "box", "leg");
}
