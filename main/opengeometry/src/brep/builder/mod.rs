mod boundary;
mod brep_builder;
#[cfg(test)]
mod fixtures;
#[cfg(test)]
mod tests;
mod wire;

pub use boundary::{boundary, plane_boundary, projected_boundary, uv_conic, uv_line, UvCurve};
pub(crate) use boundary::{dimensions, padded_uv_bounds};
pub use brep_builder::{Builder, Use};
#[cfg(test)]
pub(crate) use fixtures::{wrapped_sector, SECTOR};
pub(crate) use wire::WireUse;
