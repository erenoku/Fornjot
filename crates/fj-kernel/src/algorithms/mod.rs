//! Collection of algorithms that are used by the kernel
//!
//! Algorithmic code is collected in this module, to keep other modules focused
//! on their respective purpose.

mod approx;
mod sweep;
mod sweep2;
mod transform;
mod triangulation;

pub mod intersection;

pub use self::{
    approx::{CycleApprox, FaceApprox, InvalidTolerance, Tolerance},
    sweep::sweep_shape,
    sweep2::sweep,
    transform::transform_shape,
    triangulation::triangulate,
};
