use crate::{
    kernel::{
        geometry::{Circle, Curve, Line},
        shape::handle::Handle,
    },
    math::{Point, Scalar, Vector},
};

use super::vertices::Vertex;

/// A cycle of connected edges
///
/// The end of each edge in the cycle must connect to the beginning of the next
/// edge. The end of the last edge must connect to the beginning of the first
/// one.
#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct Cycle {
    pub edges: Vec<Handle<Edge>>,
}

/// An edge of a shape
#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct Edge {
    /// Access the curve that defines the edge's geometry
    ///
    /// The edge can be a segment of the curve that is bounded by two vertices,
    /// or if the curve is continuous (i.e. connects to itself), the edge could
    /// be defined by the whole curve, and have no bounding vertices.
    pub curve: Curve,

    /// Access the vertices that bound the edge on the curve
    ///
    /// If there are no such vertices, that means that both the curve and the
    /// edge are continuous (i.e. connected to themselves).
    ///
    /// # Implementation note
    ///
    /// Since these vertices bound the edge, they must lie on the curve. This
    /// isn't enforced at all, however. It would make sense to store 1D vertices
    /// here, and indeed, this was the case in the past.
    ///
    /// It got in the way of some work, however, so it made sense to simplify
    /// it by storing 3D vertices. It will probably make sense to revert this
    /// and store 1D vertices again, at some point.
    pub vertices: Option<[Handle<Vertex>; 2]>,
}

impl Edge {
    /// Create a line segment
    pub fn line_segment(vertices: [Handle<Vertex>; 2]) -> Self {
        Edge {
            curve: Curve::Line(Line::from_points(
                vertices.clone().map(|vertex| vertex.point()),
            )),
            vertices: Some(vertices),
        }
    }

    /// Create a circle
    pub fn circle(radius: Scalar) -> Self {
        Edge {
            curve: Curve::Circle(Circle {
                center: Point::origin(),
                radius: Vector::from([radius, Scalar::ZERO]),
            }),
            vertices: None,
        }
    }
}
