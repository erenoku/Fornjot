use fj_math::{Point, Scalar, Transform, Triangle, Vector};

use crate::{
    iter::ObjectIters,
    objects::{Curve, Cycle, Edge, Face, Surface, Vertex, VerticesOfEdge},
    shape::{LocalForm, Shape},
};

use super::{CycleApprox, Tolerance};

/// Create a solid by sweeping a sketch
pub fn sweep(
    source: Shape,
    path: impl Into<Vector<3>>,
    tolerance: Tolerance,
    color: [u8; 4],
) -> Shape {
    let path = path.into();

    let is_sweep_along_negative_direction =
        path.dot(&Vector::from([0., 0., 1.])) < Scalar::ZERO;

    let mut target = Shape::new();

    for edge in source.edge_iter() {
        if let Some(vertices) = edge.vertices() {
            create_non_continuous_side_face(
                path,
                is_sweep_along_negative_direction,
                vertices,
                color,
                &mut target,
            );
            continue;
        }

        create_continuous_side_face(edge, path, tolerance, color, &mut target);
    }

    target
}

fn create_non_continuous_side_face(
    path: Vector<3>,
    is_sweep_along_negative_direction: bool,
    vertices_bottom: [Vertex; 2],
    color: [u8; 4],
    target: &mut Shape,
) {
    let vertices = {
        let vertices_top = vertices_bottom.map(|vertex| {
            let point = vertex.point + path;
            Vertex { point }
        });

        let [[a, b], [c, d]] = [vertices_bottom, vertices_top];

        let vertices = if is_sweep_along_negative_direction {
            [b, a, c, d]
        } else {
            [a, b, d, c]
        };

        vertices.map(|vertex| target.get_handle_or_insert(vertex))
    };

    let surface = {
        let [a, b, _, c] = vertices.clone().map(|vertex| vertex.get().point);
        Surface::plane_from_points([a, b, c])
    };
    let surface = target.get_handle_or_insert(surface);

    let cycle = {
        let [a, b, c, d] = vertices;

        let mut vertices =
            vec![([0., 0.], a), ([1., 0.], b), ([1., 1.], c), ([0., 1.], d)];
        if let Some(vertex) = vertices.first().cloned() {
            vertices.push(vertex);
        }

        let mut edges = Vec::new();
        for vertices in vertices.windows(2) {
            // Can't panic, as we passed `2` to `windows`.
            //
            // Can be cleaned up, once `array_windows` is stable"
            // https://doc.rust-lang.org/std/primitive.slice.html#method.array_windows
            let [a, b] = [&vertices[0], &vertices[1]];

            let curve = {
                let local = Curve::line_from_points([a.0, b.0]);

                let global = [a, b].map(|vertex| vertex.1.get().point);
                let global = Curve::line_from_points(global);
                let global = target.get_handle_or_insert(global);

                LocalForm::new(local, global)
            };

            let vertices = VerticesOfEdge::from_vertices([
                LocalForm::new(Point::from([0.]), a.1.clone()),
                LocalForm::new(Point::from([1.]), b.1.clone()),
            ]);

            let edge = {
                let local = Edge {
                    curve: curve.clone(),
                    vertices: vertices.clone(),
                };

                let global = Edge {
                    curve: LocalForm::canonical_only(curve.canonical()),
                    vertices,
                };
                let global = target.get_handle_or_insert(global);

                LocalForm::new(local, global)
            };

            edges.push(edge);
        }

        let cycle = {
            let local = Cycle { edges };

            let global =
                Cycle::new(local.edges.iter().map(|edge| edge.canonical()));
            let global = target.get_handle_or_insert(global);

            LocalForm::new(local, global)
        };

        cycle
    };

    let face = Face::new(surface, [cycle], [], color);
    target.get_handle_or_insert(face);
}

fn create_continuous_side_face(
    edge: Edge<3>,
    path: Vector<3>,
    tolerance: Tolerance,
    color: [u8; 4],
    target: &mut Shape,
) {
    let translation = Transform::translation(path);

    let mut tmp = Shape::new();
    let edge = tmp.merge(edge);
    let cycle = Cycle::new(vec![edge]);
    let approx = CycleApprox::new(&cycle, tolerance);

    let mut quads = Vec::new();
    for segment in approx.segments() {
        let [v0, v1] = segment.points();
        let [v3, v2] = {
            let segment = translation.transform_segment(&segment);
            segment.points()
        };

        quads.push([v0, v1, v2, v3]);
    }

    let mut side_face: Vec<(Triangle<3>, _)> = Vec::new();
    for [v0, v1, v2, v3] in quads {
        side_face.push(([v0, v1, v2].into(), color));
        side_face.push(([v0, v2, v3].into(), color));
    }

    target.insert(Face::Triangles(side_face));
}

#[cfg(test)]
mod tests {
    use fj_math::{Scalar, Vector};

    use crate::{
        algorithms::Tolerance,
        iter::ObjectIters,
        objects::{Face, Surface},
        shape::Shape,
    };

    use super::sweep;

    #[test]
    fn side_faces_positive() -> anyhow::Result<()> {
        test_side_faces(
            [0., 0., 1.],
            [
                Surface::plane_from_points([
                    [0., 0., 0.],
                    [1., 0., 0.],
                    [0., 0., 1.],
                ]),
                Surface::plane_from_points([
                    [1., 0., 0.],
                    [0., 1., 0.],
                    [1., 0., 1.],
                ]),
                Surface::plane_from_points([
                    [0., 1., 0.],
                    [0., 0., 0.],
                    [0., 1., 1.],
                ]),
            ],
        )
    }

    #[test]
    fn side_faces_negative() -> anyhow::Result<()> {
        test_side_faces(
            [0., 0., -1.],
            [
                Surface::plane_from_points([
                    [0., 0., 0.],
                    [0., 1., 0.],
                    [0., 0., -1.],
                ]),
                Surface::plane_from_points([
                    [0., 1., 0.],
                    [1., 0., 0.],
                    [0., 1., -1.],
                ]),
                Surface::plane_from_points([
                    [1., 0., 0.],
                    [0., 0., 0.],
                    [1., 0., -1.],
                ]),
            ],
        )
    }

    fn test_side_faces(
        direction: impl Into<Vector<3>>,
        expected_surfaces: [Surface; 3],
    ) -> anyhow::Result<()> {
        let tolerance = Tolerance::from_scalar(Scalar::ONE)?;

        let mut shape = Shape::new();

        let surface = Surface::xy_plane();
        let _sketch = Face::builder(surface, &mut shape)
            .with_exterior_polygon([[0., 0.], [1., 0.], [0., 1.]])
            .build();

        let solid = sweep(shape, direction, tolerance, [255, 0, 0, 255]);

        let mut shape = Shape::new();
        let faces = expected_surfaces.map(|surface| {
            Face::builder(surface, &mut shape)
                .with_exterior_polygon([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
                .build()
                .get()
        });

        for face in faces {
            assert!(solid.face_iter().any(|f| f == face));
        }

        Ok(())
    }

    // TASK: Bottom face, positive direction.
    // TASK: Bottom face, negative direction.
    // TASK: Top face, positive direction.
    // TASK: Top face, negative direction.
}
