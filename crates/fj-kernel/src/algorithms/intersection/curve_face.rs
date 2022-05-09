use fj_math::Scalar;

use crate::{geometry::Curve, topology::Face};

/// Determine the intersection between a [`Curve`] and a [`Face`]
pub fn curve_face(curve: &Curve, face: &Face) -> Vec<[Scalar; 2]> {
    let line = match curve {
        Curve::Line(line) => line,
        _ => todo!("Curve-face intersection only supports lines"),
    };

    let face_as_polygon = face
        .exteriors()
        .chain(face.interiors())
        .flat_map(|cycle| {
            let edges: Vec<_> = cycle.edges().collect();
            edges
        })
        .map(|edge| {
            let line = match edge.curve() {
                Curve::Line(line) => line,
                _ => todo!("Curve-face intersection only supports polygons"),
            };

            let vertices = match edge.vertices() {
                Some(vertices) => vertices,
                None => todo!(
                    "Curve-face intersection does not support faces with \
                    continuous edges"
                ),
            };

            (line, vertices)
        });

    for (edge, vertices) in face_as_polygon {
        // TASK: Convert `line` to `Line<2>`.
        // TASK: Convert `edge` and `vertices` to `Segment<2>`.
        // TASK: Determine intersection and store it in a `Vec`.

        // TASK: Implement.
        let _ = edge;
        let _ = vertices;
    }

    // TASK: Assert that there is an even number of intersection points.
    // TASK: Sort intersection points, pair them, return.

    // TASK: Implement.
    let _ = line;
    todo!()
}

#[cfg(test)]
mod tests {
    use fj_math::{Line, Point, Scalar, Vector};

    use crate::{
        geometry::{Curve, Surface},
        shape::Shape,
        topology::Face,
    };

    #[test]
    #[ignore]
    fn curve_face() -> anyhow::Result<()> {
        let mut shape = Shape::new();

        let curve = Curve::Line(Line {
            origin: Point::from([-3., 0., 0.]),
            direction: Vector::from([1., 0., 0.]),
        });

        #[rustfmt::skip]
        let exterior = [
            [-2., -2., 0.],
            [ 2., -2., 0.],
            [ 2.,  2., 0.],
            [-2.,  2., 0.],
        ];
        #[rustfmt::skip]
        let interior = [
            [-1., -1., 0.],
            [ 1., -1., 0.],
            [ 1.,  1., 0.],
            [-1.,  1., 0.],
        ];

        let face = Face::builder(Surface::xy_plane(), &mut shape)
            .with_exterior_polygon(exterior)
            .with_interior_polygon(interior)
            .build()?
            .get();

        let expected: Vec<_> = [[1., 2.], [4., 5.]]
            .into_iter()
            .map(|interval: [f64; 2]| interval.map(Scalar::from))
            .collect();
        assert_eq!(super::curve_face(&curve, &face), expected);

        Ok(())
    }
}
