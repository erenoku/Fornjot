use fj_math::{Point, Scalar};

use crate::objects::{GlobalVertex, SurfaceVertex};

use super::{Validate, ValidationConfig, ValidationError};

impl Validate for SurfaceVertex {
    fn validate_with_config(
        &self,
        config: &ValidationConfig,
        errors: &mut Vec<ValidationError>,
    ) {
        SurfaceVertexValidationError::check_position(self, config, errors);
    }
}

impl Validate for GlobalVertex {
    fn validate_with_config(
        &self,
        _: &ValidationConfig,
        _: &mut Vec<ValidationError>,
    ) {
    }
}

/// [`SurfaceVertex`] validation error
#[derive(Clone, Debug, thiserror::Error)]
pub enum SurfaceVertexValidationError {
    /// Mismatch between position and position of global form
    #[error(
        "`SurfaceVertex` position doesn't match position of its global form\n\
        - Surface position: {surface_position:?}\n\
        - Surface position converted to global position: \
            {surface_position_as_global:?}\n\
        - Global position: {global_position:?}\n\
        - Distance between the positions: {distance}\n\
        - `SurfaceVertex`: {surface_vertex:#?}"
    )]
    PositionMismatch {
        /// The position of the surface vertex
        surface_position: Point<2>,

        /// The surface position converted into a global position
        surface_position_as_global: Point<3>,

        /// The position of the global vertex
        global_position: Point<3>,

        /// The distance between the positions
        distance: Scalar,

        /// The surface vertex
        surface_vertex: SurfaceVertex,
    },
}

impl SurfaceVertexValidationError {
    fn check_position(
        surface_vertex: &SurfaceVertex,
        config: &ValidationConfig,
        errors: &mut Vec<ValidationError>,
    ) {
        let surface_position_as_global = surface_vertex
            .surface()
            .geometry()
            .point_from_surface_coords(surface_vertex.position());
        let global_position = surface_vertex.global_form().position();

        let distance = surface_position_as_global.distance_to(&global_position);

        if distance > config.identical_max_distance {
            errors.push(
                Box::new(Self::PositionMismatch {
                    surface_position: surface_vertex.position(),
                    surface_position_as_global,
                    global_position,
                    distance,
                    surface_vertex: surface_vertex.clone(),
                })
                .into(),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use fj_math::Point;

    use crate::{
        insert::Insert,
        objects::{GlobalVertex, SurfaceVertex},
        partial::{
            Partial, PartialGlobalVertex, PartialObject, PartialSurfaceVertex,
        },
        services::Services,
        validate::Validate,
    };

    #[test]
    fn surface_vertex_position_mismatch() -> anyhow::Result<()> {
        let mut services = Services::new();

        let valid = PartialSurfaceVertex {
            position: Some([0., 0.].into()),
            surface: Partial::from(services.objects.surfaces.xy_plane()),
            global_form: Partial::from_partial(PartialGlobalVertex {
                position: Some(Point::from([0., 0., 0.])),
            }),
        }
        .build(&mut services.objects);
        let invalid = SurfaceVertex::new(
            valid.position(),
            valid.surface().clone(),
            GlobalVertex::new([1., 0., 0.]).insert(&mut services.objects),
        );

        valid.validate_and_return_first_error()?;
        assert!(invalid.validate_and_return_first_error().is_err());

        Ok(())
    }
}
