
use super::primitives::Vector2D;
use super::primitives::Point2D;
use super::primitives::ConvexPolygon;

pub struct Position2D {
    x: Point2D,
    y: Point2D,
}

pub struct RigidBody2D {
    pub polygon: ConvexPolygon,
    pub position: Position2D,
    pub linear_vel: Vector2D,
    pub angular_vel: f32,
}
