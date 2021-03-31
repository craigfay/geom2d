
use super::primitives::Vector2D;
use super::primitives::Point2D;
use super::primitives::ConvexPolygon;


pub struct RigidBody2D {
    pub polygon: ConvexPolygon,
    pub position: Point2D,
    pub linear: Vector2D,
    pub angular: f32,
}

impl RigidBody2D {
    pub fn new() -> RigidBody2D {
        RigidBody2D {
            polygon: ConvexPolygon::new(vec![]),
            position: Point2D::new(0.0, 0.0),
            linear: Vector2D::new(0.0, 0.0),
            angular: 0.0,
        }
    }
    pub fn vertex(mut self, x: f32, y: f32) -> RigidBody2D {
        self.polygon.vertices.push(Point2D::new(x, y));
        self
    }
    pub fn position(mut self, x: f32, y: f32) -> RigidBody2D {
        self.position = Point2D::new(x, y);
        self
    }
    pub fn linear(mut self, x: f32, y: f32) -> RigidBody2D {
        self.linear = Vector2D::new(x, y);
        self
    }
    pub fn angular(mut self, velocity: f32) -> RigidBody2D {
        self.angular = velocity;
        self
    }
}


#[test]
fn construct_rigid_body_2d() {
    let rb = RigidBody2D::new()
        .vertex(2.0, 2.0)
        .vertex(4.0, 2.0)
        .vertex(4.0, 4.0)
        .vertex(4.0, 2.0)
        .position(3.0, 3.0)
        .linear(1.0, 1.0)
        .angular(1.0);
}
