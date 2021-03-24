
#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Point2D {
    pub x: f32,
    pub y: f32,
}

impl Point2D {
    pub fn new(x: f32, y: f32) -> Point2D {
        Point2D { x, y }
    }
    pub fn plus(&self, other: &Point2D) -> Point2D {
        Point2D::new(self.x + other.x, self.y + other.y)
    }
    pub fn minus(&self, other: &Point2D) -> Point2D {
        Point2D::new(self.x - other.x, self.y - other.y)
    }
    pub fn translate(&self, direction: &Vector2D) -> Point2D {
        Point2D::new(self.x + direction.x, self.y + direction.y)
    }
    pub fn midpoint(&self, other: &Point2D) -> Point2D {
        let ab = Vector2D::join(&self, &other);
        self.translate(&ab.scale(0.5))
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Vector2D {
    pub x: f32,
    pub y: f32,
}

impl Vector2D {
    /// Create a Vector2D from it's x and y components 
    pub fn new(x: f32, y: f32) -> Vector2D {
        Vector2D { x, y }
    }
    /// Determine the Vector2D that connects two Point2Ds
    pub fn join(a: &Point2D, b: &Point2D) -> Vector2D {
        Vector2D::new(b.x - a.x, b.y - a.y)
    }
    /// Convert a Point2D into a Vector2D
    pub fn from_point(p: &Point2D) -> Vector2D {
        Vector2D::new(p.x, p.y)
    }
    /// Return the same Vector2D, but pointing in the opposite direction
    pub fn reverse(&self) -> Vector2D {
        Vector2D::new(-1.0 * self.x, -1.0 * self.y)
    }
    /// Calculate the cross product of self and another Vector2D
    pub fn cross(&self, other: &Vector2D) -> f32 {
        self.x * other.y - self.y * other.x
    }
    /// Calculate the sum of self and another Vector2D
    pub fn plus(&self, other: &Vector2D) -> Vector2D {
        Vector2D::new(self.x + other.x, self.y + other.y)
    }
    /// Calculate the difference of self and another Vector2D
    pub fn minus(&self, other: &Vector2D) -> Vector2D {
        Vector2D::new(self.x - other.x, self.y - other.y)
    }
    /// Calculate the dot product of self and another Vector2D
    pub fn dot(&self, other: &Vector2D) -> f32 {
        self.x * other.x + self.y * other.y
    }
    /// Multiply self by a scalar value, modifying the length
    pub fn scale(&self, scalar: f32) -> Vector2D {
        Vector2D::new(self.x * scalar, self.y * scalar)
    }
    /// Determine the length
    pub fn magnitude(&self) -> f32 {
        (self.x.powf(2.0) + self.y.powf(2.0)).sqrt()
    }
    /// Return a vector in the same direction, but with a length of 1.0
    pub fn normalize(&self) -> Vector2D {
        let magnitude = self.magnitude();
        self.scale(1.0 / magnitude)
    }
    /// Return a 90 degree clockwise rotation of self
    pub fn orthagonal(&self) -> Vector2D {
        Vector2D::new(-self.y, self.x)
    }
    /// Return self, or the opposite of self, whichever is pointing more
    /// strongly in the direction of q.
    pub fn facing(&self, q: &Point2D) -> Vector2D {
        match self.dot(&Vector2D::new(q.x, q.y)) >= 0.0 {
            false => self.reverse(),
            true => self.clone(),
        }
    }
    /// Reflect self across another Vector2D
    pub fn reflect(&self, other: &Vector2D) -> Vector2D {
        let normal_other = &other.normalize();
        let dp = self.dot(&normal_other);
        let dpx2 = 2.0 * dp;
        let v3 = normal_other.scale(dpx2);
        let r = self.minus(&v3);
        r
    }

}


#[test]
fn point2d_midpoint() {
    let a = Point2D::new(2.0, 1.0);
    let b = Point2D::new(-3.0, 3.0);
    let midpoint = a.midpoint(&b);
    assert_eq!(midpoint, Point2D::new(-0.5, 2.0))
}

#[test]
fn vector_orthagonal_facing() {
    let q = Point2D::new(-3.0, 2.0);
    let a = Point2D::new(-8.0, 3.0);
    let b = Point2D::new(-7.0, 7.0);

    let d = Vector2D::join(&a, &b).orthagonal().facing(&q.minus(&a));
    assert_eq!(d, Vector2D::new(4.0, -1.0));

    // Expect the result to be inverted when q is on the other side of ab
    let q = Point2D::new(-9.0, 6.0);
    let d = Vector2D::join(&a, &b).orthagonal().facing(&q.minus(&a));
    assert_eq!(d, Vector2D::new(-4.0, 1.0));
}

#[derive(Debug, PartialEq, Clone)]
pub struct ConvexPolygon {
    pub vertices: Vec<Point2D>,
}

impl ConvexPolygon {
    pub fn new(vertices: Vec<Point2D>) -> ConvexPolygon {
        ConvexPolygon { vertices }
    }
}
