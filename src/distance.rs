
use super::primitives::Vector2D;
use super::primitives::Point2D;
use super::primitives::ConvexPolygon;


/// Defines an interface for 2D primitives to calculate the distance to other
/// 2D primitives without having to use a separate method for each type.
pub trait Distance2D<T> {
    fn distance(&self, other: T) -> f32;
}
/// Allows a Point2D to calculate distance to a &Point2D
impl Distance2D<&Point2D> for Point2D {
    fn distance(&self, point: &Point2D) -> f32 {
        Vector2D::join(&self, &point).magnitude()
    }
}
/// Allows a Point2D to calculate distance to a &ConvexPolygon
impl Distance2D<&ConvexPolygon> for Point2D {
    fn distance(&self, polygon: &ConvexPolygon) -> f32 {
        GJK::polygon_to_point_distance(&polygon.vertices, &self)
    }
}
/// Allows a ConvexPolygon to calculate distance to a &Point2D
impl Distance2D<&Point2D> for ConvexPolygon {
    fn distance(&self, point: &Point2D) -> f32 {
        GJK::polygon_to_point_distance(&self.vertices, &point)
    }
}
/// Allows a ConvexPolygon to calculate distance to a &ConvexPolygon
impl Distance2D<&ConvexPolygon> for ConvexPolygon {
    fn distance(&self, polygon: &ConvexPolygon) -> f32 {
        GJK::polygon_to_polygon_distance(&self.vertices, &polygon.vertices)
    }
}

/// The 2D line segment formed between two points
struct Segment {
    a: Point2D,
    b: Point2D,
}

impl Segment {
    pub fn new(a: Point2D, b: Point2D) -> Segment {
        Segment { a, b }
    }

    // Express the position of q as a ratio of a and b
    // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
    pub fn barycentric(&self, q: &Point2D) -> (f32, f32) {
        // Defining a unit vector pointing from a to b
        let ab = Vector2D::join(&self.a, &self.b);
        let n: Vector2D = ab.normalize();

        // Finding the midpoints of the lines connecting a -> q,
        // and b -> q
        let q_minus_a = q.minus(&self.a);
        let b_minus_q = self.b.minus(&q);

        // Calculating the length of projections onto n
        let q_minus_a_proj_len = Vector2D::from_point(&q_minus_a).dot(&n);
        let b_minus_q_proj_len = Vector2D::from_point(&b_minus_q).dot(&n);

        let len_ab = ab.magnitude();

        // Calculating the barycentric weights of the query
        // point's projection onto n. The two coordinates, or
        // weights, will sum to 1.0.
        let v = q_minus_a_proj_len / len_ab;
        let u = b_minus_q_proj_len / len_ab;
        (v, u)
    }

    // Return the point on the segment that is closest to q, and a boolean
    // indicating whether that point is on an edge, false if it's a vertex.
    pub fn do_voronoi_analysis(&self, q: &Point2D) -> SegmentVoronoiAnalysis {
        let (v, u) = self.barycentric(&q);

        // When some point on ab is closer to query_point than
        // a or b themselves...
        if v > 0.0 && u > 0.0 {
            // Using the barycentric weights to calculate the closest
            // point: G = uA + vB
            let ua = Vector2D::from_point(&self.a).scale(u);
            let vb = Vector2D::from_point(&self.b).scale(v);

            return SegmentVoronoiAnalysis {
                point: Point2D::new(ua.x + vb.x, ua.y + vb.y),
                region: SegmentVoronoiRegion::AB,
            };
        }

        match v <= 0.0 {
            true => {
                SegmentVoronoiAnalysis {
                    point: self.a,
                    region: SegmentVoronoiRegion::A
                }
            }
            false => {
                SegmentVoronoiAnalysis {
                    point: self.b,
                    region: SegmentVoronoiRegion::B
                }
            }
        }
    }
}

struct SegmentVoronoiAnalysis {
    point: Point2D,
    region: SegmentVoronoiRegion,
}

#[derive(Debug, PartialEq)]
enum SegmentVoronoiRegion {
    A,
    B,
    AB,
}

impl SegmentVoronoiAnalysis {
    pub fn is_edge(&self) -> bool {
        match self.region {
            SegmentVoronoiRegion::AB => true,
            _ => false,
        }
    }
}

#[derive(Debug)]
struct TriangleVoronoiAnalysis {
    point: Point2D,
    region: TriangleVoronoiRegion,
}

impl TriangleVoronoiAnalysis {
    pub fn is_edge(&self) -> bool {
        match self.region {
            TriangleVoronoiRegion::AB => true,
            TriangleVoronoiRegion::BC => true,
            TriangleVoronoiRegion::CA => true,
            _ => false,
        }
    }

    pub fn is_contained(&self) -> bool {
        match self.region {
            TriangleVoronoiRegion::ABC => true,
            _ => false,
        }
    }
}

#[derive(Debug, PartialEq)]
enum TriangleVoronoiRegion {
    A,
    B,
    C,
    AB,
    BC,
    CA,
    ABC,
}

#[derive(Debug)]
struct Triangle {
    a: Point2D,
    b: Point2D,
    c: Point2D,
}

impl Triangle {
    pub fn new(a: Point2D, b: Point2D, c: Point2D) -> Triangle {
        Triangle { a, b, c }
    }

    // This is a variant of the shoelace formula that allows negative area,
    // depending on the winding order of ABC. This is helpful for cases when
    // we need to determine the Voronoi region of a point that's outside of ABC.
    pub fn signed_area(&self) -> f32 {
        let ab = Vector2D::join(&self.a, &self.b);
        let ac = Vector2D::join(&self.a, &self.c);
        0.5 * ab.cross(&ac)
    }

    // Express the position of q as a ratio of a, b, and c
    pub fn barycentric(&self, q: &Point2D) -> (f32, f32, f32) {
        let area_abc = self.signed_area();
        let u = Triangle::new(*q, self.b, self.c).signed_area() / area_abc;
        let v = Triangle::new(*q, self.c, self.a).signed_area() / area_abc;
        let w = Triangle::new(*q, self.a, self.b).signed_area() / area_abc;
        (u, v, w)
    }

    // Return the closest point on the triangle to q, and whether or not
    // that point lies on an edge
    pub fn do_voronoi_analysis(&self, q: &Point2D) -> TriangleVoronoiAnalysis {
        // Calculating the barycentric coordinates needed to find
        // the voronoi region of q
        let (uab, vab) = Segment::new(self.a, self.b).barycentric(&q);
        let (ubc, vbc) = Segment::new(self.b, self.c).barycentric(&q);
        let (uca, vca) = Segment::new(self.c, self.a).barycentric(&q);
        let (uabc, vabc, wabc) = self.barycentric(&q);

        // Region A
        if vca <= 0.0 && uab <= 0.0 {
            return TriangleVoronoiAnalysis {
                point: self.a,
                region: TriangleVoronoiRegion::A,
            };
        }

        // Region B
        else if vab <= 0.0 && ubc <= 0.0 {
            return TriangleVoronoiAnalysis {
                point: self.b,
                region: TriangleVoronoiRegion::B,
            };
        }

        // Region C
        else if uca <= 0.0 && vbc <= 0.0 {
            return TriangleVoronoiAnalysis {
                point: self.c,
                region: TriangleVoronoiRegion::C,
            };
        }

        // Region BC
        else if uabc <= 0.0 && ubc > 0.0 && vbc > 0.0 {
            let closest_point = Point2D {
                x: vbc * self.b.x + ubc * self.c.x,
                y: vbc * self.b.y + ubc * self.c.y,
            };
            return TriangleVoronoiAnalysis {
                point: closest_point,
                region: TriangleVoronoiRegion::BC,
            };
        }

        // Region CA
        else if vabc <= 0.0 && uca > 0.0 && vca > 0.0 {
            let closest_point = Point2D {
                x: vca * self.c.x + uca * self.a.x,
                y: vca * self.c.y + uca * self.a.y,
            };
            return TriangleVoronoiAnalysis {
                point: closest_point,
                region: TriangleVoronoiRegion::CA,
            };
        }

        // Region AB
        else if wabc <= 0.0 && uab > 0.0 && vab > 0.0 {
            let closest_point = Point2D {
                x: vab * self.a.x + uab * self.b.x,
                y: vab * self.a.y + uab * self.b.y,
            };
            return TriangleVoronoiAnalysis {
                point: closest_point,
                region: TriangleVoronoiRegion::AB,
            };
        }

        // Region ABC
        else {
            assert!(uabc > 0.0 && vabc > 0.0 && wabc > 0.0);
            return TriangleVoronoiAnalysis {
                point: q.clone(),
                region: TriangleVoronoiRegion::ABC,
            };
        }
    }
}

// A helper struct that contains helpful data that is generated between
// iterations of the GJK distance algorithm. Without using a helper struct,
// this data wouldn't be available to higher level functions.
#[derive(Debug)]
struct IterationResult {
    closest_point: Point2D,
    search_direction: Vector2D,
    simplex_contains_point: bool,
}

type SupportFn = dyn Fn(&Vector2D) -> Point2D;

// https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm
// https://www.youtube.com/watch?v=MDusDn8oTSE&ab_channel=Winterdev
// https://www.youtube.com/watch?v=Qupqu1xe7Io
// GJK is a distance algorithm for convex polygons. It presupposes the presence
// of several algorithms, implemented below.
struct GJK;
impl GJK {

    // Determine the closest point on a simplex s to a query point q,
    // a VectorD which can be used to find the next support point,
    // and whether or a 3-point simplex contained the query point
    pub fn iterate(s: &Vec<Point2D>, q: &Point2D) -> IterationResult {
        match s.len() {
            1 => {
                IterationResult {
                    closest_point: s[0],
                    search_direction: Vector2D::join(&s[0], &q),
                    simplex_contains_point: false,
                }
            },
            2 => {
                let a = s[1];
                let b = s[0];

                let voronoi_analysis = Segment::new(a, b)
                    .do_voronoi_analysis(&q);

                let search_direction = match voronoi_analysis.is_edge() {
                    true => Vector2D::join(&a, &b).orthagonal().facing(&q.minus(&a)),
                    false => Vector2D::join(&voronoi_analysis.point, &q),
                };

                IterationResult {
                    closest_point: voronoi_analysis.point,
                    simplex_contains_point: false,
                    search_direction,
                }
            },
            3 => {
                let a = s[2];
                let b = s[1];
                let c = s[0];

                let triangle = Triangle::new(a, b, c);
                let voronoi_analysis = triangle.do_voronoi_analysis(&q);

                // When the closest point is on an edge, a vector that
                // points directly at the query point is vulnerable to
                // floating point rounding errors. To prevent these types
                // of errors, we can use a vector that's perpendicular
                // to the edge, facing in the general direction of q.
                let search_direction = match voronoi_analysis.is_edge() {
                    true => Vector2D::join(&a, &b).orthagonal().facing(&q.minus(&a)),
                    false => Vector2D::join(&voronoi_analysis.point, &q),
                };

                IterationResult {
                    closest_point: voronoi_analysis.point,
                    simplex_contains_point: voronoi_analysis.is_contained(),
                    search_direction,
                }
            },
            _ => panic!("Unsupported simplex size"),
        }
    }


    pub fn polygon_to_polygon_distance(a: &Vec<Point2D>, b: &Vec<Point2D>) -> f32 {
        let support_fn = move |d| {
            // Calculating the support point on the convex hull 
            // of the minkowski difference of a and b
            let mut support_a = GJK::support(&a, d);
            let mut support_b = GJK::support(&b, d.reverse());
            let mut support_z = support_a.minus(&support_b);
            support_z
        };

        let origin = Point2D::new(0.0, 0.0);
        GJK::abstract_distance(&origin, support_fn)
    }


    pub fn polygon_to_point_distance(polygon: &Vec<Point2D>, query_point: &Point2D) -> f32 {
        let p2 = polygon.clone();

        let support_fn = move |d| { GJK::support(&p2, d) };
        GJK::abstract_distance(&query_point, support_fn)
    }


    pub fn abstract_distance<T: Fn(Vector2D) -> Point2D> (
        query_point: &Point2D,
        support: T,
    ) -> f32 {
        let arbitrary_point = support(Vector2D::new(1.0, 1.0));
        let mut simplex: Vec<Point2D> = vec![arbitrary_point];

        // Terminating if an unusual number of iterations have been done.
        // This indicates that the query point is very close to an edge or
        // vertex, and failing because of floating point imprecision.
        // What is the correct amount of iterations to allow?? 
        for _ in 0..24 {

            let iteration = GJK::iterate(&simplex, &query_point);
            let closest_point = iteration.closest_point;
            let search_direction = iteration.search_direction;
            let simplex_contains_point = iteration.simplex_contains_point;

            // Terminating early if a 3-point simplex contains query_point
            if simplex_contains_point {
                return 0.0
            }

            // Terminating early if query_point is also a point on the simplex
            if search_direction == Vector2D::new(0.0, 0.0) {
                return 0.0;
            }

            // Culling non-contributing vertices from the simplex
            while simplex.len() >= 3 {
                simplex.remove(0);
            }

            // Determining the vertex furthest in the search direction
            let support_point = support(search_direction);

            // Terminating if the support point is already in the simplex
            for point in &simplex {
                if *point == support_point {
                    return Vector2D::join(&closest_point, &query_point).magnitude();
                }
            }

            // Merging the support point with the current simplex
            simplex.push(support_point);
        }

        0.0
    }

    // Determine the point on polygon p which is furthest in the search
    // direction d.
    pub fn support(p: &Vec<Point2D>, d: Vector2D) -> Point2D {
        let mut furthest_point = &p[0];
        let mut highest_dp = Vector2D::from_point(&furthest_point).dot(&d);

        for i in 1..p.len() {
            let dp = Vector2D::from_point(&p[i]).dot(&d);
            if dp > highest_dp {
                highest_dp = dp;
                furthest_point = &p[i];
            }
        }

        furthest_point.clone()
    }
}

// Conservative Advancement can find the time of impact between two convex
// polygons, even with both linear and angular velocity. It requires a radius
// value, which is the distance of the furthest point from the centroid. This
// means that a circle centered at the polygon's centroid with the given radius
// would contain the polygon. Conservative Advancement also requires a method
// for computing the distance between the polygons.

#[test]
fn polygon_to_polygon_distance() {
    let a = vec![
        Point2D::new(2.0, 2.0),
        Point2D::new(2.0, 4.0),
        Point2D::new(4.0, 4.0),
        Point2D::new(4.0, 2.0),
    ];

    let b = vec![
        Point2D::new(-1.0, 2.0),
        Point2D::new(-4.0, 2.0),
        Point2D::new(-2.0, 4.0),
    ];

    let distance = GJK::polygon_to_polygon_distance(&a, &b);
    assert_eq!(distance, 3.0)
}

#[test]
fn segment_closest_point() {
    /*
    * * * * * * * * * * * * * * *
    * * * * * Q * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * O * A * B * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    */
    let query_point = Point2D::new(-2.0, 4.0);

    let segment = Segment {
        a: Point2D::new(4.0, 0.0),
        b: Point2D::new(2.0, 0.0),
    };

    let voronoi_analysis = segment.do_voronoi_analysis(&query_point);
    assert_eq!(voronoi_analysis.point, Point2D::new(2.0, 0.0));
    assert_eq!(voronoi_analysis.is_edge(), false);
}

#[test]
fn gjk_search_direction_2_point_simplex() {
    /*
    * * * * * * * * * * * * * * *
    * * * * * Q * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * O * A * B * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    */
    let query_point = Point2D::new(-2.0, 4.0);

    let simplex = vec![
        Point2D::new(4.0, 0.0),
        Point2D::new(2.0, 0.0),
    ];

    let iteration = GJK::iterate(&simplex, &query_point);
    assert_eq!(iteration.closest_point, Point2D::new(2.0, 0.0));
    assert_eq!(iteration.search_direction, Vector2D::new(-4.0, 4.0));
}



#[test]
fn gjk_search_direction_2_point_simplex_2() {
    /*
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * Q * * *
    * * * * * * * O * * * * * * *
    * * * * * * * * * * A P * B *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    * * * * * * * * * * * * * * *
    */
    let query_point = Point2D::new(4.0, 1.0);

    let simplex = vec![
        Point2D::new(6.0, -1.0),
        Point2D::new(3.0, -1.0),
    ];

    let iteration = GJK::iterate(&simplex, &query_point);
    assert_eq!(iteration.closest_point, Point2D::new(4.0, -1.0));
    assert_eq!(iteration.search_direction, Vector2D::new(0.0, 3.0));
}

#[test]
fn gjk_search_direction_2_point_simplex_3() {
    let query_point = Point2D::new(2.0, 0.0);

    let simplex = vec![
        Point2D::new(6.0, -1.0),
        Point2D::new(3.0, -1.0),
    ];

    let iteration = GJK::iterate(&simplex, &query_point);
    assert_eq!(iteration.closest_point, Point2D::new(3.0, -1.0));
    assert_eq!(iteration.search_direction, Vector2D::new(-1.0, 1.0));
}

#[test]
fn gjk_support_point_square() {
    let search_direction = Vector2D::new(-4.0, 4.0);

    let polygon = vec![
        Point2D::new(2.0, 2.0),
        Point2D::new(4.0, 2.0),
        Point2D::new(4.0, 0.0),
        Point2D::new(2.0, 0.0),
    ];

    let support_point = GJK::support(&polygon, search_direction);
    assert_eq!(support_point, Point2D::new(2.0, 2.0));
}

#[test]
fn area_of_right_triangle() {
    let a = Point2D::new(2.0, 1.0);
    let b = Point2D::new(2.0, 3.0);
    let c = Point2D::new(5.0, 1.0);

    let area = Triangle::new(a, b, c).signed_area();
    assert_eq!(area, -3.0);
}

#[test]
fn area_of_scalene_triangle() {
    let a = Point2D::new(-1.0, 2.0);
    let b = Point2D::new(-1.0, -3.0);
    let c = Point2D::new(-3.0, 0.0);
    let area = Triangle::new(a, b, c).signed_area();
    assert_eq!(area, -5.0);
}

#[test]
fn segment_barycentric() {
    let a = Point2D::new(0.0, 0.0);
    let b = Point2D::new(4.0, 0.0);
    let q = Point2D::new(1.0, 0.0);

    let ab = Segment::new(a, b);
    let (uab, vab) = ab.barycentric(&q);

    assert!(uab == 0.25);
    assert!(vab == 0.75);
}

#[test]
fn segment_barycentric_reversed() {
    let a = Point2D::new(0.0, 0.0);
    let b = Point2D::new(4.0, 0.0);
    let q = Point2D::new(1.0, 0.0);

    let ba = Segment::new(b, a);
    let (uba, vba) = ba.barycentric(&q);

    assert!(uba == 0.75);
    assert!(vba == 0.25);
}

#[test]
fn segment_barycentric_outside() {
    let a = Point2D::new(0.0, 0.0);
    let b = Point2D::new(4.0, 0.0);
    let q = Point2D::new(5.0, 0.0);

    let ab = Segment::new(a, b);
    let (uab, vab) = ab.barycentric(&q);

    assert!(uab == 1.25);
    assert!(vab == -0.25);
}

#[test]
fn segment_barycentric_outside_reversed() {
    let a = Point2D::new(0.0, 0.0);
    let b = Point2D::new(4.0, 0.0);
    let q = Point2D::new(5.0, 0.0);

    let ba = Segment::new(b, a);
    let (uba, vba) = ba.barycentric(&q);

    assert!(uba == -0.25);
    assert!(vba == 1.25);
}

#[test]
fn segment_voronoi_analysis_inside() {
    let a = Point2D::new(1.0, 1.0);
    let b = Point2D::new(7.0, 1.0);
    let q = Point2D::new(3.0, 3.0);

    let ab = Segment::new(a, b);
    let voronoi_analysis = ab.do_voronoi_analysis(&q);

    let expected_point = Point2D::new(3.0, 1.0);
    let difference = expected_point.minus(&voronoi_analysis.point);

    assert!(difference.x.abs() < 0.01);
    assert!(difference.y.abs() < 0.01);
    assert_eq!(voronoi_analysis.region, SegmentVoronoiRegion::AB);


    let q = Point2D::new(4.0, -1.0);
    let voronoi_analysis = ab.do_voronoi_analysis(&q);

    let expected_point = Point2D::new(4.0, 1.0);
    let difference = expected_point.minus(&voronoi_analysis.point);

    assert!(difference.x.abs() < 0.01);
    assert!(difference.y.abs() < 0.01);
    assert_eq!(voronoi_analysis.region, SegmentVoronoiRegion::AB);
}


#[test]
fn triangle_voronoi_regions() {
    let triangle = Triangle {
        a: Point2D::new(3.0, 3.0),
        b: Point2D::new(3.0, 1.0),
        c: Point2D::new(5.0, 3.0),
    };

    let q = Point2D::new(2.0, 5.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::A);

    let q = Point2D::new(2.0, 0.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::B);

    let q = Point2D::new(6.0, 7.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::C);

    let q = Point2D::new(1.0, 2.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::AB);

    let q = Point2D::new(4.0, 1.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::BC);

    let q = Point2D::new(4.0, 5.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::CA);

    let q = Point2D::new(3.5, 2.0);
    let voronoi_analysis = triangle.do_voronoi_analysis(&q);
    assert_eq!(voronoi_analysis.region, TriangleVoronoiRegion::ABC);
}

#[test]
fn triangle_voronoi_closest_points() {
    let triangle = Triangle {
        a: Point2D::new(3.0, -2.0),
        b: Point2D::new(-1.0, -2.0),
        c: Point2D::new(0.0, 0.0),
    };

    // AB
    let q = Point2D::new(2.0, -6.0);
    let va = triangle.do_voronoi_analysis(&q);
    let distance = va.point.distance(&Point2D::new(2.0, -2.0));
    assert!(distance < 0.001);


    let triangle = Triangle {
        a: Point2D::new(0.0, 0.0),
        b: Point2D::new(3.0, -2.0),
        c: Point2D::new(-1.0, -2.0),
    };

    // BC
    let q = Point2D::new(2.0, -6.0);
    let va = triangle.do_voronoi_analysis(&q);
    let distance = va.point.distance(&Point2D::new(2.0, -2.0));
    assert!(distance < 0.001);

    let triangle = Triangle {
        a: Point2D::new(-1.0, -2.0),
        b: Point2D::new(0.0, 0.0),
        c: Point2D::new(3.0, -2.0),
    };

    // CA
    let q = Point2D::new(2.0, -6.0);
    let va = triangle.do_voronoi_analysis(&q);
    let distance = va.point.distance(&Point2D::new(2.0, -2.0));
    assert!(distance < 0.001);
}

#[test]
fn all_positive_polygon_exterior_point_distance() {
    let square = vec![
        Point2D::new(5.0, 3.0),
        Point2D::new(5.0, 1.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(3.0, 3.0),
    ];

    let query_point = Point2D::new(1.0, 2.0);

    let distance = GJK::polygon_to_point_distance(&square, &query_point);
    assert_eq!(distance, 2.0);
}


#[test]
fn half_positive_polygon_interior_point_distance() {
    let triangle = vec![
        Point2D::new(-1.0, 3.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(-2.0, 1.0),
    ];

    let query_point = Point2D::new(-1.0, 2.0);

    let distance = GJK::polygon_to_point_distance(&triangle, &query_point);
    assert_eq!(distance, 0.0);
}

#[test]
fn all_negative_polygon_edge_point_distance() {
    let triangle = vec![
        Point2D::new(-1.0, -1.0),
        Point2D::new(-4.0, -4.0),
        Point2D::new(-5.0, -2.0),
    ];

    let query_point = Point2D::new(-3.0, -3.0);

    let distance = GJK::polygon_to_point_distance(&triangle, &query_point);
    assert_eq!(distance, 0.0);
}

#[test]
fn all_negative_polygon_vertex_point_distance() {
    let triangle = vec![
        Point2D::new(-1.0, -1.0),
        Point2D::new(-4.0, -4.0),
        Point2D::new(-5.0, -2.0),
    ];

    let query_point = Point2D::new(-5.0, -2.0);

    let distance = GJK::polygon_to_point_distance(&triangle, &query_point);
    assert_eq!(distance, 0.0);
}

#[test]
fn all_negative_polygon_almost_vertex_point_distance() {
    let triangle = vec![
        Point2D::new(-1.0, -1.0),
        Point2D::new(-4.0, -4.0),
        Point2D::new(-5.0, -2.0),
    ];

    let query_point = Point2D::new(-5.001, -2.002);

    let distance = GJK::polygon_to_point_distance(&triangle, &query_point);
    assert!(distance < 0.01);
}

#[test]
fn polygon_to_point_distance() {
    let polygon = ConvexPolygon::new(
        vec![
            Point2D::new(-2.0, 2.0),
            Point2D::new( -2.0, 0.0),
            Point2D::new(-4.0, 0.0),
            Point2D::new(-4.0, 2.0),
        ],
    );

    let query_point = Point2D::new(0.0, 2.0);


    let distance = polygon.distance(&query_point);
    assert_eq!(distance, 2.0);
}

