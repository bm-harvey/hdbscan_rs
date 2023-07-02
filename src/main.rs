#![allow(dead_code)]
#![allow(unused_variables)]

use std::cmp::min;

fn main() {
    println!("Hello, world!");

    let mut data: Vec<Point> = vec![];
    for _ in 0..10 {
        data.push(Point::from(vec![1., 2., 3.]));
    }

    let mut clust = Clusterer::new();
    clust.set_leaf_size(45);
}

struct Point {
    coordinate: Vec<f64>,
}

impl Point {
    fn new() -> Point {
        Point { coordinate: vec![] }
    }

    fn from(coordinate: Vec<f64>) -> Point {
        Point { coordinate }
    }

    fn num_dimensions(&self) -> usize {
        self.coordinate.len()
    }

    fn distance_to(&self, other: &Point) -> f64 {
        let dims = min(self.num_dimensions(), other.num_dimensions());

        (0..dims)
            .map(|idx| self.coordinate[idx] - other.coordinate[idx])
            .map(|diff| diff.powi(2))
            .sum::<f64>()
            .sqrt()
    }

    fn distance_to_sqr(&self, other: &Point) -> f64 {
        let dims = min(self.num_dimensions(), other.num_dimensions());

        (0..dims)
            .map(|idx| self.coordinate[idx] - other.coordinate[idx])
            .map(|diff| diff.powi(2))
            .sum::<f64>()
    }
}

struct ClusteredPoint<'a> {
    point: &'a Point,
    cluster_id: u16,
}

impl<'a> ClusteredPoint<'a> {
    fn from(point: &'a Point) -> ClusteredPoint<'a> {
        ClusteredPoint {
            point,
            cluster_id: 0,
        }
    }

    fn distance_to(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }

    fn distance_to_sqr(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }
}

struct Clusterer<'a> {
    points: Vec<ClusteredPoint<'a>>,

    leaf_size: usize,
}

impl<'a> Clusterer<'a> {
    fn new() -> Clusterer<'a> {
        Clusterer {
            points: vec![],
            leaf_size: 50,
        }
    }

    fn from_data(data: &'a Vec<Point>) -> Clusterer {
        Clusterer {
            points: data
                .iter()
                .map(|point| ClusteredPoint::from(point))
                .collect(),

            leaf_size: 50,
        }
    }

    fn set_leaf_size(&'a mut self, leaf_size: usize) {
        self.leaf_size = leaf_size;
    }

    fn count(&self) -> usize {
        self.points.len()
    }
}

//
//
// Ball Tree
//
//
//

struct BranchData<'a> {
    children: (Box<BallTree<'a>>, Box<BallTree<'a>>),
    radius: f64,
    pivot: Point,
}
struct LeafData<'a> {
    member_data: Vec<ClusteredPoint<'a>>,
    radius: f64,
    pivot: Point,
}

enum BallTree<'a> {
    Branch(BranchData<'a>),
    Leaf(LeafData<'a>),
}

impl<'a> BallTree<'a> {
    fn size(&self) -> usize {
        match self {
            BallTree::Leaf(leaf) => leaf.member_data.len(),
            BallTree::Branch(branch) => branch.children.0.size() + branch.children.1.size(),
        }
    }

    fn radius(&self) -> f64 {
        match self {
            BallTree::Leaf(leaf) => leaf.radius,
            BallTree::Branch(branch) => branch.radius,
        }
    }

    fn create_ball_tree(data: &mut Vec<&'a Point>, leaf_size: usize) -> Box<BallTree<'a>> {
        if data.len() < min(leaf_size, 3) {
            let result = Box::new(BallTree::Leaf(LeafData {
                member_data: data.iter().map(|p| ClusteredPoint::from(&p)).collect(),
                pivot: Point::new(),
                radius: 0.,
            }));
            result
        } else {
            let random_point = data[0].clone();

            let mut point_1 = data[1];
            let mut max_distance = point_1.distance_to(&random_point);

            for point in data.iter() {
                let distance: f64 = point.distance_to(&random_point);
                if distance > max_distance {
                    max_distance = distance;
                    point_1 = point;
                }
            }

            let mut point_2 = data[0];
            let mut max_distance = point_1.distance_to(&random_point);

            for point in data.iter() {
                let distance: f64 = point.distance_to(&random_point);

                if distance > max_distance {
                    max_distance = distance;
                    point_2 = point;
                }
            }
            
            data.sort_unstable_by(|x, y| {
                let xr_1 = x.distance_to(&point_1);
                let xr_2 = x.distance_to(&point_2);

                let yr_1 = y.distance_to(&point_1);
                let yr_2 = y.distance_to(&point_2);

                let v_1 = (xr_1 - xr_2) / (xr_1 + xr_2);
                let v_2 = (yr_1 - yr_2) / (yr_1 + yr_2);

                v_2.partial_cmp(&v_2).unwrap_or(std::cmp::Ordering::Equal)
            });

            let half_index = data.len() / 2;

            let mut vec_0: Vec<&Point> = data[..half_index].iter().map(|&x| x.clone()).collect();
            let mut vec_1: Vec<&Point> = data[half_index..].iter().map(|&x| x.clone()).collect();

            let result = Box::new(BallTree::Branch(BranchData {
                children: (
                    BallTree::create_ball_tree(&mut vec_0, leaf_size),
                    BallTree::create_ball_tree(&mut vec_1, leaf_size),
                ),
                radius: 0.,
                pivot: Point::new(),
            }));

            return result;
        }
    }
}
