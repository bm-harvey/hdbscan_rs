use std::cmp::min;

use crate::point::Point;
use crate::clusterer::ClusteredPoint;

pub struct BranchData<'a> {
    children: (Box<BallTree<'a>>, Box<BallTree<'a>>),
    radius: f64,
    pivot: Point,
}
pub struct LeafData<'a> {
    member_data: Vec<ClusteredPoint<'a>>,
    radius: f64,
    pivot: Point,
}

pub enum BallTree<'a> {
    Branch(BranchData<'a>),
    Leaf(LeafData<'a>),
}

impl<'a> BallTree<'a> {
    
    // access fns
    pub fn size(&self) -> usize {
        match self {
            BallTree::Leaf(leaf) => leaf.member_data.len(),
            BallTree::Branch(branch) => branch.children.0.size() + branch.children.1.size(),
        }
    }

    pub fn radius(&self) -> f64 {
        match self {
            BallTree::Leaf(leaf) => leaf.radius,
            BallTree::Branch(branch) => branch.radius,
        }
    }

    // ctor fns
    pub fn create_ball_tree(data: &mut Vec<&'a Point>, leaf_size: usize) -> Box<BallTree<'a>> {
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
