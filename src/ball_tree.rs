use std::cmp::min;

use crate::clusterer::ClusteredPoint;
use crate::point::Point;

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
    pub fn new(data: &mut Vec<&'a Point>, leaf_size: usize) -> Box<BallTree<'a>> {
        if data.len() < min(leaf_size, 3) {
            Box::new(BallTree::Leaf(LeafData {
                member_data: data
                    .iter()
                    .map(|point| ClusteredPoint::from(point))
                    .collect(),
                pivot: Point::new(),
                radius: 0.,
            }))
        } else {
            let random_point = data[0];

            let mut point_1 = data[1];
            let mut max_distance = point_1.distance_to(random_point);

            for point in data.iter() {
                let distance: f64 = point.distance_to(random_point);
                if distance > max_distance {
                    max_distance = distance;
                    point_1 = point;
                }
            }

            let mut point_2 = data[0];
            let mut max_distance = point_1.distance_to(random_point);

            for point in data.iter() {
                let distance: f64 = point.distance_to(random_point);

                if distance > max_distance {
                    max_distance = distance;
                    point_2 = point;
                }
            }

            data.sort_unstable_by(|x, y| {
                let xr_1 = x.distance_to(point_1);
                let xr_2 = x.distance_to(point_2);

                let yr_1 = y.distance_to(point_1);
                let yr_2 = y.distance_to(point_2);

                let v_1 = (xr_1 - xr_2) / (xr_1 + xr_2);
                let v_2 = (yr_1 - yr_2) / (yr_1 + yr_2);

                v_2.partial_cmp(&v_2).unwrap_or(std::cmp::Ordering::Equal)
            });

            let half_index = data.len() / 2;

            let mut vec_0: Vec<&Point> = data[..half_index].iter().map(<&Point>::clone).collect();
            let mut vec_1: Vec<&Point> = data[half_index..].iter().map(<&Point>::clone).collect();

            Box::new(BallTree::Branch(BranchData {
                children: (
                    BallTree::new(&mut vec_0, leaf_size),
                    BallTree::new(&mut vec_1, leaf_size),
                ),
                radius: 0.,
                pivot: Point::new(),
            }))
        } // else
    } // new
    
    pub fn iter(&'a self) -> BallTreeItr<'a> {
        match self {
            BallTree::Branch(tree) => BallTreeItr::Branch(BranchItrData {
                data: tree,
                child_is_left: true,
                child_itr: Box::new(tree.children.0.iter()),
            }),
            BallTree::Leaf(tree) => BallTreeItr::Leaf(LeafItrData {
                data: tree,
                counter: 0,
            }),
        }
    }
} // impl BallTree

pub struct BranchItrData<'a> {
    data: &'a BranchData<'a>,
    child_is_left: bool,
    child_itr: Box<BallTreeItr<'a>>,
}
pub struct LeafItrData<'a> {
    data: &'a LeafData<'a>,
    counter: usize,
}
pub enum BallTreeItr<'a> {
    Branch(BranchItrData<'a>),
    Leaf(LeafItrData<'a>),
}

impl<'a> Iterator for BallTreeItr<'a> {
    type Item = &'a ClusteredPoint<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            BallTreeItr::Branch(itr) => {
                let result = itr.child_itr.next();
                match result {
                    None => {
                        if itr.child_is_left {
                            itr.child_is_left = false;
                            itr.child_itr = Box::new(itr.data.children.1.iter());
                            itr.child_itr.next()
                        } else {
                            None
                        }
                    }
                    Some(res) => Some(res),
                }
            }
            BallTreeItr::Leaf(itr) => {
                let index = itr.counter;
                let pc = &itr.data.member_data;
                if index < pc.len() {
                    itr.counter += 1;
                    Some(&pc[index])
                } else {
                    None
                }
            }
        } // match
    } // fn next
} // impl Iterator for BallTreeIter
