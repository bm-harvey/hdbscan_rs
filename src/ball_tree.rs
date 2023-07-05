//use std::borrow::Borrow;
//use std::borrow::BorrowMut;
//use std::cell::Ref;
use std::cell::RefCell;
//use std::cell::RefMut;
use std::cmp::max;
use std::rc::Rc;

use crate::clusterer::ClusteredPoint;
use crate::point::Point;

pub struct BranchData<'a> {
    children: (Box<BallTree<'a>>, Box<BallTree<'a>>),
    radius: f64,
    pivot: Point,
}

type Rcc<T> = Rc<RefCell<T>>;

pub struct LeafData<'a> {
    member_data: Vec<RefCell<ClusteredPoint<'a>>>,
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

    pub fn set_k_neareset_neigbors(&'a mut self, param_k: usize) {
        for point_ref in self.iter() {
            point_ref
                .borrow_mut()
                .clear_neighbors()
                .set_neighbor_capacity(param_k);
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
        let pivot = BallTree::pivot_from_data(data);
        let radius = BallTree::radius_from_data(data, &pivot);

        if data.len() < max(leaf_size, 3) {
            Box::new(BallTree::Leaf(LeafData {
                member_data: data
                    .iter()
                    .map(|point| RefCell::new(ClusteredPoint::from(point)))
                    .collect(),
                pivot,
                radius,
            }))
        } else {
            let random_point = data[0];

            let point_1: &Point = data
                .iter()
                .max_by(|x, y| {
                    x.distance_to_sqr(random_point)
                        .total_cmp(&y.distance_to_sqr(random_point))
                })
                .unwrap();

            let point_2: &Point = data
                .iter()
                .max_by(|x, y| {
                    x.distance_to_sqr(point_1)
                        .total_cmp(&y.distance_to_sqr(point_1))
                })
                .unwrap();

            data.sort_unstable_by(|x, y| {
                let xr_1 = x.distance_to_sqr(point_1);
                let xr_2 = x.distance_to_sqr(point_2);

                let yr_1 = y.distance_to_sqr(point_1);
                let yr_2 = y.distance_to_sqr(point_2);

                let v_1 = (xr_1 - xr_2) / (xr_1 + xr_2);
                let v_2 = (yr_1 - yr_2) / (yr_1 + yr_2);

                v_1.partial_cmp(&v_2).unwrap_or(std::cmp::Ordering::Equal)
            });

            let half_index = data.len() / 2;

            let mut vec_0: Vec<&Point> = data[..half_index].iter().map(<&Point>::clone).collect();
            let mut vec_1: Vec<&Point> = data[half_index..].iter().map(<&Point>::clone).collect();

            Box::new(BallTree::Branch(BranchData {
                children: (
                    BallTree::new(&mut vec_0, leaf_size),
                    BallTree::new(&mut vec_1, leaf_size),
                ),
                radius,
                pivot,
            }))
        } // else
    } // new

    fn pivot_from_data(data: &[&Point]) -> Point {
        let dims = data[0].num_dimensions();
        let mut pivot = Point::from(vec![0.; dims]);

        for idx in 0..dims {
            pivot.set(
                idx,
                data.iter().map(|point| point.get(idx)).sum::<f64>() / (data.len() as f64),
            );
        }
        pivot
    }

    fn radius_from_data(data: &[&Point], pivot: &Point) -> f64 {
        data.iter()
            .map(|x| x.distance_to(pivot))
            .max_by(|x, y| x.total_cmp(y))
            .unwrap_or(std::f64::MAX)
    }

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
    type Item = &'a RefCell<ClusteredPoint<'a>>;

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
                let point_cloud = &itr.data.member_data;
                if index < point_cloud.len() {
                    itr.counter += 1;
                    Some(&point_cloud[index])
                } else {
                    None
                }
            }
        } // match
    } // fn next
} // impl Iterator for BallTreeIter
