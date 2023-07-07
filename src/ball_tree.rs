use std::cell::RefCell;
use std::cmp::max;
use std::rc::Rc;

use crate::clusterer::ClusteredPoint;
use crate::point::Point;

pub struct BranchData {
    children: (Box<BallTree>, Box<BallTree>),
    radius: f64,
    pivot: Point,
}

pub struct LeafData {
    member_data: Vec<Rc<RefCell<ClusteredPoint>>>,
    radius: f64,
    pivot: Point,
}

pub enum BallTree {
    Branch(Rc<BranchData>),
    Leaf(Rc<LeafData>),
}

impl BallTree {
    // access fns
    pub fn size(&self) -> usize {
        match self {
            BallTree::Leaf(leaf) => leaf.member_data.len(),
            BallTree::Branch(branch) => branch.children.0.size() + branch.children.1.size(),
        }
    }

    pub fn find_k_nearest_neighbors(
        &self,
        //neighbors: &mut Vec<&'a ClusteredPoint<'a>>,
        neighbors: &mut Vec<Rc<RefCell<ClusteredPoint>>>,
        target_point: &Point,
        param_k: usize,
    ) {
        let at_capacity: bool = param_k == neighbors.len();

        // the points list already has k points, and nothing in this ball is close enough
        // to the target point to benifit, so we should return before going further with this node
        // this is the step that saves most of time and makes the query O(klog(n))
        if at_capacity
            && self.pivot().distance_to(target_point) - self.radius()
                > target_point.distance_to(neighbors.last().unwrap().borrow().point())
        {
            return;
        }

        // there are points in this ball that are closer than some of the points in `neighbors`
        match self {
            // branches always defer the heavy lifting to the children
            BallTree::Branch(branch_data) => {
                let child_0_closer = branch_data.children.0.pivot().distance_to(target_point)
                    < branch_data.children.1.pivot().distance_to(target_point);

                if child_0_closer {
                    branch_data.children.0.find_k_nearest_neighbors(
                        neighbors,
                        target_point,
                        param_k,
                    );
                    branch_data.children.1.find_k_nearest_neighbors(
                        neighbors,
                        target_point,
                        param_k,
                    );
                } else {
                    branch_data.children.1.find_k_nearest_neighbors(
                        neighbors,
                        target_point,
                        param_k,
                    );
                    branch_data.children.0.find_k_nearest_neighbors(
                        neighbors,
                        target_point,
                        param_k,
                    );
                }
            }
            // loop through the points in this leaf and insert them sorted into neighbors keeping
            // the vector contained to param_k in length
            BallTree::Leaf(leaf_data) => {
                // nothing better to do here than brute force through the vector of points

                for point in leaf_data.member_data.iter() {
                    if std::ptr::eq(point.borrow().point(), target_point) {
                        continue;
                    }

                    let distance_to_target = point.borrow().point().distance_to(target_point);

                    let pos = if (neighbors.len() < param_k)
                        || (distance_to_target
                            < target_point.distance_to(neighbors.last().unwrap().borrow().point()))
                    {
                        Some(neighbors.partition_point(|x| {
                            x.borrow().point().distance_to(target_point) < distance_to_target
                        }))
                    } else {
                        None
                    };
                    if let Some(idx) = pos {
                        if neighbors.len() >= param_k {
                            neighbors.pop();
                        }
                        neighbors.insert(idx, Rc::clone(point));
                    }
                }
            }
        }
    }

    pub fn set_k_neareset_neighbors(&self, param_k: usize) {
        for point_ref in self.iter() {
            let mut neighbors: Vec<Rc<RefCell<ClusteredPoint>>> = Vec::new();
            self.find_k_nearest_neighbors(&mut neighbors, point_ref.borrow().point(), param_k);
            point_ref.borrow_mut().set_neighbors(neighbors);
        }
    }

    pub fn radius(&self) -> f64 {
        match self {
            BallTree::Leaf(leaf) => leaf.radius,
            BallTree::Branch(branch) => branch.radius,
        }
    }

    pub fn pivot(&self) -> &Point {
        match self {
            BallTree::Leaf(leaf) => &leaf.pivot,
            BallTree::Branch(branch) => &branch.pivot,
        }
    }

    // ctor fns
    pub fn new(data: &mut Vec<Rc<Point>>, leaf_size: usize) -> Box<BallTree> {
        let pivot = BallTree::pivot_from_data(data);
        let radius = BallTree::radius_from_data(data, &pivot);

        if data.len() < max(leaf_size, 3) {
            Box::new(BallTree::Leaf(Rc::new(LeafData {
                member_data: data
                    .iter()
                    .map(|point| Rc::new(RefCell::new(ClusteredPoint::from(Rc::clone(point)))))
                    .collect(),
                pivot,
                radius,
            })))
        } else {
            let random_point = Rc::clone(&data[0]);

            let point_1 = Rc::clone(
                data.iter()
                    .max_by(|x, y| {
                        x.distance_to_sqr(&random_point)
                            .total_cmp(&y.distance_to_sqr(&random_point))
                    })
                    .unwrap(),
            );

            //let point_2 = Rc::clone(
            //data.iter()
            //.max_by(|x, y| {
            //x.distance_to_sqr(&point_1)
            //.total_cmp(&y.distance_to_sqr(&point_1))
            //})
            //.unwrap(),
            //);

            data.sort_unstable_by(|x, y| {
                let xr_1 = x.distance_to_sqr(&point_1);
                //let xr_2 = x.distance_to(&point_2);

                let yr_1 = y.distance_to_sqr(&point_1);
                //let yr_2 = y.distance_to(&point_2);

                //let relative_diff_1 = (xr_1 - xr_2) / (xr_1 + xr_2);
                //let relative_diff_2 = (yr_1 - yr_2) / (yr_1 + yr_2);

                xr_1.partial_cmp(&yr_1).unwrap_or(std::cmp::Ordering::Equal)
            });

            let half_index = data.len() / 2;

            let mut vec_0: Vec<Rc<Point>> = data[..half_index].iter().map(Rc::clone).collect();
            let mut vec_1: Vec<Rc<Point>> = data[half_index..].iter().map(Rc::clone).collect();

            Box::new(BallTree::Branch(Rc::new(BranchData {
                children: (
                    BallTree::new(&mut vec_0, leaf_size),
                    BallTree::new(&mut vec_1, leaf_size),
                ),
                radius,
                pivot,
            })))
        } // else
    } // new

    fn pivot_from_data(data: &[Rc<Point>]) -> Point {
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

    fn radius_from_data(data: &[Rc<Point>], pivot: &Point) -> f64 {
        data.iter()
            .map(|x| x.distance_to(pivot))
            .max_by(|x, y| x.total_cmp(y))
            .unwrap()
    }

    pub fn iter(&self) -> BallTreeItr {
        match self {
            BallTree::Branch(tree) => BallTreeItr::Branch(BranchItrData {
                data: Rc::clone(tree),
                child_is_left: true,
                child_itr: Box::new(tree.children.0.iter()),
            }),
            BallTree::Leaf(tree) => BallTreeItr::Leaf(LeafItrData {
                data: Rc::clone(tree),
                counter: 0,
            }),
        }
    }
} // impl BallTree

pub struct BranchItrData {
    data: Rc<BranchData>,
    child_is_left: bool,
    child_itr: Box<BallTreeItr>,
}
pub struct LeafItrData {
    data: Rc<LeafData>,
    counter: usize,
}
pub enum BallTreeItr {
    Branch(BranchItrData),
    Leaf(LeafItrData),
}

impl Iterator for BallTreeItr {
    type Item = Rc<RefCell<ClusteredPoint>>;

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
                    Some(Rc::clone(&point_cloud[index]))
                } else {
                    None
                }
            }
        } // match
    } // fn next
} // impl Iterator for BallTreeIter
