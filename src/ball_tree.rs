use std::cell::RefCell;
use std::cmp::max;
use std::rc::Rc;

use std::cmp::min;

use crate::clusterer::ClusteredPoint;
use crate::point::Point;

type ClusterPntRef = Rc<RefCell<ClusteredPoint>>;

#[derive(Copy, Clone)]
pub enum Metric {
    Euclidean,
}
pub enum MutualReachability {
    Yes(Rc<BallTree>, usize),
    No,
}

pub struct BallTreeBuilder<'a> {
    data: &'a mut Vec<ClusterPntRef>,
    metric: Metric,
    mutual_reachability: MutualReachability,
    leaf_size: usize,
}

impl<'a> BallTreeBuilder<'a> {
    pub fn new(data: &'a mut Vec<ClusterPntRef>) -> Self {
        Self {
            data,
            metric: Metric::Euclidean,
            mutual_reachability: MutualReachability::No,
            leaf_size: 15,
        }
    }
    pub fn with_metric(mut self, metric: Metric) -> Self {
        self.metric = metric;
        self
    }

    pub fn with_mutual_reachability(mut self, mutual_reachability: MutualReachability) -> Self {
        self.mutual_reachability = mutual_reachability;
        self
    }

    pub fn build(mut self) -> Rc<BallTree> {
        BallTree::new(
            &mut self.data,
            self.metric,
            Rc::new(self.mutual_reachability),
            self.leaf_size,
        )
    }
}

pub struct BranchData {
    children: (Rc<BallTree>, Rc<BallTree>),
    radius: f64,
    pivot: ClusterPntRef,
    metric: Metric,
    mutual_reachability: Rc<MutualReachability>,
}

pub struct LeafData {
    member_data: Vec<ClusterPntRef>,
    radius: f64,
    pivot: ClusterPntRef,
    metric: Metric,
    mutual_reachability: Rc<MutualReachability>,
}

pub enum BallTree {
    Branch(Rc<BranchData>),
    Leaf(Rc<LeafData>),
}

impl BallTree {
    /// Counts all of the `ClusterPoints` in the tree starting from the `self` node
    pub fn size(&self) -> usize {
        match self {
            BallTree::Leaf(leaf) => leaf.member_data.len(),
            BallTree::Branch(branch) => branch.children.0.size() + branch.children.1.size(),
        }
    }

    fn metric(&self) -> Metric {
        match self {
            BallTree::Leaf(leaf_data) => leaf_data.metric,
            BallTree::Branch(branch_data) => branch_data.metric,
        }
    }

    fn mutual_reachability(&self) -> &MutualReachability {
        match self {
            BallTree::Leaf(leaf_data) => &leaf_data.mutual_reachability,
            BallTree::Branch(branch_data) => &branch_data.mutual_reachability,
        }
    }

    pub fn core_distance(point: &ClusterPntRef, metric: Metric) -> f64 {
        let point = point.borrow();

        BallTree::distance_between_points(
            point.furthest_neighbor().unwrap().borrow().point(),
            point.point(),
            metric,
        )
    }

    pub fn distance_between_cluster_points(
        point_1: &ClusterPntRef,
        point_2: &ClusterPntRef,
        metric: Metric,
        mutual_reachability: &MutualReachability,
    ) -> f64 {
        let distance = BallTree::distance_between_points(
            &point_1.borrow().point(),
            &point_2.borrow().point(),
            metric,
        );

        match mutual_reachability {
            MutualReachability::No => distance,
            MutualReachability::Yes(ball_tree, _) => {
                let point_1 = point_1.borrow();
                let point_2 = point_2.borrow();

                let core_1 = BallTree::distance_between_points(
                    point_1.furthest_neighbor().unwrap().borrow().point(),
                    point_1.point(),
                    ball_tree.metric(),
                );

                let core_2 = BallTree::distance_between_points(
                    point_2.furthest_neighbor().unwrap().borrow().point(),
                    point_2.point(),
                    ball_tree.metric(),
                );

                let max_core = f64::max(core_1, core_2);
                let mutual_reachability_distance = f64::max(max_core, distance);
                mutual_reachability_distance
            }
        }
    } // pub fn distance_betwee_cluster_points

    pub fn distance_between_points(point_1: &Point, point_2: &Point, metric: Metric) -> f64 {
        let dims = min(point_1.num_dimensions(), point_2.num_dimensions());

        match metric {
            Metric::Euclidean => (0..dims)
                .map(|idx| point_1.get(idx) - point_2.get(idx))
                .map(|diff| diff.powi(2))
                .sum::<f64>()
                .sqrt(),
        }
    }

    /// Iterates over the contained data and sets the neighbors
    /// This function is called from a clusterer and passes the heavy lifting to the
    /// `find_nearest_neighbors` function
    pub fn set_k_neareset_neighbors(&self, param_k: usize) {
        for point_ref in self.cluster_points() {
            let mut neighbors: Vec<ClusterPntRef> = Vec::new();
            neighbors.reserve(param_k);
            self.find_k_nearest_neighbors(&mut neighbors, &point_ref, param_k);
            point_ref.borrow_mut().set_neighbors(neighbors);
        }
    }

    /// Iterates over the pivots and sets the neighbors
    /// This function is called from a clusterer and passes the heavy lifting to the
    /// `find_nearest_neighbors` function
    pub fn set_k_neareset_neighbors_for_pivots(&self, param_k: usize) {
        for point_ref in self.pivots() {
            let mut neighbors: Vec<ClusterPntRef> = Vec::new();
            neighbors.reserve(param_k);
            self.find_k_nearest_neighbors(&mut neighbors, &point_ref, param_k);
            point_ref.borrow_mut().set_neighbors(neighbors);
        }
    }

    /// The result of the search is stored in `neighbors`. The top level call of this function
    /// expects that `neighbors` is empty, and that the top level call originated from the root
    /// node. The algorithm is efficient for param_k << self.size(), with O( k**2 log (n)).
    fn find_k_nearest_neighbors(
        &self,
        neighbors: &mut Vec<ClusterPntRef>,
        target: &ClusterPntRef,
        param_k: usize,
    ) {
        let at_capacity: bool = param_k == neighbors.len();

        let target_ref = target.borrow();
        let target_pnt = target_ref.point();

        // the points list already has k points, and nothing in this ball is closer than the
        // already found kth nearest neighbor, so skip this node entirely.
        if at_capacity {
            let distance_to_ball = BallTree::distance_between_cluster_points(
                &self.pivot(),
                target,
                self.metric(),
                &self.mutual_reachability(),
            ) - self.radius();

            let distance_to_kth_nearest = BallTree::distance_between_cluster_points(
                &target,
                neighbors.last().unwrap(),
                self.metric(),
                &self.mutual_reachability(),
            );

            let ball_is_too_far = distance_to_ball > distance_to_kth_nearest;

            if ball_is_too_far {
                return;
            }
        }

        // there are points in this ball that are closer than some of the points in `neighbors`
        match self {
            // branches always defer the heavy lifting to the children
            BallTree::Branch(branch_data) => {
                // It is more likely that the nearer neighbors are in the nearer balls.
                // This check makes more of the "at_capacity" skips happen.
                let distance_to_0 = BallTree::distance_between_cluster_points(
                    &branch_data.children.0.pivot(),
                    target,
                    self.metric(),
                    &self.mutual_reachability(),
                );
                let distance_to_1 = BallTree::distance_between_cluster_points(
                    &branch_data.children.1.pivot(),
                    target,
                    self.metric(),
                    &self.mutual_reachability(),
                );

                if distance_to_0 < distance_to_1 {
                    branch_data
                        .children
                        .0
                        .find_k_nearest_neighbors(neighbors, target, param_k);
                    branch_data
                        .children
                        .1
                        .find_k_nearest_neighbors(neighbors, target, param_k);
                } else {
                    branch_data
                        .children
                        .1
                        .find_k_nearest_neighbors(neighbors, target, param_k);
                    branch_data
                        .children
                        .0
                        .find_k_nearest_neighbors(neighbors, target, param_k);
                }
            }
            BallTree::Leaf(leaf_data) => {
                // there is nothing better to do here than brute force through the vector of points
                // inserts keep the vector ordered, so the last point is always the furthest away
                for point in leaf_data.member_data.iter() {
                    // don't count the point as it's own neighbor
                    if std::ptr::eq(point.borrow().point(), target_pnt) {
                        continue;
                    }

                    let distance_to_target = BallTree::distance_between_cluster_points(
                        &point,
                        &target,
                        self.metric(),
                        &self.mutual_reachability(),
                    );

                    // pos is where we are about to insert `point` into the neighbors. None means
                    // that there are already k neighbors, and `point` isn't closer than any of
                    // them.
                    let pos = if (neighbors.len() < param_k)
                        || (distance_to_target
                            < BallTree::distance_between_cluster_points(
                                &neighbors.last().unwrap(),
                                &target,
                                self.metric(),
                                &self.mutual_reachability(),
                            )) {
                        Some(neighbors.partition_point(|x| {
                            BallTree::distance_between_cluster_points(
                                x,
                                &target,
                                self.metric(),
                                &self.mutual_reachability(),
                            ) < distance_to_target
                        }))
                    } else {
                        None
                    };
                    if let Some(idx) = pos {
                        // pop before insert to avoid copying the vector for more allocation
                        if neighbors.len() >= param_k {
                            neighbors.pop();
                        }
                        neighbors.insert(idx, Rc::clone(point));
                    }
                }
            }
        }
    }

    /// returns the stored value of radius
    fn radius(&self) -> f64 {
        match self {
            BallTree::Leaf(leaf) => leaf.radius,
            BallTree::Branch(branch) => branch.radius,
        }
    }

    /// returns reference to the stored pivot
    fn pivot(&self) -> &ClusterPntRef {
        match self {
            BallTree::Leaf(leaf) => &leaf.pivot,
            BallTree::Branch(branch) => &branch.pivot,
        }
    }

    /// returns a pivot calculated from the data
    fn pivot_from_data(data: &[ClusterPntRef]) -> ClusterPntRef {
        let dims = data[0].borrow().point().num_dimensions();
        let mut pivot_coord = vec![0.; dims];

        for idx in 0..dims {
            pivot_coord[idx] = data
                .iter()
                .map(|point| point.borrow().point().get(idx))
                .sum::<f64>()
                / (data.len() as f64);
        }

        ClusteredPoint::from_as_rcc(pivot_coord)
    }

    /// returns a radious calculated from the data
    fn radius_from_data(
        data: &[ClusterPntRef],
        pivot: &ClusterPntRef,
        metric: Metric,
        mutual_reachability: &Rc<MutualReachability>,
    ) -> f64 {
        data.iter()
            .map(|x| {
                BallTree::distance_between_cluster_points(x, &pivot, metric, &mutual_reachability)
            })
            .max_by(|x, y| x.total_cmp(y))
            .unwrap()
    }

    /// recursive new. When called externally from this function will return the root node of the
    /// ball tree (boxed). No mutation of the `Point`s actually occurs here, but the data is
    /// sorted in each step, so the `Vec` must be a mutable ref. If the length of data is larger
    /// than leaf_size, a branch is returned containing 2 children ball trees, each with half of
    /// the data.  
    pub fn new(
        data: &mut Vec<ClusterPntRef>,
        metric: Metric,
        mutual_reachability: Rc<MutualReachability>,
        leaf_size: usize,
    ) -> Rc<BallTree> {
        let pivot = BallTree::pivot_from_data(data);

        match &*mutual_reachability {
            MutualReachability::Yes(ball_tree, param_k) => {
                let mut neighbors: Vec<ClusterPntRef> = Vec::new();
                neighbors.reserve(*param_k);
                ball_tree.find_k_nearest_neighbors(&mut neighbors, &pivot, *param_k);
                pivot.borrow_mut().set_neighbors(neighbors);
            }
            MutualReachability::No => {}
        }

        let radius = BallTree::radius_from_data(data, &pivot, metric, &mutual_reachability);

        if data.len() < max(leaf_size, 3) {
            Rc::new(BallTree::Leaf(Rc::new(LeafData {
                member_data: data.iter().map(Rc::clone).collect(),
                pivot,
                radius,
                metric,
                mutual_reachability,
            })))
        } else {
            let random_point = Rc::clone(&data[0]);

            let point_1 = Rc::clone(
                data.iter()
                    .max_by(|x, y| {
                        let x_dist = BallTree::distance_between_cluster_points(
                            x,
                            &random_point,
                            metric,
                            &mutual_reachability,
                        );
                        let y_dist = BallTree::distance_between_cluster_points(
                            y,
                            &random_point,
                            metric,
                            &mutual_reachability,
                        );
                        x_dist.total_cmp(&y_dist)
                    })
                    .unwrap(),
            );

            data.sort_unstable_by(|x, y| {
                let x_dist = BallTree::distance_between_cluster_points(
                    x,
                    &point_1,
                    metric,
                    &mutual_reachability,
                );
                let y_dist = BallTree::distance_between_cluster_points(
                    y,
                    &point_1,
                    metric,
                    &mutual_reachability,
                );

                x_dist
                    .partial_cmp(&y_dist)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            let half_index = data.len() / 2;

            let mut vec_0 = data[..half_index].iter().map(Rc::clone).collect();
            let mut vec_1 = data[half_index..].iter().map(Rc::clone).collect();

            Rc::new(BallTree::Branch(Rc::new(BranchData {
                children: (
                    BallTree::new(
                        &mut vec_0,
                        metric,
                        Rc::clone(&mutual_reachability),
                        leaf_size,
                    ),
                    BallTree::new(
                        &mut vec_1,
                        metric,
                        Rc::clone(&mutual_reachability),
                        leaf_size,
                    ),
                ),
                radius,
                pivot,
                metric,
                mutual_reachability,
            })))
        } // else
    } // new

    /// Returns an iterator over the stored data via a recursive depth first search.
    pub fn cluster_points(&self) -> BallTreeDataItr {
        match self {
            BallTree::Branch(tree) => BallTreeDataItr::Branch(BranchItrData {
                data: Rc::clone(tree),
                child_is_left: true,
                child_itr: Box::new(tree.children.0.cluster_points()),
            }),
            BallTree::Leaf(tree) => BallTreeDataItr::Leaf(LeafItrData {
                data: Rc::clone(tree),
                counter: 0,
            }),
        }
    }

    /// iterate over the pivot points of the ball trees
    pub fn pivots(&self) -> BallTreePivotItr {
        match self {
            BallTree::Branch(tree) => BallTreePivotItr::Branch(BranchPivotItrData {
                data: Rc::clone(tree),
                status: PivotItrStatus::Left,
                child_itr: Box::new(tree.children.0.pivots()),
            }),
            BallTree::Leaf(tree) => BallTreePivotItr::Leaf(LeafPivotItrData {
                data: Rc::clone(tree),
                previously_returned: false,
            }),
        }
    }
} // impl BallTree

enum PivotItrStatus {
    Left,
    Right,
    Here,
    Done,
}

pub struct BranchPivotItrData {
    data: Rc<BranchData>,
    status: PivotItrStatus,
    child_itr: Box<BallTreePivotItr>,
}

pub struct LeafPivotItrData {
    data: Rc<LeafData>,
    previously_returned: bool,
}

pub enum BallTreePivotItr {
    Branch(BranchPivotItrData),
    Leaf(LeafPivotItrData),
}

impl Iterator for BallTreePivotItr {
    type Item = ClusterPntRef;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            BallTreePivotItr::Branch(data) => match data.status {
                PivotItrStatus::Done => None,
                PivotItrStatus::Here => {
                    data.status = PivotItrStatus::Done;
                    Some(Rc::clone(&data.data.pivot))
                }
                PivotItrStatus::Left => {
                    let result = data.child_itr.next();
                    match result {
                        Some(itr) => Some(itr),
                        None => {
                            data.status = PivotItrStatus::Right;
                            data.child_itr = Box::new(data.data.children.1.pivots());
                            data.child_itr.next()
                        }
                    }
                }
                PivotItrStatus::Right => {
                    let result = data.child_itr.next();
                    match result {
                        Some(itr) => Some(itr),
                        None => {
                            data.status = PivotItrStatus::Here;
                            data.child_itr = Box::new(data.data.children.1.pivots());
                            data.child_itr.next()
                        }
                    }
                }
            },
            BallTreePivotItr::Leaf(data) => {
                if data.previously_returned {
                    None
                } else {
                    data.previously_returned = true;
                    Some(Rc::clone(&data.data.pivot))
                }
            }
        } // match
    } // fn next
}

pub struct BranchItrData {
    data: Rc<BranchData>,
    child_is_left: bool,
    child_itr: Box<BallTreeDataItr>,
}

pub struct LeafItrData {
    data: Rc<LeafData>,
    counter: usize,
}

pub enum BallTreeDataItr {
    Branch(BranchItrData),
    Leaf(LeafItrData),
}
/// Iterator for the stored `ClusterPoint`s  
impl Iterator for BallTreeDataItr {
    type Item = Rc<RefCell<ClusteredPoint>>;

    /// For branches, deffer to child 0 until that returns None, then deffer to child 1. When child
    /// 1 returns None, return None. For leaves, keep track of the index and return a clone of the
    /// appropriate `ClusterPoint`
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            BallTreeDataItr::Branch(itr) => {
                let result = itr.child_itr.next();
                match result {
                    None => {
                        if itr.child_is_left {
                            itr.child_is_left = false;
                            itr.child_itr = Box::new(itr.data.children.1.cluster_points());
                            itr.child_itr.next()
                        } else {
                            None
                        }
                    }
                    Some(res) => Some(res),
                }
            }
            BallTreeDataItr::Leaf(itr) => {
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

#[cfg(test)]
mod tests {

    use crate::ball_tree::BallTree;
    use crate::ClusteredPoint;
    use rand::prelude::*;
    use std::rc::Rc;

    #[test]
    fn nearest_neighbors_correct() {
        let num_points: usize = 1_000;
        let param_k = 15;
        let leaf_size = 15;

        let mut rng = rand::thread_rng();

        let mut data = vec![];
        data.reserve(num_points);
        for _ in 0..num_points {
            data.push(ClusteredPoint::from_as_rcc(vec![
                rng.gen::<f64>(),
                rng.gen::<f64>(),
            ]))
        }

        let ball_tree = BallTree::new(
            &mut data,
            crate::Metric::Euclidean,
            Rc::new(crate::MutualReachability::No),
            leaf_size,
        );
        ball_tree.set_k_neareset_neighbors(param_k);

        for target_point_ref in ball_tree.cluster_points() {
            let core_distance =
                BallTree::core_distance(&target_point_ref, crate::Metric::Euclidean);

            let mut num_points_within_core = 0;

            for other_point_ref in ball_tree.cluster_points() {
                let distance = BallTree::distance_between_cluster_points(
                    &other_point_ref,
                    &target_point_ref,
                    crate::Metric::Euclidean,
                    &Rc::new(crate::MutualReachability::No),
                );

                if distance <= core_distance {
                    num_points_within_core += 1;
                }
            }
            assert_eq!(param_k + 1, num_points_within_core);
        }
    }
}
