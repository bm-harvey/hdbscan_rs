use crate::ball_tree::BallTree;
use crate::point::Point;

use std::cell::RefCell;
use std::rc::Rc;

type PointRef = Rc<Point>;

pub struct ClusteredPoint {
    point: PointRef,
    //neighbors: Vec<&'a ClusteredPoint<'a>>,
    neighbors: Vec<Rc<RefCell<ClusteredPoint>>>,
    cluster_id: u16,
}

impl ClusteredPoint {
    pub fn from(point: PointRef) -> ClusteredPoint {
        ClusteredPoint {
            point,
            neighbors: vec![],
            cluster_id: 0,
        }
    }
    pub fn point(&self) -> &Point {
        &self.point
    }
    //pub fn neighbors(&'a self) -> &'a Vec<&'a ClusteredPoint<'a>> {
    //&self.neighbors
    //}

    pub fn neighbors(&self) -> &Vec<Rc<RefCell<ClusteredPoint>>> {
        &self.neighbors
    }

    pub fn set_neighbors(&mut self, neighbors: Vec<Rc<RefCell<ClusteredPoint>>>) -> &mut Self {
        self.neighbors = neighbors;
        self
    }

    pub fn set_neighbor_capacity(&mut self, capacity: usize) -> &mut Self {
        self.neighbors.reserve(capacity);
        self
    }

    pub fn clear_neighbors(&mut self) -> &mut Self {
        self.neighbors.clear();
        self
    }

    pub fn core_distance(&self) -> f64 {
        return self
            .point()
            .distance_to(self.neighbors.last().unwrap().borrow().point());
    }

    pub fn distance_to(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point.as_ref())
    }

    pub fn distance_to_sqr(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point.as_ref())
    }
}

pub struct Clusterer<'a> {
    // data
    point_cloud: &'a [Rc<Point>],
    // parameters
    leaf_size: usize,
    param_k: usize,
}

impl<'a> Clusterer<'a> {
    pub fn new(data: &'a [Rc<Point>]) -> Clusterer {
        Clusterer {
            point_cloud: data,
            leaf_size: 50,
            param_k: 15,
        }
    }

    pub fn with_leaf_size(mut self, leaf_size: usize) -> Clusterer<'a> {
        self.leaf_size = leaf_size;
        self
    }

    pub fn with_param_k(mut self, param_k: usize) -> Clusterer<'a> {
        self.param_k = param_k;
        self
    }

    pub fn fit(self) -> ClusterResult {
        let mut data_refs = self.point_cloud.iter().map(Rc::clone).collect();

        println!("Generating spatial indexing tree");
        let spatial_index_root = BallTree::new(&mut data_refs, self.leaf_size);

        spatial_index_root.set_k_neareset_neighbors(self.param_k);

        let result = ClusterResult { spatial_index_root };

        let _size = result.spatial_index_root.size();

        println!("Finding {}-nearest neighbors for all data", self.param_k);

        result
    }
}

pub struct ClusterResult {
    spatial_index_root: Box<BallTree>,
}

impl ClusterResult {
    pub fn builder(data: &[Rc<Point>]) -> Clusterer {
        Clusterer::new(data)
    }

    pub fn ball_tree(&self) -> &BallTree {
        &self.spatial_index_root
    }

    //pub fn ball_tree_iter(&self) -> BallTreeItr {
    //self.spatial_index_root.iter()
    //}
}
