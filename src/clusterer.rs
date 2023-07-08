use crate::ball_tree::BallTree;
use crate::ball_tree::Metric;
use crate::ball_tree::MutualReachability;

use crate::point::Point;

use std::cell::RefCell;
use std::rc::Rc;

use std::time::Instant;

type PointRef = Rc<Point>;

pub struct ClusteredPoint {
    point: PointRef,
    neighbors: Vec<Rc<RefCell<ClusteredPoint>>>,
}

impl ClusteredPoint {
    pub fn from(point: PointRef) -> ClusteredPoint {
        ClusteredPoint {
            point,
            neighbors: vec![],
        }
    }
    pub fn point(&self) -> &Point {
        &self.point
    }

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

pub struct Clusterer {
    // data
    point_cloud: Vec<Rc<Point>>,
    // parameters
    leaf_size: usize,
    param_k: usize,
}

impl Clusterer {
    pub fn new(data: Vec<Rc<Point>>) -> Clusterer {
        Clusterer {
            point_cloud: data,
            leaf_size: 50,
            param_k: 15,
        }
    }

    pub fn with_leaf_size(mut self, leaf_size: usize) -> Clusterer {
        self.leaf_size = leaf_size;
        self
    }

    pub fn with_param_k(mut self, param_k: usize) -> Clusterer {
        self.param_k = param_k;
        self
    }

    pub fn fit(self) -> ClusterResult {
        let mut data_refs = self.point_cloud.iter().map(Rc::clone).collect();

        let algorithm_timer_start = Instant::now();

        println!("Generating spatial indexing tree ...");
        let start = Instant::now();
        let spatial_index_root =
            BallTree::new(&mut data_refs, Metric::Euclidean, MutualReachability::No, self.leaf_size);
        println!("\tdone in {} s", start.elapsed().as_secs_f32());

        println!("Finding {}-nearest neighbors for all data...", self.param_k);
        let start = Instant::now();
        spatial_index_root.set_k_neareset_neighbors(self.param_k);
        println!("\tdone in {} s", start.elapsed().as_secs_f32());

        let result = ClusterResult { spatial_index_root };

        println!(
            "The algorithm took {} s all together",
            algorithm_timer_start.elapsed().as_secs_f32(),
        );

        result
    }
}

pub struct ClusterResult {
    spatial_index_root: Box<BallTree>,
}

impl ClusterResult {
    pub fn ball_tree(&self) -> &BallTree {
        &self.spatial_index_root
    }
}
