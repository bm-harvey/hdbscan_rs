use crate::ball_tree::BallTree;
use crate::point::Point;

pub struct ClusteredPoint<'a> {
    point: &'a Point,
    neighbors: Vec<&'a ClusteredPoint<'a>>,
    cluster_id: u16,
}

impl<'a> ClusteredPoint<'a> {
    pub fn from(point: &'a Point) -> ClusteredPoint<'a> {
        ClusteredPoint {
            point,
            neighbors: vec![],
            cluster_id: 0,
        }
    }

    pub fn neighbors(&'a self) -> &'a Vec<&'a ClusteredPoint<'a>> {
        &self.neighbors
    }

    pub fn set_neighbor_capacity(&'a mut self, capacity: usize) -> &'a mut ClusteredPoint {
        self.neighbors.reserve(capacity);
        self
    }

    pub fn clear_neighbors(&'a mut self) -> &'a mut ClusteredPoint {
        self.neighbors.clear();
        self
    }

    pub fn core_distance(&self) -> f64 {
        return self.distance_to(self.neighbors.last().unwrap());
    }

    pub fn distance_to(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }

    pub fn distance_to_sqr(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }
}

pub struct ClustererBuilder<'a> {
    // data
    point_cloud: &'a [Point],
    // parameters
    leaf_size: usize,
}

impl<'a> ClustererBuilder<'a> {
    pub fn new(data: &'a [Point]) -> ClustererBuilder {
        ClustererBuilder {
            point_cloud: data,
            leaf_size: 50,
        }
    }

    pub fn build<'b>(&'b mut self) -> Clusterer<'a> {
        let mut data_refs = self.point_cloud.iter().collect();
        Clusterer {
            spatial_index_root: BallTree::new(&mut data_refs, self.leaf_size),
        }
    }

    pub fn with_leaf_size(&mut self, leaf_size: usize) -> &'a mut ClustererBuilder {
        self.leaf_size = leaf_size;
        self
    }
}

pub struct Clusterer<'a> {
    spatial_index_root: Box<BallTree<'a>>,
}

impl<'a> Clusterer<'a> {
    pub fn builder(data: &'a [Point]) -> ClustererBuilder {
        ClustererBuilder::new(data)
    }
}
