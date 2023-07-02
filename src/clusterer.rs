use crate::point::Point;

pub struct ClusteredPoint<'a> {
    point: &'a Point,
    cluster_id: u16,
}

impl<'a> ClusteredPoint<'a> {
    pub fn from(point: &'a Point) -> ClusteredPoint<'a> {
        ClusteredPoint {
            point,
            cluster_id: 0,
        }
    }

    pub fn distance_to(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }

    pub fn distance_to_sqr(&self, other: &ClusteredPoint) -> f64 {
        self.point.distance_to(other.point)
    }
}

pub struct Clusterer<'a> {
    points: Vec<ClusteredPoint<'a>>,

    leaf_size: usize,
}

impl<'a> Clusterer<'a> {
    pub fn new() -> Clusterer<'a> {
        Clusterer {
            points: vec![],
            leaf_size: 50,
        }
    }

    pub fn from_data(data: &'a Vec<Point>) -> Clusterer {
        Clusterer {
            points: data
                .iter()
                .map(|point| ClusteredPoint::from(point))
                .collect(),

            leaf_size: 50,
        }
    }

    pub fn set_leaf_size(&'a mut self, leaf_size: usize) {
        self.leaf_size = leaf_size;
    }

    pub fn count(&self) -> usize {
        self.points.len()
    }
}
