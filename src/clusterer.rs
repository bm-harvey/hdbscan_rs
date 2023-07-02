use crate::ball_tree::BallTree;
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
    spatial_index_root: Box<BallTree<'a>>,
}

impl<'a> Clusterer<'a> {
    pub fn new(data: &'a [Point], leaf_size: usize) -> Clusterer {
    //pub fn new(data: &'a Vec<Point>, leaf_size: usize) -> Clusterer {
        
        let mut data_ref : Vec<&Point> = vec![];
        for point in data.iter() {
            data_ref.push(point);
        }
            
        Clusterer {
            points: data
                .iter()
                .map(ClusteredPoint::from)
                .collect(),
            
            spatial_index_root: BallTree::create_ball_tree(&mut data_ref, leaf_size),
        }
    }

    pub fn count(&self) -> usize {
        self.points.len()
    }
}
