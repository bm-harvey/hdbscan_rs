use std::cmp::min;
use std::rc::Rc;

#[derive(Debug)]
pub struct Point {
    coordinate: Vec<f64>,
}

impl Point {
    pub fn new() -> Point {
        Point { coordinate: vec![] }
    }

    pub fn from(coordinate: Vec<f64>) -> Point {
        Point { coordinate }
    }

    pub fn num_dimensions(&self) -> usize {
        self.coordinate.len()
    }

    pub fn get(&self, idx: usize) -> f64 {
        self.coordinate[idx]
    }

    pub fn set(&mut self, idx: usize, value: f64) {
        self.coordinate[idx] = value;
    }

    pub fn distance_to(&self, other: &Point) -> f64 {
        let dims = min(self.num_dimensions(), other.num_dimensions());

        (0..dims)
            .map(|idx| self.coordinate[idx] - other.coordinate[idx])
            .map(|diff| diff.powi(2))
            .sum::<f64>()
            .sqrt()
    }

    pub fn distance_to_sqr(&self, other: &Rc<Point>) -> f64 {
        let dims = min(self.num_dimensions(), other.num_dimensions());

        (0..dims)
            .map(|idx| self.coordinate[idx] - other.coordinate[idx])
            .map(|diff| diff.powi(2))
            .sum::<f64>()
    }
}
