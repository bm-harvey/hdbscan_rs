#![allow(dead_code)]
#![allow(unused_variables)]

mod clusterer;
mod point;

use clusterer::Clusterer;
use point::Point;

fn main() {
    println!("Hello, world!");

    let mut data: Vec<Point> = vec![];
    for _ in 0..10 {
        data.push(Point::from(vec![1., 2., 3.]));
    }

    let mut clust = Clusterer::new();
    clust.set_leaf_size(45);
}
