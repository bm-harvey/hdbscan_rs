#![allow(dead_code)]
#![allow(unused_variables)]

mod clusterer;
mod ball_tree;
mod point;

use clusterer::Clusterer;
use point::Point;

fn main() {
    println!("Hello, world!");

    let mut data: Vec<Point> = vec![];
    for _ in 0..10 {
        data.push(Point::from(vec![1., 2., 3.]));
    }

    let clust = Clusterer::new(&data, 50);
}
