#![allow(dead_code)]
#![allow(unused_variables)]

mod ball_tree;
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

    let clust = Clusterer::builder(&data).with_leaf_size(30).build();
}
