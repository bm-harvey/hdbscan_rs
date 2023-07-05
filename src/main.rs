#![allow(dead_code)]

use rand::prelude::*;

use std::time::Instant;

mod ball_tree;
mod clusterer;
mod point;

use clusterer::Clusterer;
use point::Point;

fn main() {
    let num_points: usize = 10000000;
    let mut rng = rand::thread_rng();

    let mut data: Vec<Point> = vec![];
    data.reserve(num_points);
    for _ in 0..num_points {
        data.push(Point::from(vec![rng.gen::<f64>(), rng.gen::<f64>()]));
    }

    let start = Instant::now();
    let leaf_size = 40* (num_points as f64).log2() as usize;
    println!("leaf_size = {}", leaf_size);
    let cluster_result = Clusterer::new(&data).with_leaf_size(leaf_size).fit();

    let elapsed_sec = start.elapsed().as_secs();

    let ball_tree = cluster_result.ball_tree();
    let mut counter = 0;
    for _ in ball_tree.iter() {
        counter += 1;
    }

    println!("setting up the ball tree for {} points took {} s", counter, elapsed_sec);
}
