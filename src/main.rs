#![allow(dead_code)]

use rand::prelude::*;

use std::time::Instant;

mod ball_tree;
mod clusterer;
mod point;

use clusterer::Clusterer;
use point::Point;

use std::rc::Rc;

use crate::clusterer::ClusterResult;

fn main() {
    let num_points: usize = 100_000;
    let mut rng = rand::thread_rng();

    let mut data: Vec<Rc<Point>> = vec![];
    data.reserve(num_points);

    for _ in 0..num_points {
        data.push(Rc::new(Point::from(vec![
            rng.gen::<f64>(),
            rng.gen::<f64>(),
        ])));
    }

    let start = Instant::now();
    let leaf_size = 100;
    let param_k = 15;

    let clusterer = Clusterer::new(&data)
        .with_leaf_size(leaf_size)
        .with_param_k(param_k);

    let cluster_result: ClusterResult = clusterer.fit();

    let elapsed_sec = start.elapsed().as_secs();

    let ball_tree = cluster_result.ball_tree();
    let mut counter = 0;

    for _point_ref in ball_tree.iter() {
        counter += 1;
    }

    let average_core_distance = ball_tree
        .iter()
        .map(|x| x.borrow().core_distance())
        .sum::<f64>()
        / (counter as f64);

        
    

    println!(
        "setting up the ball tree for {} points took {} s",
        counter, elapsed_sec
    );

    println!(
        "the average core distance for param k = {} is {}.",
        param_k, average_core_distance
    );

}
