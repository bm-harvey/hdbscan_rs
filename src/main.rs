use hdbscan_rs::{BallTree, ClusteredPoint};
use rand::prelude::*;

use std::time::Instant;
use std::usize;

use hdbscan_rs::Clusterer;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(
    name = "hdbscan",
    author = "Bryan Harvey",
    about = "a clustering routine"
)]
struct Arg {
    ///Number of points to generate
    #[structopt(short, long, default_value = "1000")]
    num_samples: usize,

    ///Number of nearest neighbors to find
    #[arg(short = 'k', long = "param_k", default_value = "5")]
    param_k: usize,

    ///Size of leaf nodes in ball trees
    #[arg(short, long, default_value = "50")]
    leaf_size: usize,
}

fn main() {
    let args = Arg::parse();

    let mut rng = rand::thread_rng();

    let mut data = vec![];

    for _ in 0..args.num_samples / 2 {
        data.push(ClusteredPoint::from_as_rcc(vec![
            rng.gen::<f64>() * 10.,
            rng.gen::<f64>(),
        ]));
    }

    for _ in 0..args.num_samples / 2 {
        data.push(ClusteredPoint::from_as_rcc(vec![
            rng.gen::<f64>() + 5.,
            rng.gen::<f64>() + 5.,
        ]));
    }

    let start = Instant::now();

    let clusterer = Clusterer::new(data)
        .with_leaf_size(args.leaf_size)
        .with_param_k(args.param_k);

    let cluster_result = clusterer.fit();

    let elapsed_sec = start.elapsed().as_secs_f32();

    let ball_tree = cluster_result.ball_tree();
    let mut counter = 0;

    for _point_ref in ball_tree.cluster_points() {
        counter += 1;
    }

    let average_core_distance = ball_tree
        .cluster_points()
        .map(|x| BallTree::core_distance(&x, hdbscan_rs::Metric::Euclidean))
        .sum::<f64>()
        / (counter as f64);

    println!(
        "setting up the ball tree for {} points took {} s",
        counter, elapsed_sec
    );

    println!(
        "the average core distance for param k = {} is {}.",
        args.param_k, average_core_distance
    );
}
