// use colored::Colorize;
// use math::round;
// use std::cmp;
use std::convert::TryFrom;
// use std::fmt;
use structopt::StructOpt;

mod cli_args;
mod requirements;
mod solutions;

use cli_args::Opt;
use requirements::Requirements;
use solutions::SolutionSet;

//----
// Main

fn main() {
    if let Ok(reqs) = Requirements::try_from(Opt::from_args()) {
        println!("{}", SolutionSet::from(reqs));
    }
}
