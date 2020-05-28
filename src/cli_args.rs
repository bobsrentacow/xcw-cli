use structopt::StructOpt;

//----
// Command Line Parsing

#[derive(Debug, StructOpt)]
#[structopt(
    name = "xilinx clock wizard",
    about = "How to calulate and display results"
)]
pub struct Opt {
    /// Select UltraScale MMCM mode
    #[structopt(short = "m", long = "mmcm")]
    pub use_mmcm: bool,
    /// Select UltraScale PLL mode
    #[structopt(short = "p", long = "pll")]
    pub use_pll: bool,

    /// Sort by root_mean_square_error
    #[structopt(short = "r", long = "rmse")]
    pub sort_by_rmse: bool,
    /// Sort by worst_ratiometric_error
    #[structopt(short = "w", long = "worst")]
    pub sort_by_worst: bool,

    /// Specify number of displayed solutions
    #[structopt(short = "n", long = "max_solutions", default_value = "32")]
    pub max_solutions: usize,

    /// Input frequency
    #[structopt(name = "inp_MHz")]
    pub inp_megahz: f64,

    /// Output frequencies
    #[structopt(name = "out_MHz")]
    pub output_specifiers: Vec<String>,
}

