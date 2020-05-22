use colored::Colorize;
use math::round;
use std::cmp;
use std::convert::TryFrom;
use std::fmt;
use structopt::StructOpt;

// TODO: Learn how to use the doc tool and add special comments

#[derive(Debug, StructOpt)]
#[structopt(name = "xilinx clock wizard", about = "How to calulate and display results")]
struct Opt {
    /// Select UltraScale MMCM mode
    #[structopt(short = "m", long = "mmcm")]
    use_mmcm: bool,
    /// Select UltraScale PLL mode
    #[structopt(short = "p", long = "pll")]
    use_pll: bool,

    /// Sort by root_mean_square_error
    #[structopt(short = "r", long = "rmse")]
    sort_by_rmse: bool,
    /// Sort by worst_ratiometric_error
    #[structopt(short = "w", long = "worst")]
    sort_by_worst: bool,

    /// Specify number of displayed solutions
    #[structopt(short = "n", long = "max_solutions", default_value = "32")]
    max_solutions: usize,

    /// Input frequency
    #[structopt(name = "inp_MHz")]
    inp_megahz: f64,

    /// Output frequencies
    #[structopt(name = "out_MHz")]
    output_specifiers: Vec<String>,
}

//----
// OutputConstraint

enum OutputConstraint {
    Normal(f64),
    Range { min: f64, max: f64 },
    LessThan(f64),
    LessThanOrEqual(f64),
    Equal(f64),
    GreaterThanOrEqual(f64),
    GreaterThan(f64),
}
impl OutputConstraint {
    fn try_parse(the_string: &str) -> Result<Self, &'static str> {
        // text_io parsers
        let try_scan_less_than_or_equal =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "lte{}", target);
                Ok(OutputConstraint::LessThanOrEqual(target))
            };
        let try_scan_less_than =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "lt{}", target);
                Ok(OutputConstraint::LessThan(target))
            };
        let try_scan_equal =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "eq{}", target);
                Ok(OutputConstraint::Equal(target))
            };
        let try_scan_greater_than_or_equal =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "gte{}", target);
                Ok(OutputConstraint::GreaterThanOrEqual(target))
            };
        let try_scan_greater_than =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "gt{}", target);
                Ok(OutputConstraint::GreaterThan(target))
            };
        let try_scan_plus_minus_ppm =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                let tolerance: f64;
                text_io::try_scan!(string.bytes() => "{}+-{}ppm", target, tolerance);
                Ok(OutputConstraint::Range {
                    min: target * (1_f64 - 1e-6 * tolerance),
                    max: target * (1_f64 + 1e-6 * tolerance),
                })
            };
        let try_scan_plus_minus_percent =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                let tolerance: f64;
                text_io::try_scan!(string.bytes() => "{}+-{}pct", target, tolerance);
                Ok(OutputConstraint::Range {
                    min: target * (1_f64 - 1e-2 * tolerance),
                    max: target * (1_f64 + 1e-2 * tolerance),
                })
            };
        let try_scan_plus_minus =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                let tolerance: f64;
                text_io::try_scan!(string.bytes() => "{}+-{}", target, tolerance);
                Ok(OutputConstraint::Range {
                    min: target - tolerance,
                    max: target + tolerance,
                })
            };
        let try_scan_range = |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
            let min: f64;
            let max: f64;
            text_io::try_scan!(string.bytes() => "{}-{}", min, max);
            Ok(OutputConstraint::Range { min, max })
        };
        let try_scan_normal = |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
            let target: f64;
            text_io::try_scan!(string.bytes() => "{}", target);
            Ok(OutputConstraint::Normal(target))
        };

        if let Ok(res) = try_scan_less_than_or_equal(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_less_than(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_equal(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_greater_than_or_equal(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_greater_than(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_plus_minus_ppm(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_plus_minus_percent(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_plus_minus(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_range(&the_string) {
            Ok(res)
        } else if let Ok(res) = try_scan_normal(&the_string) {
            Ok(res)
        } else {
            Err("invalid output specifier")
        }
    }
}

//----
// Requirements

#[derive(Debug, Copy)]
struct Fraction {
    num: u16,
    den: u16,
}
impl Fraction {
    fn into_f64(self) -> f64 {
        (self.num as f64) / (self.den as f64)
    }
}
impl Clone for Fraction {
    fn clone(&self) -> Fraction {
        *self
    }
}

#[derive(Debug)]
#[allow(dead_code)] // Absolute is never constructed
enum ErrorType {
    Absolute,
    Ratiometric,
}
impl fmt::Display for ErrorType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorType::Absolute    => write!(f, "Absolute"),
            ErrorType::Ratiometric => write!(f, "Ratiometric"),
        }
    }
}

#[derive(Debug)]
enum SortOrder {
    RootMeanSquareError,
    RatiometricErrorWorstChannel,
    RatiometricErrorOnChannel(u8),
}
impl fmt::Display for SortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SortOrder::RootMeanSquareError           => write!(f, "RootMeanSquareError"),
            SortOrder::RatiometricErrorWorstChannel  => write!(f, "RatiometricErrorWorstChannel"),
            SortOrder::RatiometricErrorOnChannel(ch) => write!(f, "RatiometricErrorOnChannel({})", ch),
        }
    }
}

struct Requirements {
    error_type: ErrorType,
    sort_order: SortOrder,

    max_solutions: usize,

    vco_divider_num_min: Fraction,
    vco_divider_num_max: Fraction,
    vco_divider_den_min: u16,
    vco_divider_den_max: u16,
    vco_megahz_max: f64,
    vco_megahz_min: f64,
    chan_divider_min: Vec<Fraction>,
    chan_divider_max: Vec<Fraction>,

    inp_megahz: f64,
    output_constraints: Vec<OutputConstraint>,
}
impl TryFrom<Opt> for Requirements {
    type Error = &'static str;

    fn try_from(opt: Opt) -> Result<Self, Self::Error> {
        let mut error_type = ErrorType::Ratiometric;
        let mut sort_order = SortOrder::RootMeanSquareError;

        let mut max_outputs = 0_usize;
        let mut vco_divider_num_min = Fraction { num: 0, den: 0 };
        let mut vco_divider_num_max = Fraction { num: 0, den: 0 };
        let mut vco_divider_den_min = 0_u16;
        let mut vco_divider_den_max = 0_u16;
        let mut vco_megahz_max = 0_f64;
        let mut vco_megahz_min = 0_f64;
        let mut chan_divider_min = Vec::<Fraction>::new();
        let mut chan_divider_max = Vec::<Fraction>::new();

        if opt.sort_by_rmse && opt.sort_by_worst {
            return Err("Can't specify two different sort orders");
        } else if opt.sort_by_rmse {
            sort_order = SortOrder::RootMeanSquareError;
        } else if opt.sort_by_worst {
            sort_order = SortOrder::RatiometricErrorWorstChannel;
        }

        // What hardware are we targetting?  Apply limits associated with that hardware.
        if opt.use_mmcm && opt.use_pll {
            return Err("must specify exactly one target");
        } else if opt.use_mmcm {
            max_outputs = 8;
            vco_divider_num_min = Fraction { num:  2 * 8, den: 8 };
            vco_divider_num_max = Fraction { num: 64 * 8, den: 8 };
            vco_divider_den_min =   1;
            vco_divider_den_max = 106;
            vco_megahz_max = 1200_f64;
            vco_megahz_min =  600_f64;
            chan_divider_min.push(Fraction { num:   1 * 8, den: 8});
            chan_divider_max.push(Fraction { num: 128 * 8, den: 8});
            for _ in 1..max_outputs {
                chan_divider_min.push(Fraction { num:   1, den: 1});
                chan_divider_max.push(Fraction { num: 128, den: 1});
            }
        } else if opt.use_pll {
            max_outputs = 2;
            vco_divider_num_min = Fraction { num:  2 * 8, den: 8 };
            vco_divider_num_max = Fraction { num: 64 * 8, den: 8 };
            vco_divider_den_min =   1;
            vco_divider_den_max = 106;
            vco_megahz_max = 1200_f64;
            vco_megahz_min =  600_f64;
            for _ in 0..max_outputs {
                chan_divider_min.push(Fraction { num: 128, den: 1});
                chan_divider_max.push(Fraction { num: 128, den: 1});
            }
        } else {
            return Err("must specify exactly one target");
        }

        // Check limits associated with target.
        if opt.output_specifiers.len() > max_outputs {
            return Err("too many outputs requested");
        }

        let mut output_constraints = Vec::<OutputConstraint>::new();

        // Figure out what constraints to apply for each requested output.
        for (ii, the_string) in opt.output_specifiers.iter().enumerate() {
            match OutputConstraint::try_parse(&the_string) {
                Ok(constraint) => {
                    if let OutputConstraint::Equal(_) = constraint {
                        if opt.sort_by_rmse || opt.sort_by_worst {
                            return Err("Can't constrain an output to eq when a sort order is specified");
                        }
                        sort_order = SortOrder::RatiometricErrorOnChannel(ii as u8);
                    }
                    output_constraints.push(constraint)
                },
                Err(err_str) => return Err(&err_str),
            }
        }

        Ok(Requirements {
            error_type,
            sort_order,
            max_solutions: opt.max_solutions,

            vco_divider_num_min,
            vco_divider_num_max,
            vco_divider_den_min,
            vco_divider_den_max,
            vco_megahz_max,
            vco_megahz_min,
            chan_divider_min,
            chan_divider_max,

            inp_megahz: opt.inp_megahz,
            output_constraints,
        })
    }
}
impl fmt::Display for Requirements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "error_type         : {}", self.error_type)?;
        writeln!(f, "sort_order         : {}", self.sort_order)?;
        writeln!(f, "max_solutions      : {}", self.max_solutions)?;
        writeln!(f, "vco_divider_num_min: {}", self.vco_divider_num_min.into_f64())?;
        writeln!(f, "vco_divider_num_max: {}", self.vco_divider_num_max.into_f64())?;
        writeln!(f, "vco_divider_den_min: {}", self.vco_divider_den_min)?;
        writeln!(f, "vco_divider_den_max: {}", self.vco_divider_den_max)?;
        writeln!(f, "vco_megahz_max     : {}", self.vco_megahz_max)?;
        writeln!(f, "vco_megahz_min     : {}", self.vco_megahz_min)?;
        writeln!(f, "inp_megahz         : {}", self.inp_megahz)?;
        for (ii, constraint) in self.output_constraints.iter().enumerate() {
            match constraint {
                OutputConstraint::Normal(target) =>
                    writeln!(f, "output{:02}:   {}", ii, *target),
                OutputConstraint::Range { min, max } =>
                    writeln!(f, "output{:02}:   {}-{}", ii, *min, *max),
                OutputConstraint::LessThan(target) =>
                    writeln!(f, "output{:02}: < {}", ii, *target),
                OutputConstraint::LessThanOrEqual(target) =>
                    writeln!(f, "output{:02}: <={}", ii, *target),
                OutputConstraint::Equal(target) =>
                    writeln!(f, "output{:02}:  ={}", ii, *target),
                OutputConstraint::GreaterThanOrEqual(target) =>
                    writeln!(f, "output{:02}: >={}", ii, *target),
                OutputConstraint::GreaterThan(target) =>
                    writeln!(f, "output{:02}: > {}", ii, *target),
            }?;
        }
        Ok(())
    }
}

//----
// SolutionSet

#[derive(Debug)]
struct VcoSolution {
    input: f64,
    output: f64,
    numerator: Fraction,
    denominator: u16,
}
impl VcoSolution {
    fn get_solutions(reqs: &Requirements) -> Vec<VcoSolution> {
        let vco_divider_num_den = reqs.vco_divider_num_min.den;
        let vco_divider_num_den_f64 = vco_divider_num_den as f64;

        let in_num_min = cmp::max(
            reqs.vco_divider_num_min.num,
            round::ceil(reqs.vco_megahz_min * (reqs.vco_divider_den_min as f64) * vco_divider_num_den_f64 / reqs.inp_megahz, 0) as u16,
        );
        let in_num_max = cmp::min(
            reqs.vco_divider_num_max.num,
            round::floor(reqs.vco_megahz_max * (reqs.vco_divider_den_max as f64) * vco_divider_num_den_f64 / reqs.inp_megahz, 0) as u16,
        );
        //println!("in_num_min {}, in_num_max {}", in_num_min, in_num_max);

        let mut vco_solns = Vec::<VcoSolution>::new();
        for in_num in in_num_min..=in_num_max {
            let vco_divider_num_f64 = (in_num as f64) / vco_divider_num_den_f64;

            let in_den_min = cmp::max(
                reqs.vco_divider_den_min,
                round::ceil(
                    reqs.inp_megahz * vco_divider_num_f64 / reqs.vco_megahz_max,
                    0,
                ) as u16,
            );
            let in_den_max = cmp::min(
                reqs.vco_divider_den_max,
                round::floor(
                    reqs.inp_megahz * vco_divider_num_f64 / reqs.vco_megahz_min,
                    0,
                ) as u16,
            );
            //println!("in_num {}, in_den_min {}, in_den_max {}", in_num, in_den_min, in_den_max);

            for in_den in in_den_min..=in_den_max {
                let vco = reqs.inp_megahz * vco_divider_num_f64 / (in_den as f64);

                // check for out of bounds vco freq
                if vco < reqs.vco_megahz_min {
                    println!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco lo break\n",
                        reqs.inp_megahz,
                        vco_divider_num_f64,
                        in_den,
                        vco
                    );
                    break;
                }
                if vco > reqs.vco_megahz_max {
                    println!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco hi continue\n",
                        reqs.inp_megahz,
                        vco_divider_num_f64,
                        in_den,
                        vco
                    );
                    continue;
                }

                let thresh = 1_f64 + 1e-9; // magic number for vco frequency equality threshold
                let found = vco_solns.iter().find(|&x| ((x.output / vco) < thresh) && ((vco / x.output) < thresh));
                match found {
                    Some(_) => (),
                    None => {
                        vco_solns.push(VcoSolution{
                            input: reqs.inp_megahz,
                            output: vco,
                            numerator: Fraction{ num: in_num, den: vco_divider_num_den },
                            denominator: in_den,
                        });
                    }
                }
            }
        }

        // sort from high to low vco frequencies to reduce output jitter
        vco_solns.sort_by(|a, b| b.output.partial_cmp(&a.output).unwrap());

        vco_solns
    }
}

#[derive(Debug)]
struct ChannelSolution {
    input: f64,
    chan_idx: u8,
    divider: Fraction,
    target: f64,
    actual: f64,
    absolute_error: f64,
    ratiometric_error: f64,
}
impl ChannelSolution {
    fn try_solve(
        vco: f64,
        chan_idx: u8,
        constraint: &OutputConstraint,
        divider_min: Fraction,
        divider_max: Fraction
    ) -> Result<ChannelSolution, String> {
        let den = divider_min.den;
        let den_f64 = den as f64;
        let target = match constraint {
            OutputConstraint::Normal(target) => *target,
            OutputConstraint::Range { min, max } => (*min + *max) / 2_f64,
            OutputConstraint::LessThan(target) => *target,
            OutputConstraint::LessThanOrEqual(target) => *target,
            OutputConstraint::Equal(target) => *target,
            OutputConstraint::GreaterThanOrEqual(target) => *target,
            OutputConstraint::GreaterThan(target) => *target,
        };

        // get closest integer solution
        let num = (vco * den_f64 / target).round() as u16;
        // check num + -1..=+1 solutions
        let num_candidates = vec![num-1, num, num+1];

        let mut num_tuples = Vec::<(u16, f64, f64, f64, f64)>::new();
        for num in num_candidates {
            if (num < divider_min.num) || (num > divider_max.num) {
                //println!("vco {}, dev {:.3} -- dq {} < {} || {} > {}", vco, Fraction {num, den}.into_f64(), num, divider_min.num, num, divider_max.num);
                continue;
            }

            let actual = vco * (den_f64 / (num as f64));
            let error = actual - target;
            let absolute_error = error.abs();
            let ratiometric_error = (error / target).abs();

            //println!("vco {}, dev {:.3}, out {}", vco, Fraction {num, den}.into_f64(), actual);

            // check output range constraints
            match constraint {
                OutputConstraint::Range { min, max } => {
                    if (actual < *min) || (actual > *max) {
                        continue;
                    }
                }
                OutputConstraint::LessThan(_) => {
                    if actual >= target {
                        continue;
                    }
                }
                OutputConstraint::LessThanOrEqual(_) => {
                    if actual > target {
                        continue;
                    }
                }
                OutputConstraint::GreaterThanOrEqual(_) => {
                    if actual < target {
                        continue;
                    }
                }
                OutputConstraint::GreaterThan(_) => {
                    if actual <= target {
                        continue;
                    }
                }
                _ => ()
            }
            num_tuples.push((num, actual, error, absolute_error, ratiometric_error));
        }

        if !num_tuples.is_empty() {
            num_tuples.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap());
            let (num, actual, error, absolute_error, ratiometric_error) = num_tuples[0];

            Ok(ChannelSolution {
                input: vco,
                chan_idx,
                divider: Fraction { num, den },
                target,
                actual,
                absolute_error,
                ratiometric_error,
            })
        } else {
            Err("No solutions found".to_string())
        }
    }
}

#[derive(Debug)]
struct Solution {
    vco_solution: VcoSolution,
    channel_solutions: Vec<ChannelSolution>,
    root_mean_square_error: f64,
    worst_error: f64,
    channel_with_worst_error: u8,
}

#[derive(Debug)]
struct SolutionSet {
    error_type: ErrorType,
    sort_order: SortOrder,
    solutions: Vec<Solution>,
}

impl SolutionSet {
    fn from(reqs: Requirements) -> Self {
        let vco_solns = VcoSolution::get_solutions(&reqs);
        let mut solutions = Vec::<Solution>::new();
         'vco: for VcoSolution {
            input,
            output: vco_freq,
            numerator,
            denominator,
         } in &vco_solns {
            //println!("vco {:4.6}, numerator {}, denominator {}, ", vco_freq, numerator.into_f64(), denominator);

            // compute nearest output frequencies
            let mut channel_solutions = Vec::<ChannelSolution>::new();
            let mut mse = 0_f64;
            let mut max_err = -1_f64;
            let mut max_err_chan = 0; 
            for (chan, constraint) in reqs.output_constraints.iter().enumerate() {
                match ChannelSolution::try_solve(
                    *vco_freq,
                    chan as u8,
                    &constraint,
                    reqs.chan_divider_min[chan],
                    reqs.chan_divider_max[chan]
                ) {
                    Ok(soln) => {
                        mse += soln.absolute_error * soln.absolute_error;
                        let err = match reqs.error_type {
                            ErrorType::Absolute => soln.absolute_error,
                            ErrorType::Ratiometric => soln.ratiometric_error,
                        };
                        if max_err < err {
                            max_err = err;
                            max_err_chan = chan;
                        }
                        channel_solutions.push(soln);
                    },
                    Err(_) => continue 'vco,
                }
            }

            //---- Solution is valid ----

            solutions.push(Solution {
                vco_solution: VcoSolution {
                    input: *input,
                    output: *vco_freq,
                    numerator: *numerator,
                    denominator: *denominator,
                },
                channel_solutions,
                root_mean_square_error: mse.sqrt(),
                worst_error: max_err,
                channel_with_worst_error: max_err_chan as u8,
            });
        }

        //---- Sort and trim ----
        match reqs.sort_order {
            SortOrder::RootMeanSquareError => {
                solutions.sort_by(|a, b| a.root_mean_square_error.partial_cmp(&b.root_mean_square_error).unwrap());
            }
            SortOrder::RatiometricErrorWorstChannel => {
                solutions.sort_by(|a, b| a.worst_error.partial_cmp(&b.worst_error).unwrap());
            }
            SortOrder::RatiometricErrorOnChannel(ch) => {
                solutions.sort_by(|a, b|
                    a.channel_solutions[ch as usize].ratiometric_error.partial_cmp(
                        &b.channel_solutions[ch as usize].ratiometric_error
                    ).unwrap()
               );
            }
        };
        solutions.truncate(reqs.max_solutions);

        SolutionSet {
            error_type: reqs.error_type,
            sort_order: reqs.sort_order,
            solutions,
        }
    }
}

impl fmt::Display for SolutionSet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.solutions.is_empty() {
            write!(f, "No solutions found")
        } else {
            let num_outputs = self.solutions[0].channel_solutions.len();

            //-- Table Annotation
            match self.sort_order {
                SortOrder::RootMeanSquareError =>
                    writeln!(f, "Sorting in order of increasing root_mean_square_error"),
                SortOrder::RatiometricErrorWorstChannel =>
                    writeln!(f, "Sorting in order of increasing ppm_err_max"),
                SortOrder::RatiometricErrorOnChannel(ch) =>
                    writeln!(f, "Sorting in order of increasing error on channel {}", ch),
            }?;
            match self.error_type {
                ErrorType::Absolute => writeln!(f, "Worst absolute output in {}", "red".red()),
                ErrorType::Ratiometric => writeln!(f, "Worst ratiomnetric output in {}", "red".red()),
            }?;
            writeln!(f)?;

            //-- Header
            write!(f, "{:>5} {:>12}", "sln#", "vco")?;
            for ii in 0..num_outputs {
                write!(f, " {:>12}{:1}", "MHz", ii)?;
            }
            match self.error_type {
                ErrorType::Absolute => write!(f, " {:>13}", "MHz_err_max"),
                ErrorType::Ratiometric => write!(f, " {:>13}", "ppm_err_max"),
            }?;
            write!(f, " {:>13}", "rms_err(MHz)")?;
            write!(f, " {:>6}", "clkfb")?;
            write!(f, " {:>6}", "divclk")?;
            for ii in 0..num_outputs {
                write!(f, " {:>5}{:1}", "odiv", ii)?;
            }

            //-- Solutions
            for (ii, soln) in self.solutions.iter().enumerate() {
                writeln!(f)?;
                write!(f, "{:>5} {:>12.6}", ii, soln.vco_solution.output)?;
                for chan_soln in &soln.channel_solutions {
                    let mut str_megahz: String = format!(" {:>13.6}", chan_soln.actual);
                    if ii == (soln.channel_with_worst_error as usize) {
                        str_megahz = str_megahz.red().to_string();
                    }
                    write!(f, "{}", str_megahz)?;
                }
                write!(f, " {:>13.3}", 1e6 * soln.worst_error)?;
                write!(f, " {:>13.6}", soln.root_mean_square_error)?;
                write!(f, " {:>6.3}", soln.vco_solution.numerator.into_f64())?;
                write!(f, " {:>6}", soln.vco_solution.denominator as f64)?;
                for chan_soln in &soln.channel_solutions {
                    write!(f, " {:6.3}", chan_soln.divider.into_f64())?;
                }
            }

            write!(f, "")
        }
    }
}

fn main() {
    if let Ok(reqs) = Requirements::try_from(Opt::from_args()) {
        println!("{}", SolutionSet::from(reqs));
    }
}

// TODO: fix tests for the refactor that happened a few commits ago
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_requirements_from_opt_use_mmcm_and_pll() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: true,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("gte187.5")],
        };
        let reqs = Requirements::from(opt);
        assert_eq!(reqs.valid, false);
        Ok(())
    }

    #[test]
    fn test_requirements_from_opt_use_mmcm() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("gte187.5")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_constraints.len();

        let reqs = Requirements::from(opt);

        assert_eq!(reqs.valid           , true);
        assert_eq!(reqs.max_solutions   , max_solutions);
        assert_eq!(reqs.max_outputs     ,        8);
        assert_eq!(reqs.vco_megahz_max  , 1200_f64);
        assert_eq!(reqs.vco_megahz_min  ,  600_f64);
        assert_eq!(reqs.inp_megahz      , inp_megahz);
        assert_eq!(reqs.output_constraints.len(), num_outputs);
        Ok(())
    }

    #[test]
    fn test_requirements_from_opt_use_pll() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("gte187.5")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_constraints.len();

        let reqs = Requirements::from(opt);

        assert_eq!(reqs.valid           , true);
        assert_eq!(reqs.max_solutions   , max_solutions);
        assert_eq!(reqs.max_outputs     ,        2);
        assert_eq!(reqs.vco_megahz_max  , 1200_f64);
        assert_eq!(reqs.vco_megahz_min  ,  600_f64);
        assert_eq!(reqs.inp_megahz      , inp_megahz);
        assert_eq!(reqs.output_constraints.len(), num_outputs);
        Ok(())
    }

    #[test]
    fn test_requirements_from_opt_range() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("181.1-188.8")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_constraints.len();

        let reqs = Requirements::from(opt);

        assert_eq!(reqs.valid           , true);
        assert_eq!(reqs.max_solutions   , max_solutions);
        assert_eq!(reqs.max_outputs     ,        2);
        assert_eq!(reqs.vco_megahz_max  , 1200_f64);
        assert_eq!(reqs.vco_megahz_min  ,  600_f64);
        assert_eq!(reqs.inp_megahz      , inp_megahz);
        assert_eq!(reqs.output_constraints.len(), num_outputs);
        if let OutputConstraint::Range { min, max, } = reqs.output_constraints[1] {
            assert!((min - 181.1).abs() < 1e-6);
            assert!((max - 188.8).abs() < 1e-6);
            Ok(())
        } else {
            Err(String::from("parsing error"))
        }
    }

    #[test]
    fn test_requirements_from_opt_range_pct() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("181.1+-0.9pct")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_constraints.len();

        let reqs = Requirements::from(opt);

        assert_eq!(reqs.valid           , true);
        assert_eq!(reqs.max_solutions   , max_solutions);
        assert_eq!(reqs.max_outputs     ,        2);
        assert_eq!(reqs.vco_megahz_max  , 1200_f64);
        assert_eq!(reqs.vco_megahz_min  ,  600_f64);
        assert_eq!(reqs.inp_megahz      , inp_megahz);
        assert_eq!(reqs.output_constraints.len(), num_outputs);
        if let OutputConstraint::Range { min, max, } = reqs.output_constraints[1] {
            assert!((min - (181.1 * 0.991)).abs() < 1e-6);
            assert!((max - (181.1 * 1.009)).abs() < 1e-6);
            Ok(())
        } else {
            Err(String::from("parsing error"))
        }
    }

    #[test]
    fn test_requirements_from_opt_range_ppm() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_constraints: vec![String::from("166.6"), String::from("181.1+-3.5ppm")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_constraints.len();

        let reqs = Requirements::from(opt);

        assert_eq!(reqs.valid           , true);
        assert_eq!(reqs.max_solutions   , max_solutions);
        assert_eq!(reqs.max_outputs     ,        2);
        assert_eq!(reqs.vco_megahz_max  , 1200_f64);
        assert_eq!(reqs.vco_megahz_min  ,  600_f64);
        assert_eq!(reqs.inp_megahz      , inp_megahz);
        assert_eq!(reqs.output_constraints.len(), num_outputs);
        if let OutputConstraint::Range { min, max, } = reqs.output_constraints[1] {
            assert!((min - (181.1 * (1_f64 - 3.5e-6))).abs() < 1e-6);
            assert!((max - (181.1 * (1_f64 + 3.5e-6))).abs() < 1e-6);
            Ok(())
        } else {
            Err(String::from("parsing error"))
        }
    }

    // TODO: add many tests
}
