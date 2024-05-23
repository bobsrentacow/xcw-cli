use colored::Colorize;
use math::round;
use std::cmp;
use std::fmt;
use structopt::StructOpt;

// TODO: Learn how to use the doc tool and add special comments

// TODO: make the --help option useful
#[derive(Debug, StructOpt)]
#[structopt(name = "cli options", about = "How to calulate and display results")]
struct Opt {
    /// Select MMCM mode
    #[structopt(short = "m", long = "mmcm")]
    use_mmcm: bool,

    /// Select PLL mode
    #[structopt(short = "p", long = "pll")]
    use_pll: bool,

    /// Specify number of displayed solutions
    #[structopt(short = "s", long = "num_solutions", default_value = "32")]
    num_solutions: usize,

    /// Input frequency
    #[structopt(name = "inp_MHz")]
    inp_megahz: f64,

    /// Output frequencies
    #[structopt(name = "out_MHz")]
    out_megahz: Vec<String>,
}

enum OutputConstraint {
    Normal(f64),
    Range { min: f64, max: f64 },
    LessThan(f64),
    LessThanOrEqual(f64),
    Equal(f64),
    GreaterThanOrEqual(f64),
    GreaterThan(f64),
}

struct Requirements {
    valid: bool,
    max_solutions: usize,

    max_outputs: u8,
    vco_megahz_max: f64,
    vco_megahz_min: f64,

    inp_megahz: f64,
    out_megahz: Vec<OutputConstraint>,
}

impl From<Opt> for Requirements {
    fn from(opt: Opt) -> Self {
        let mut reqs = Requirements {
            valid: false,
            max_solutions: opt.num_solutions,

            max_outputs: 0,
            vco_megahz_max: 0_f64,
            vco_megahz_min: 0_f64,

            inp_megahz: opt.inp_megahz,
            out_megahz: Vec::<OutputConstraint>::new(),
        };
        if opt.use_mmcm {
            reqs.valid = true;
            reqs.max_outputs = 8;
            reqs.vco_megahz_max = 1200_f64;
            reqs.vco_megahz_min =  600_f64;
        } else if opt.use_pll {
            reqs.valid = true;
            reqs.max_outputs = 2;
            reqs.vco_megahz_max = 1200_f64;
            reqs.vco_megahz_min =  600_f64;
        } else {
            reqs.valid = false;
            return reqs;
        }
        if opt.out_megahz.len() > reqs.max_outputs as usize {
            reqs.valid = false;
            return reqs;
        }
        for the_string in opt.out_megahz.into_iter() {
            if the_string.starts_with("lte") {
                if let Ok(megahz) = the_string[3..].parse::<f64>() {
                    reqs.out_megahz
                        .push(OutputConstraint::LessThanOrEqual(megahz));
                } else {
                    reqs.valid = false;
                    return reqs;
                }
            } else if the_string.starts_with("lt") {
                if let Ok(megahz) = the_string[2..].parse::<f64>() {
                    reqs.out_megahz.push(OutputConstraint::LessThan(megahz));
                } else {
                    reqs.valid = false;
                    return reqs;
                }
            } else if the_string.starts_with("eq") {
                if let Ok(megahz) = the_string[2..].parse::<f64>() {
                    reqs.out_megahz.push(OutputConstraint::Equal(megahz));
                } else {
                    reqs.valid = false;
                    return reqs;
                }
            } else if the_string.starts_with("gte") {
                if let Ok(megahz) = the_string[3..].parse::<f64>() {
                    reqs.out_megahz
                        .push(OutputConstraint::GreaterThanOrEqual(megahz));
                } else {
                    reqs.valid = false;
                    return reqs;
                }
            } else if the_string.starts_with("gt") {
                if let Ok(megahz) = the_string[2..].parse::<f64>() {
                    reqs.out_megahz.push(OutputConstraint::GreaterThan(megahz));
                } else {
                    reqs.valid = false;
                    return reqs;
                }
            } else {
                let scan_plus_minus_ppm =
                    |string: &str| -> Result<(f64, f64), Box<dyn std::error::Error>> {
                        let target: f64;
                        let tolerance: f64;
                        text_io::try_scan!(string.bytes() => "{}+-{}ppm", target, tolerance);
                        Ok((target, tolerance))
                    };
                let scan_plus_minus_percent =
                    |string: &str| -> Result<(f64, f64), Box<dyn std::error::Error>> {
                        let target: f64;
                        let tolerance: f64;
                        text_io::try_scan!(string.bytes() => "{}+-{}pct", target, tolerance);
                        Ok((target, tolerance))
                    };
                let scan_plus_minus =
                    |string: &str| -> Result<(f64, f64), Box<dyn std::error::Error>> {
                        let target: f64;
                        let tolerance: f64;
                        text_io::try_scan!(string.bytes() => "{}+-{}", target, tolerance);
                        Ok((target, tolerance))
                    };
                let scan_range = |string: &str| -> Result<(f64, f64), Box<dyn std::error::Error>> {
                    let min: f64;
                    let max: f64;
                    text_io::try_scan!(string.bytes() => "{}-{}", min, max);
                    Ok((min, max))
                };
                let scan_normal = |string: &str| -> Result<f64, Box<dyn std::error::Error>> {
                    let target: f64;
                    text_io::try_scan!(string.bytes() => "{}", target);
                    Ok(target)
                };
                if let Ok((target, tolerance)) = scan_plus_minus_ppm(&the_string) {
                    reqs.out_megahz.push(OutputConstraint::Range {
                        min: target * (1_f64 - 1e-6 * tolerance),
                        max: target * (1_f64 + 1e-6 * tolerance),
                    });
                } else if let Ok((target, tolerance)) = scan_plus_minus_percent(&the_string) {
                    reqs.out_megahz.push(OutputConstraint::Range {
                        min: target * (1_f64 - 1e-2 * tolerance),
                        max: target * (1_f64 + 1e-2 * tolerance),
                    });
                } else if let Ok((target, tolerance)) = scan_plus_minus(&the_string) {
                    reqs.out_megahz.push(OutputConstraint::Range {
                        min: target - tolerance,
                        max: target + tolerance,
                    });
                } else if let Ok((min, max)) = scan_range(&the_string) {
                    reqs.out_megahz.push(OutputConstraint::Range { min, max });
                } else if let Ok(target) = scan_normal(&the_string) {
                    reqs.out_megahz.push(OutputConstraint::Normal(target));
                } else {
                    println!("ERROR: invalid output specifier: {}", the_string);
                    reqs.valid = false;
                    return reqs;
                }
            }
        }
        reqs
    }
}

impl fmt::Display for Requirements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "valid         : {}", self.valid)?;
        writeln!(f, "max_solutions : {}", self.max_solutions)?;
        writeln!(f, "max_outputs   : {}", self.max_outputs)?;
        writeln!(f, "vco_megahz_max: {}", self.vco_megahz_max)?;
        writeln!(f, "vco_megahz_min: {}", self.vco_megahz_min)?;
        writeln!(f, "inp_megahz    : {}", self.inp_megahz)?;
        for constraint in &self.out_megahz {
            match constraint {
                OutputConstraint::Normal(target) => writeln!(f, "out_megahz    :   {}", *target),
                OutputConstraint::Range {
                    min: target_min,
                    max: target_max,
                } => writeln!(f, "out_megahz    :   {}-{}", *target_min, *target_max),
                OutputConstraint::LessThan(target) => writeln!(f, "out_megahz    : < {}", *target),
                OutputConstraint::LessThanOrEqual(target) => {
                    writeln!(f, "out_megahz    : <={}", *target)
                }
                OutputConstraint::Equal(target) => writeln!(f, "out_megahz    :  ={}", *target),
                OutputConstraint::GreaterThanOrEqual(target) => {
                    writeln!(f, "out_megahz    : >={}", *target)
                }
                OutputConstraint::GreaterThan(target) => {
                    writeln!(f, "out_megahz    : > {}", *target)
                }
            }?;
        }
        Ok(())
    }
}

#[derive(Debug)]
#[allow(dead_code)] // Absolute is never constructed
enum ErrorType {
    Absolute,
    Relative,
}

#[derive(Debug)]
enum SortOrder {
    IncreasingMaxAbsErr,
    IncreasingAbsErrOnChannel(u8),
}

#[derive(Debug)]
struct Solution {
    clkfbout_mult_f_x8: u16,
    divclk_divide: u8,
    clkout_divide: Vec<u16>,

    vco_megahz: f64,
    out_megahz: Vec<(f64, f64, bool)>, // abs, err, print_red
    rmse: f64,
    max_abs_err_megahz: f64,
    max_abs_err_ppm: f64,
}

#[derive(Debug)]
struct SolutionSet {
    error_type: ErrorType,
    sort_order: SortOrder,
    solutions: Vec<Solution>,
}

// TODO: Understand lifetimes enough to know how to clean this up.
fn find_close_vco_freq<'a>(
    vco: f64,
    tolerance: f64,
    vec: &'a [(f64, u16, u8)],
) -> Option<&'a (f64, u16, u8)> {
    vec.iter()
        .find(|&x| ((x.0 / vco) > (1_f64 - tolerance)) && ((x.0 / vco) < (1_f64 + tolerance)))
}

impl From<Requirements> for SolutionSet {
    // TODO: split this up into something more digestable
    fn from(reqs: Requirements) -> Self {
        let mut set = SolutionSet {
            error_type: ErrorType::Relative,
            sort_order: SortOrder::IncreasingMaxAbsErr,
            solutions: Vec::<Solution>::new(),
        };

        let in_num_min = cmp::max(
            16,
            round::ceil(reqs.vco_megahz_min * 1_f64 * 8_f64 / reqs.inp_megahz, 0) as u16,
        );
        let in_num_max = cmp::min(
            512,
            round::floor(reqs.vco_megahz_max * 106_f64 * 8_f64 / reqs.inp_megahz, 0) as u16,
        );
        //println!("in_num_min {}, in_num_max {}", in_num_min, in_num_max);

        for (idx, constraint) in reqs.out_megahz.iter().enumerate() {
            if let OutputConstraint::Equal(_) = constraint {
                set.sort_order = SortOrder::IncreasingAbsErrOnChannel(idx as u8);
            }
        }

        let mut vco_solns = Vec::<(f64, u16, u8)>::new();
        for in_num in in_num_min..=in_num_max {
            let in_den_min = cmp::max(
                1,
                round::ceil(
                    (reqs.inp_megahz * (in_num as f64) / 8_f64) / reqs.vco_megahz_max,
                    0,
                ) as u8,
            );
            let in_den_max = cmp::min(
                106,
                round::floor(
                    (reqs.inp_megahz * (in_num as f64) / 8_f64) / reqs.vco_megahz_min,
                    0,
                ) as u8,
            );
            //println!("in_num {}, in_den_min {}, in_den_max {}", in_num, in_den_min, in_den_max);

            for in_den in in_den_min..=in_den_max {
                let vco = (reqs.inp_megahz * 0.125 * (in_num as f64)) / (in_den as f64);

                // check for out of bounds vco freq
                if vco < reqs.vco_megahz_min {
                    println!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco lo break\n",
                        reqs.inp_megahz,
                        0.125 * (in_num as f64),
                        in_den,
                        vco
                    );
                    break;
                }
                if vco > reqs.vco_megahz_max {
                    println!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco hi continue\n",
                        reqs.inp_megahz,
                        0.125 * (in_num as f64),
                        in_den,
                        vco
                    );
                    continue;
                }

                match find_close_vco_freq(vco, 1e-9, &vco_solns) {
                    Some(_) => (),
                    None => {
                        vco_solns.push((vco, in_num, in_den));
                    }
                }
            }
        }

        // sort from high to low vco frequencies to reduce output jitter
        vco_solns.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

        for (vco, in_num, in_den) in &vco_solns {
            //println!("vco {:4.6}, in_num {}, in_den {}, ", vco, in_num, in_den);

            // compute nearest output frequencies
            let mut out_div = vec![0_u16; reqs.out_megahz.len()];
            let mut out_megahz = vec![0_f64; reqs.out_megahz.len()];
            let mut out_err = vec![0_f64; reqs.out_megahz.len()];
            let mut sse = 0_f64;

            let mut max_abs_err_megahz = 0_f64;
            let mut _max_abs_err_chan = -1_i16;
            let mut max_err_ppm = 0_f64;
            let mut max_err_ppm_chan = -1_i16;
            for (chan, constraint) in reqs.out_megahz.iter().enumerate() {
                let target_mean = match constraint {
                    OutputConstraint::Normal(target) => *target,
                    OutputConstraint::Range {
                        min: target_min,
                        max: target_max,
                    } => (*target_min + *target_max) / 2_f64,
                    OutputConstraint::LessThan(target) => *target,
                    OutputConstraint::LessThanOrEqual(target) => *target,
                    OutputConstraint::Equal(target) => *target,
                    OutputConstraint::GreaterThanOrEqual(target) => *target,
                    OutputConstraint::GreaterThan(target) => *target,
                };

                out_div[chan] = match chan {
                    0 => (8_f64 * vco / target_mean).round() as u16,
                    _ => 8 * (vco / target_mean).round() as u16,
                };

                if out_div[chan] < 7 {                        
                    println!("ERROR: out_div[{}]: {} - too low continue\n", chan, out_div[chan]);
                    continue;
                } else if out_div[chan] < 8 {
                    out_div[chan] = 8;
                }

                if out_div[chan] > 1025 {
                    println!("ERROR: out_div[{}]: {} - too high continue\n", chan, out_div[chan]);
                    continue;
                } else if out_div[chan] > 1024 {
                    out_div[chan] = 1024;
                }

                out_megahz[chan] = (8_f64 * vco) / (out_div[chan] as f64);

                out_err[chan] = out_megahz[chan] - target_mean;
                let abs_err_megahz = out_err[chan].abs();
                let abs_err_ppm = (out_err[chan] / target_mean).abs();

                sse += out_err[chan] * out_err[chan];

                if abs_err_megahz > max_abs_err_megahz {
                    max_abs_err_megahz = abs_err_megahz;
                    _max_abs_err_chan = chan as i16;
                }
                if abs_err_ppm > max_err_ppm {
                    max_err_ppm = abs_err_ppm;
                    max_err_ppm_chan = chan as i16;
                }

                // check output range constraints
                match constraint {
                    OutputConstraint::LessThan(_) => {
                        if out_err[chan] >= 0_f64 {
                            continue;
                        }
                    }
                    OutputConstraint::LessThanOrEqual(_) => {
                        if out_err[chan] > 0_f64 {
                            continue;
                        }
                    }
                    OutputConstraint::GreaterThanOrEqual(_) => {
                        if out_err[chan] < 0_f64 {
                            continue;
                        }
                    }
                    OutputConstraint::GreaterThan(_) => {
                        if out_err[chan] <= 0_f64 {
                            continue;
                        }
                    }
                    OutputConstraint::Range {
                        min: target_min,
                        max: target_max,
                    } => {
                        if (out_megahz[chan] < *target_min) || (*target_max < out_megahz[chan]) {
                            continue;
                        }
                    }
                    _ => (),
                }
            }

            //---- Solution is valid ----

            let mut soln = Solution {
                clkfbout_mult_f_x8: *in_num,
                divclk_divide: *in_den,
                vco_megahz: *vco,

                clkout_divide: out_div,
                out_megahz: Vec::<(f64, f64, bool)>::new(), // abs, err, print_red
                rmse: sse.sqrt(),
                max_abs_err_megahz,
                max_abs_err_ppm: max_err_ppm,
            };
            for chan in 0..reqs.out_megahz.len() {
                let print_red = chan == (max_err_ppm_chan as usize);
                soln.out_megahz
                    .push((out_megahz[chan], out_err[chan], print_red));
            }
            set.solutions.push(soln);
        }

        //---- Sort ----

        match set.sort_order {
            SortOrder::IncreasingMaxAbsErr => {
                set.solutions
                    .sort_by(|a, b| a.max_abs_err_ppm.partial_cmp(&b.max_abs_err_ppm).unwrap());
            }
            SortOrder::IncreasingAbsErrOnChannel(ch) => {
                set.solutions.sort_by(|a, b| {
                    a.out_megahz[ch as usize]
                        .1
                        .abs()
                        .partial_cmp(&b.out_megahz[ch as usize].1.abs())
                        .unwrap()
                });
            }
        };

        //---- Trim to the requested number of solutions ----

        set.solutions.truncate(reqs.max_solutions);
        set
    }
}

impl fmt::Display for SolutionSet {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let num_outputs = self.solutions[0].out_megahz.len();

        //-- Table Annotation
        match self.sort_order {
            SortOrder::IncreasingMaxAbsErr => {
                writeln!(f, "Sorting in order of increasing ppm_err_max")
            }
            SortOrder::IncreasingAbsErrOnChannel(ch) => {
                writeln!(f, "Sorting in order of increasing error on channel {}", ch)
            }
        }?;
        match self.error_type {
            ErrorType::Absolute => writeln!(f, "Worst absolute output in {}", "red".red()),
            ErrorType::Relative => writeln!(f, "Worst ratiomnetric output in {}", "red".red()),
        }?;
        writeln!(f)?;

        //-- Header
        write!(f, "{:>5} {:>12}", "sln#", "vco")?;
        for ii in 0..num_outputs {
            write!(f, " {:>12}{:1}", "MHz", ii)?;
        }
        write!(f, " {:>13}", "ppm_err_max")?;
        write!(f, " {:>13}", "rms_err(MHz)")?;
        write!(f, " {:>6}", "clkfb")?;
        write!(f, " {:>6}", "divclk")?;
        for ii in 0..num_outputs {
            write!(f, " {:>5}{:1}", "odiv", ii)?;
        }

        //-- Solutions
        for (ii, soln) in self.solutions.iter().enumerate() {
            writeln!(f)?;
            write!(f, "{:>5} {:>12.6}", ii, soln.vco_megahz)?;
            for (megahz, _, is_red) in &soln.out_megahz {
                let mut str_megahz: String = format!(" {:>13.6}", megahz);
                if *is_red {
                    str_megahz = str_megahz.red().to_string();
                }
                write!(f, "{}", str_megahz)?;
            }
            write!(f, " {:>13.3}", 1e6 * soln.max_abs_err_ppm)?;
            write!(f, " {:>13.6}", soln.rmse)?;
            write!(f, " {:>6.3}", soln.clkfbout_mult_f_x8 as f64 / 8_f64)?;
            write!(f, " {:>6}", soln.divclk_divide)?;
            for (ii, val) in soln.clkout_divide.iter().enumerate() {
                match ii {
                    0 => write!(f, " {:6.3}", *val as f64 / 8_f64),
                    _ => write!(f, " {:6}", *val / 8),
                }?;
            }
        }

        write!(f, "")
    }
}

fn main() {
    println!("{}", SolutionSet::from(Requirements::from(Opt::from_args())));
}

// TODO: add a test or 20
