use colored::Colorize;
use math::round;
use std::cmp;
use std::fmt;
use crate::requirements::ErrorType;
use crate::requirements::Fraction;
use crate::requirements::OutputConstraint;
use crate::requirements::Requirements;
use crate::requirements::SortOrder;

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

        let virtual_vco_pre_mult = reqs.vco_megahz_min * vco_divider_num_den_f64 / reqs.inp_megahz;

        let min_num_in_vco_range =
            round::ceil(virtual_vco_pre_mult * (reqs.vco_divider_den_min as f64), 0) as u16;
        let max_num_in_vco_range =
            round::floor(virtual_vco_pre_mult * (reqs.vco_divider_den_max as f64), 0) as u16;

        let in_num_min = cmp::max(reqs.vco_divider_num_min.num, min_num_in_vco_range);
        let in_num_max = cmp::min(reqs.vco_divider_num_max.num, max_num_in_vco_range);
        log::trace!("in_num_min {}, in_num_max {}", in_num_min, in_num_max);

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
            log::trace!(
                "in_num {}, in_den_min {}, in_den_max {}",
                in_num,
                in_den_min,
                in_den_max
            );

            for in_den in in_den_min..=in_den_max {
                let vco = reqs.inp_megahz * vco_divider_num_f64 / (in_den as f64);

                // check for out of bounds vco freq
                if vco < reqs.vco_megahz_min {
                    log::error!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco lo break\n",
                        reqs.inp_megahz,
                        vco_divider_num_f64,
                        in_den,
                        vco
                    );
                    break;
                }
                if vco > reqs.vco_megahz_max {
                    log::error!(
                        "ERROR: {:11.6} * {:5.3}/{} = {:11.6} - vco hi continue\n",
                        reqs.inp_megahz,
                        vco_divider_num_f64,
                        in_den,
                        vco
                    );
                    continue;
                }

                let thresh = 1_f64 + 1e-9; // magic number for vco frequency equality threshold
                let found = vco_solns
                    .iter()
                    .find(|&x| ((x.output / vco) < thresh) && ((vco / x.output) < thresh));
                match found {
                    Some(_) => (),
                    None => {
                        vco_solns.push(VcoSolution {
                            input: reqs.inp_megahz,
                            output: vco,
                            numerator: Fraction {
                                num: in_num,
                                den: vco_divider_num_den,
                            },
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

impl PartialEq for VcoSolution {
    fn eq(&self, other: &Self) -> bool {
        if (self.input - other.input).abs() > 1e-6 {
            return false;
        }
        if (self.output - other.output).abs() > 1e-6 {
            return false;
        }
        if !self.numerator.eq(&other.numerator) {
            return false;
        }
        if self.denominator != other.denominator {
            return false;
        }
        return true;
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
        divider_max: Fraction,
    ) -> Result<ChannelSolution, String> {
        let den = divider_min.den;
        let den_f64 = den as f64;
        let target = match constraint {
            OutputConstraint::Normal(target) => *target,
            OutputConstraint::RangeInclusive { min, max } => (*min + *max) / 2_f64,
            OutputConstraint::LessThan(target) => *target,
            OutputConstraint::LessThanOrEqual(target) => *target,
            OutputConstraint::Equal(target) => *target,
            OutputConstraint::GreaterThanOrEqual(target) => *target,
            OutputConstraint::GreaterThan(target) => *target,
        };

        // get closest integer solution
        let num = (vco * den_f64 / target).round() as u16;
        // check num + -1..=+1 solutions
        let num_candidates = vec![num - 1, num, num + 1];

        // No. Just no. Make a struct lol.
        let mut num_tuples = Vec::<(u16, f64, f64, f64, f64)>::new();
        for num in num_candidates {
            if (num < divider_min.num) || (num > divider_max.num) {
                log::trace!(
                    "vco {}, dev {:.3} -- disqualified {} < {} || {} > {}",
                    vco,
                    Into::<f64>::into(Fraction { num, den }),
                    num,
                    divider_min.num,
                    num,
                    divider_max.num
                );
                continue;
            }

            let actual = vco * (den_f64 / (num as f64));
            let error = actual - target;
            let absolute_error = error.abs();
            let ratiometric_error = (error / target).abs();

            log::trace!(
                "vco {}, dev {:.3}, out {}",
                vco,
                Into::<f64>::into(Fraction { num, den }),
                actual
            );

            // check output range constraints
            match constraint {
                OutputConstraint::RangeInclusive { min, max } => {
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
                _ => (),
            }
            num_tuples.push((num, actual, error, absolute_error, ratiometric_error));
        }

        if !num_tuples.is_empty() {
            num_tuples.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap());
            let (num, actual, _error, absolute_error, ratiometric_error) = num_tuples[0];

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

impl PartialEq for ChannelSolution {
    fn eq(&self, other: &Self) -> bool {
        if (self.input - other.input) > 1e-6 {
            return false;
        }
        if self.chan_idx != other.chan_idx {
            return false;
        }
        if !self.divider.eq(&other.divider) {
            return false;
        }
        if (self.target - other.target) > 1e-6 {
            return false;
        }
        if (self.actual - other.actual) > 1e-6 {
            return false;
        }
        if (self.absolute_error - other.absolute_error) > 1e-6 {
            return false;
        }
        if (self.ratiometric_error - other.ratiometric_error) > 1e-6 {
            return false;
        }
        return true;
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

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:>12.6}", self.vco_solution.output)?;
        for (ii, chan_soln) in self.channel_solutions.iter().enumerate() {
            let mut str_megahz: String = format!(" {:>13.6}", chan_soln.actual);
            if ii == (self.channel_with_worst_error as usize) {
                str_megahz = str_megahz.red().to_string();
            }
            write!(f, "{}", str_megahz)?;
        }
        write!(f, " {:>13.3}", 1e6 * self.worst_error)?;
        write!(f, " {:>13.6}", self.root_mean_square_error)?;
        write!(
            f,
            " {:>6.3}",
            Into::<f64>::into(self.vco_solution.numerator)
        )?;
        write!(f, " {:>6}", self.vco_solution.denominator as f64)?;
        for chan_soln in &self.channel_solutions {
            write!(f, " {:6.3}", Into::<f64>::into(chan_soln.divider))?;
        }

        write!(f, "")
    }
}

fn solution_header(num_outputs: usize, error_type: &ErrorType) -> String {
    let mut res = format!("{:>12}", "vco");
    for ii in 0..num_outputs {
        res.push_str(&format!(" {:>12}{:1}", "MHz", ii));
    }
    match error_type {
        ErrorType::Absolute => res.push_str(&format!(" {:>13}", "MHz_err_max")),
        ErrorType::Ratiometric => res.push_str(&format!(" {:>13}", "ppm_err_max")),
    };
    res.push_str(&format!(" {:>13} {:>6} {:>6}", "rms_err(MHz)", "clkfb", "divclk"));
    for ii in 0..num_outputs {
        res.push_str(&format!(" {:>5}{:1}", "odiv", ii));
    }

    res
}

impl PartialEq for Solution {
    fn eq(&self, other: &Self) -> bool {
        if !self.vco_solution.eq(&other.vco_solution) {
            return false;
        }
        if !(self.channel_solutions.len() == other.channel_solutions.len()) {
            return false;
        }
        for (ii, chan_soln) in self.channel_solutions.iter().enumerate() {
            if !chan_soln.eq(&other.channel_solutions[ii]) {
                return false;
            }
        }
        if (self.root_mean_square_error - other.root_mean_square_error).abs() > 1e-6 {
            return false;
        }
        if (self.worst_error - other.worst_error).abs() > 1e-6 {
            return false;
        }
        if self.channel_with_worst_error != other.channel_with_worst_error {
            return false;
        }
        return true;
    }
}

#[derive(Debug)]
pub struct SolutionSet {
    error_type: ErrorType,
    sort_order: SortOrder,
    solutions: Vec<Solution>,
}

impl SolutionSet {
    pub fn from(reqs: Requirements) -> Self {
        // Generate candidate solutions
        let vco_solns = VcoSolution::get_solutions(&reqs);
        let mut solutions = Vec::<Solution>::new();

        // Already using labels, look at you, so fancy!
        'vco: for VcoSolution {
            input,
            output: vco_freq,
            numerator,
            denominator,
        } in &vco_solns
        {
            // Solve each output channel
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
                    reqs.chan_divider_max[chan],
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
                    }
                    Err(_) => continue 'vco,
                }
            }

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

        // Sort and trim
        match reqs.sort_order {
            SortOrder::RootMeanSquareError => {
                solutions.sort_by(|a, b| {
                    a.root_mean_square_error
                        .partial_cmp(&b.root_mean_square_error)
                        .unwrap()
                });
            }
            SortOrder::RatiometricErrorWorstChannel => {
                solutions.sort_by(|a, b| a.worst_error.partial_cmp(&b.worst_error).unwrap());
            }
            SortOrder::RatiometricErrorOnChannel(ch) => {
                solutions.sort_by(|a, b| {
                    a.channel_solutions[ch as usize]
                        .ratiometric_error
                        .partial_cmp(&b.channel_solutions[ch as usize].ratiometric_error)
                        .unwrap()
                });
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
                SortOrder::RootMeanSquareError => {
                    writeln!(f, "Sorting in order of increasing root_mean_square_error")
                }
                SortOrder::RatiometricErrorWorstChannel => {
                    writeln!(f, "Sorting in order of increasing ppm_err_max")
                }
                SortOrder::RatiometricErrorOnChannel(ch) => {
                    writeln!(f, "Sorting in order of increasing error on channel {}", ch)
                }
            }?;
            match self.error_type {
                ErrorType::Absolute => writeln!(f, "Worst absolute output in {}", "red".red()),
                ErrorType::Ratiometric => {
                    writeln!(f, "Worst ratiomnetric output in {}", "red".red())
                }
            }?;
            writeln!(f)?;

            //-- Header
            writeln!(f, "{:>5} {}", "sln#", solution_header(num_outputs, &self.error_type))?;

            //-- Solutions
            for (ii, soln) in self.solutions.iter().enumerate() {
                writeln!(f, "{:>5} {}", ii, soln)?;
            }

            write!(f, "")
        }
    }
}

//----
// Test

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryFrom;
    use crate::cli_args::Opt;

    #[test]
    fn test_mmcm_156p25_164p4_210p1() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: false,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("164.4"), String::from("210.1")],
        };
        let max_solutions = opt.max_solutions;

        match Requirements::try_from(opt) {
            Err(_) => Err(String::from("failed to parse requirements from opt")),
            Ok(reqs) => {
                let set = SolutionSet::from(reqs);
                assert_eq!(set.error_type, ErrorType::Ratiometric);
                assert_eq!(set.sort_order, SortOrder::RootMeanSquareError);
                if set.solutions.len() > max_solutions {
                    return Err(String::from("more solutions returned than requested"));
                }
                if !set.solutions[0].eq(&Solution {
                    vco_solution: VcoSolution {
                        input: 156.25,
                        output: 1049.8046875,
                        numerator: Fraction { num: 215, den: 8 },
                        denominator: 4,
                    },
                    channel_solutions: vec![
                        ChannelSolution {
                            input: 1049.8046875,
                            chan_idx: 0,
                            divider: Fraction { num: 51, den: 8 },
                            target: 164.4,
                            actual: 164.67524509803923,
                            absolute_error: 0.27524509803922115,
                            ratiometric_error: 0.001674240255713024,
                        },
                        ChannelSolution {
                            input: 1049.8046875,
                            chan_idx: 1,
                            divider: Fraction { num: 5, den: 1 },
                            target: 210.1,
                            actual: 209.9609375,
                            absolute_error: 0.13906249999999432,
                            ratiometric_error: 0.0006618871965730335,
                        },
                    ],
                    root_mean_square_error: 0.3083800299968675,
                    worst_error: 0.001674240255713024,
                    channel_with_worst_error: 0,
                }) {
                    println!("{}", solution_header(set.solutions[0].channel_solutions.len(), &set.error_type));
                    println!("{}", set.solutions[0]);
                    println!("{:?}", set.solutions[0].vco_solution);
                    println!("{:?}", set.solutions[0].channel_solutions);
                    println!("{:?}", set.solutions[0].root_mean_square_error);
                    println!("{:?}", set.solutions[0].worst_error);
                    println!("{:?}", set.solutions[0].channel_with_worst_error);
                    return Err(format!("Solution is different than when test was written"));
                }
                Ok(())
            }
        }
    }

    #[test]
    fn test_mmcm_156p25_164p4_210p1_to_210p9() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: false,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("164.4"), String::from("210.1-210.9")],
        };
        let max_solutions = opt.max_solutions;

        match Requirements::try_from(opt) {
            Err(_) => Err(String::from("failed to parse requirements from opt")),
            Ok(reqs) => {
                let set = SolutionSet::from(reqs);
                assert_eq!(set.error_type, ErrorType::Ratiometric);
                assert_eq!(set.sort_order, SortOrder::RootMeanSquareError);
                if set.solutions.len() > max_solutions {
                    return Err(String::from("more solutions returned than requested"));
                }
                if !set.solutions[0].eq(&Solution {
                    vco_solution: VcoSolution {
                        input: 156.25,
                        output: 842.28515625,
                        numerator: Fraction { num: 345, den: 8 },
                        denominator: 8,
                    },
                    channel_solutions: vec![
                        ChannelSolution {
                            input: 842.28515625,
                            chan_idx: 0,
                            divider: Fraction { num: 41, den: 8 },
                            target: 164.4,
                            actual: 164.3483231707317,
                            absolute_error: 0.05167682926830253,
                            ratiometric_error: 0.0003143359444543949,
                        },
                        ChannelSolution {
                            input: 842.28515625,
                            chan_idx: 1,
                            divider: Fraction { num: 4, den: 1 },
                            target: 210.5,
                            actual: 210.5712890625,
                            absolute_error: 0.0712890625,
                            ratiometric_error: 0.00033866538004750593,
                        },
                    ],
                    root_mean_square_error: 0.08804899269925917,
                    worst_error: 0.00033866538004750593,
                    channel_with_worst_error: 1,
                }) {
                    println!("{}", solution_header(set.solutions[0].channel_solutions.len(), &set.error_type));
                    println!("{}", set.solutions[0]);
                    println!("{:?}", set.solutions[0].vco_solution);
                    println!("{:?}", set.solutions[0].channel_solutions);
                    println!("{:?}", set.solutions[0].root_mean_square_error);
                    println!("{:?}", set.solutions[0].worst_error);
                    println!("{:?}", set.solutions[0].channel_with_worst_error);
                    return Err(format!("Solution is different than when test was written"));
                }
                Ok(())
            }
        }
    }

    // TODO: test other OutputConstraints
}
