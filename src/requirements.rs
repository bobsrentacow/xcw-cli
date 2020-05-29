use std::convert::TryFrom;
use std::fmt;
use crate::cli_args::Opt;

//----
// OutputConstraint
//   constrain the solution derived for each output

pub enum OutputConstraint {
    Normal(f64),
    RangeInclusive { min: f64, max: f64 },
    LessThan(f64),
    LessThanOrEqual(f64),
    // Use equal to sort on a single channel of output error
    Equal(f64),
    GreaterThanOrEqual(f64),
    GreaterThan(f64),
}

impl OutputConstraint {
    fn try_parse(the_string: &str) -> Result<Self, &'static str> {
        // You would be better off using something like
        // https://docs.rs/scan-rules/0.2.0/scan_rules/#rule-syntax here.

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
                Ok(OutputConstraint::RangeInclusive {
                    min: target * (1_f64 - 1e-6 * tolerance),
                    max: target * (1_f64 + 1e-6 * tolerance),
                })
            };
        let try_scan_plus_minus_percent =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                let tolerance: f64;
                text_io::try_scan!(string.bytes() => "{}+-{}pct", target, tolerance);
                Ok(OutputConstraint::RangeInclusive {
                    min: target * (1_f64 - 1e-2 * tolerance),
                    max: target * (1_f64 + 1e-2 * tolerance),
                })
            };
        let try_scan_plus_minus =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                let tolerance: f64;
                text_io::try_scan!(string.bytes() => "{}+-{}", target, tolerance);
                Ok(OutputConstraint::RangeInclusive {
                    min: target - tolerance,
                    max: target + tolerance,
                })
            };
        let try_scan_range =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let min: f64;
                let max: f64;
                text_io::try_scan!(string.bytes() => "{}-{}", min, max);
                Ok(OutputConstraint::RangeInclusive { min, max })
            };
        let try_scan_normal =
            |string: &str| -> Result<OutputConstraint, Box<dyn std::error::Error>> {
                let target: f64;
                text_io::try_scan!(string.bytes() => "{}", target);
                Ok(OutputConstraint::Normal(target))
            };

        // TODO: is there a more idiomatic way to do this in rust?
        // Nope, this seems about as idiotic as it gets ;) You crushed it.
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
//   This should fully specify the problem we are trying to solve.

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Fraction {
    pub num: u16,
    pub den: u16,
}

impl Into<f64> for Fraction {
    fn into(self) -> f64 {
        (self.num as f64) / (self.den as f64)
    }
}

#[derive(Debug, PartialEq)]
#[allow(dead_code)] // Absolute is never constructed
pub enum ErrorType {
    Absolute,
    Ratiometric,
}
impl fmt::Display for ErrorType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorType::Absolute => write!(f, "Absolute"),
            ErrorType::Ratiometric => write!(f, "Ratiometric"),
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum SortOrder {
    RootMeanSquareError,
    RatiometricErrorWorstChannel,
    RatiometricErrorOnChannel(u8),
}

impl fmt::Display for SortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SortOrder::RootMeanSquareError => write!(f, "RootMeanSquareError"),
            SortOrder::RatiometricErrorWorstChannel => write!(f, "RatiometricErrorWorstChannel"),
            SortOrder::RatiometricErrorOnChannel(ch) => {
                write!(f, "RatiometricErrorOnChannel({})", ch)
            }
        }
    }
}

pub struct Requirements {
    pub error_type: ErrorType,
    pub sort_order: SortOrder,
    pub max_solutions: usize,
    pub inp_megahz: f64,
    pub output_constraints: Vec<OutputConstraint>,
    // target hardware description
    pub vco_divider_num_min: Fraction,
    pub vco_divider_num_max: Fraction,
    pub vco_divider_den_min: u16,
    pub vco_divider_den_max: u16,
    pub vco_megahz_max: f64,
    pub vco_megahz_min: f64,
    pub chan_divider_min: Vec<Fraction>,
    pub chan_divider_max: Vec<Fraction>,
}

impl TryFrom<Opt> for Requirements {
    type Error = &'static str;

    fn try_from(opt: Opt) -> Result<Self, Self::Error> {
        if opt.sort_by_rmse && opt.sort_by_worst {
            return Err("Can't specify two different sort orders");
        }
        if opt.use_mmcm && opt.use_pll {
            return Err("must specify exactly one target");
        }

        let error_type = ErrorType::Ratiometric;
        let mut sort_order = {
            if opt.sort_by_rmse {
                SortOrder::RootMeanSquareError
            } else if opt.sort_by_worst {
                SortOrder::RatiometricErrorWorstChannel
            } else {
                SortOrder::RootMeanSquareError
            }
        };

        let max_outputs;
        let vco_divider_num_min;
        let vco_divider_num_max;
        let vco_divider_den_min;
        let vco_divider_den_max;
        let vco_megahz_max;
        let vco_megahz_min;
        let mut chan_divider_min = Vec::<Fraction>::new();
        let mut chan_divider_max = Vec::<Fraction>::new();

        // What hardware are we targetting?  Apply associated limits.
        if opt.use_mmcm {
            max_outputs = 8;
            vco_divider_num_min = Fraction { num: 2 * 8, den: 8 };
            vco_divider_num_max = Fraction {
                num: 64 * 8,
                den: 8,
            };
            vco_divider_den_min = 1;
            vco_divider_den_max = 106;
            vco_megahz_max = 1200_f64;
            vco_megahz_min = 600_f64;
            chan_divider_min.push(Fraction { num: 1 * 8, den: 8 });
            chan_divider_max.push(Fraction {
                num: 128 * 8,
                den: 8,
            });
            for _ in 1..max_outputs {
                chan_divider_min.push(Fraction { num: 1, den: 1 });
                chan_divider_max.push(Fraction { num: 128, den: 1 });
            }
        } else if opt.use_pll {
            max_outputs = 2;
            vco_divider_num_min = Fraction { num: 2 * 8, den: 8 };
            vco_divider_num_max = Fraction {
                num: 64 * 8,
                den: 8,
            };
            vco_divider_den_min = 1;
            vco_divider_den_max = 106;
            vco_megahz_max = 1200_f64;
            vco_megahz_min = 600_f64;
            for _ in 0..max_outputs {
                chan_divider_min.push(Fraction { num: 1, den: 1 });
                chan_divider_max.push(Fraction { num: 128, den: 1 });
            }
        } else {
            return Err("must specify exactly one target");
        }

        // Check limits associated with target.
        if opt.output_specifiers.len() > max_outputs {
            return Err("too many outputs requested");
        }

        // Figure out what constraints to apply for each requested output.
        let mut output_constraints = Vec::<OutputConstraint>::new();
        for (ii, the_string) in opt.output_specifiers.iter().enumerate() {
            match OutputConstraint::try_parse(&the_string) {
                Ok(constraint) => {
                    if let OutputConstraint::Equal(_) = constraint {
                        if opt.sort_by_rmse || opt.sort_by_worst {
                            return Err(
                                "Can't constrain an output to eq when a sort order is specified",
                            );
                        }
                        sort_order = SortOrder::RatiometricErrorOnChannel(ii as u8);
                    }
                    output_constraints.push(constraint)
                }
                Err(err_str) => return Err(&err_str),
            }
        }

        Ok(Requirements {
            error_type,
            sort_order,
            max_solutions: opt.max_solutions,
            inp_megahz: opt.inp_megahz,
            output_constraints,

            vco_divider_num_min,
            vco_divider_num_max,
            vco_divider_den_min,
            vco_divider_den_max,
            vco_megahz_max,
            vco_megahz_min,
            chan_divider_min,
            chan_divider_max,
        })
    }
}

impl fmt::Display for Requirements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "error_type         : {}", self.error_type)?;
        writeln!(f, "sort_order         : {}", self.sort_order)?;
        writeln!(f, "max_solutions      : {}", self.max_solutions)?;
        writeln!(
            f,
            "vco_divider_num_min: {}",
            Into::<f64>::into(self.vco_divider_num_min)
        )?;
        writeln!(
            f,
            "vco_divider_num_max: {}",
            Into::<f64>::into(self.vco_divider_num_max)
        )?;
        writeln!(f, "vco_divider_den_min: {}", self.vco_divider_den_min)?;
        writeln!(f, "vco_divider_den_max: {}", self.vco_divider_den_max)?;
        writeln!(f, "vco_megahz_max     : {}", self.vco_megahz_max)?;
        writeln!(f, "vco_megahz_min     : {}", self.vco_megahz_min)?;
        writeln!(f, "inp_megahz         : {}", self.inp_megahz)?;
        for (ii, constraint) in self.output_constraints.iter().enumerate() {
            match constraint {
                OutputConstraint::Normal(target) => writeln!(f, "output{:02}:   {}", ii, *target),
                OutputConstraint::RangeInclusive { min, max } => {
                    writeln!(f, "output{:02}:   {}-{}", ii, *min, *max)
                }
                OutputConstraint::LessThan(target) => writeln!(f, "output{:02}: < {}", ii, *target),
                OutputConstraint::LessThanOrEqual(target) => {
                    writeln!(f, "output{:02}: <={}", ii, *target)
                }
                OutputConstraint::Equal(target) => writeln!(f, "output{:02}:  ={}", ii, *target),
                OutputConstraint::GreaterThanOrEqual(target) => {
                    writeln!(f, "output{:02}: >={}", ii, *target)
                }
                OutputConstraint::GreaterThan(target) => {
                    writeln!(f, "output{:02}: > {}", ii, *target)
                }
            }?;
        }
        Ok(())
    }
}

//----
// Test

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_requirements_from_opt_no_target_specified() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: false,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gte187.5")],
        };

        if let Ok(_) = Requirements::try_from(opt) {
            Err(String::from("Requirements::try_from(opt) Should have failed"))
        } else {
            Ok(())
        }
    }

    #[test]
    fn test_requirements_from_opt_use_mmcm() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: false,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gte187.5")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 8);
            assert_eq!(reqs.chan_divider_max.len(), 8);
            Ok(())
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_use_pll() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gte187.5")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);
            Ok(())
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_use_mmcm_and_pll() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: true,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gte187.5")],
        };

        if let Ok(_) = Requirements::try_from(opt) {
            Err(String::from("Requirements::try_from(opt) Should have failed"))
        } else {
            Ok(())
        }
    }

    #[test]
    fn test_requirements_from_opt_normal() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::Normal(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_range() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("181.1-188.8")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::RangeInclusive { min, max } = reqs.output_constraints[1] {
                assert!((min - 181.1).abs() < 1e-6);
                assert!((max - 188.8).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_range_pct() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("181.1+-0.9pct")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::RangeInclusive { min, max } = reqs.output_constraints[1] {
                assert!((min - (181.1 * 0.991)).abs() < 1e-6);
                assert!((max - (181.1 * 1.009)).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_range_ppm() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("181.1+-3.5ppm")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::RangeInclusive { min, max } = reqs.output_constraints[1] {
                assert!((min - (181.1 * (1_f64 - 3.5e-6))).abs() < 1e-6);
                assert!((max - (181.1 * (1_f64 + 3.5e-6))).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_less_than() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("lt181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::LessThan(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_less_than_or_equal() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("lte181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::LessThanOrEqual(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_equal() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("eq181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RatiometricErrorOnChannel(1));
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::Equal(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_greater_than_or_equal() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gte181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::GreaterThanOrEqual(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

    #[test]
    fn test_requirements_from_opt_greater_than() -> Result<(), String> {
        let opt = Opt {
            use_mmcm: false,
            use_pll: true,
            sort_by_rmse: false,
            sort_by_worst: false,
            max_solutions: 32,
            inp_megahz: 156.25,
            output_specifiers: vec![String::from("166.6"), String::from("gt181.1")],
        };
        let max_solutions = opt.max_solutions;
        let inp_megahz = opt.inp_megahz;
        let num_outputs = opt.output_specifiers.len();

        if let Ok(reqs) = Requirements::try_from(opt) {
            assert_eq!(reqs.error_type, ErrorType::Ratiometric);
            assert_eq!(reqs.sort_order, SortOrder::RootMeanSquareError);
            assert_eq!(reqs.max_solutions, max_solutions);
            assert_eq!(reqs.inp_megahz, inp_megahz);
            assert_eq!(reqs.output_constraints.len(), num_outputs);
            // target hardware description
            assert_eq!(reqs.vco_divider_num_min, Fraction{ num: 2 * 8, den: 8});
            assert_eq!(reqs.vco_divider_num_max, Fraction{ num: 64 * 8, den: 8});
            assert_eq!(reqs.vco_divider_den_min, 1);
            assert_eq!(reqs.vco_divider_den_max, 106);
            assert_eq!(reqs.vco_megahz_max, 1200_f64);
            assert_eq!(reqs.vco_megahz_min, 600_f64);
            assert_eq!(reqs.chan_divider_min.len(), 2);
            assert_eq!(reqs.chan_divider_max.len(), 2);

            if let OutputConstraint::GreaterThan(target) = reqs.output_constraints[1] {
                assert!((target- 181.1).abs() < 1e-6);
                Ok(())
            } else {
                Err(String::from("parsing error"))
            }
        } else {
            Err(String::from("failed to parse requirements"))
        }
    }

}
