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

// TODO: Couldn't figure out how to get Into<f64> working.
// Yeah it's a little tricky because you can't call `into()` in the print! macro as it takes in
// lots of types. So you have to explicitly call it with Into::<T>::into(...)
#[derive(Debug, Copy, Clone)]
pub struct Fraction {
    pub num: u16,
    pub den: u16,
}

impl Into<f64> for Fraction {
    fn into(self) -> f64 {
        (self.num as f64) / (self.den as f64)
    }
}

#[derive(Debug)]
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

#[derive(Debug)]
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

    // Try to refactor this to use no `mut` variables, or at most a few.
    fn try_from(opt: Opt) -> Result<Self, Self::Error> {
        // This isn't C. Also mutable variables are almost always a mistake, they should be needed
        // very, very rarely (see sort_order example below).
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

        // Don't forget that Rust is expression oriented, so you can do cool stuff like this.
        let mut sort_order = {
            if opt.sort_by_rmse {
                SortOrder::RootMeanSquareError
            } else if opt.sort_by_worst {
                SortOrder::RatiometricErrorWorstChannel
            } else {
                SortOrder::RootMeanSquareError
            }
        };

        // What hardware are we targetting?  Apply associated limits.
        if opt.use_mmcm && opt.use_pll {
            return Err("must specify exactly one target");
        }

        // You can use the same trick to turn this entire block into an expression that creates the
        // correct Requirements instance.
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
                chan_divider_min.push(Fraction { num: 128, den: 1 });
                chan_divider_max.push(Fraction { num: 128, den: 1 });
            }
        } else {
            return Err("must specify exactly one target");
        }

        // Check limits associated with target.
        if opt.output_specifiers.len() > max_outputs {
            return Err("too many outputs requested");
        }

        let mut output_constraints = Vec::<OutputConstraint>::new();

        // TODO: reduce number of indentations
        // Figure out what constraints to apply for each requested output.
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

