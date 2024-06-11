//! This module contains code from <https://github.com/statrs-dev/statrs>
//!
//! The statrs crate contains way more functionality than needed for hpo
//! so it only contains the logic neccessary for the hypergeometric
//! enrichment
#![allow(clippy::excessive_precision)]
#![allow(clippy::unreadable_literal)]

use std::cmp;

use crate::stats::f64_from_usize;

/// Auxiliary variable when evaluating the `gamma_ln` function
const GAMMA_R: f64 = 10.900_511;

/// Polynomial coefficients for approximating the `gamma_ln` function
const GAMMA_DK: &[f64] = &[
    2.48574089138753565546e-5,
    1.05142378581721974210,
    -3.45687097222016235469,
    4.51227709466894823700,
    -2.98285225323576655721,
    1.05639711577126713077,
    -1.95428773191645869583e-1,
    1.70970543404441224307e-2,
    -5.71926117404305781283e-4,
    4.63399473359905636708e-6,
    -2.71994908488607703910e-9,
];
pub const LN_2_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;

/// Constant value for `ln(pi)`
pub const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;

/// The maximum factorial representable
/// by a 64-bit floating point without
/// overflowing
pub const MAX_FACTORIAL: usize = 170;

// Initialization for pre-computed cache of 171 factorial
// values 0!...170!
#[allow(clippy::cast_precision_loss)]
const FCACHE: [f64; MAX_FACTORIAL + 1] = {
    let mut fcache = [1.0; MAX_FACTORIAL + 1];

    let mut i = 1;
    while i < MAX_FACTORIAL + 1 {
        fcache[i] = fcache[i - 1] * i as f64;

        i += 1;
    }
    fcache
};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Hypergeometric {
    population: u64,
    successes: u64,
    draws: u64,
}

impl Hypergeometric {
    /// Constructs a new hypergeometric distribution
    /// with a population (N) of `population`, number
    /// of successes (K) of `successes`, and number of draws
    /// (n) of `draws`
    ///
    /// # Errors
    ///
    /// If `successes > population` or `draws > population`
    pub fn new(population: u64, successes: u64, draws: u64) -> Result<Hypergeometric, String> {
        if successes > population || draws > population {
            Err("Invalid params".to_string())
        } else {
            Ok(Hypergeometric {
                population,
                successes,
                draws,
            })
        }
    }

    /// Returns the minimum value in the domain of the
    /// hypergeometric distribution representable by a 64-bit
    /// integer
    ///
    /// # Formula
    ///
    /// ```text
    /// max(0, n + K - N)
    /// ```
    ///
    /// where `N` is population, `K` is successes, and `n` is draws
    fn min(&self) -> u64 {
        (self.draws + self.successes).saturating_sub(self.population)
    }

    /// Returns the maximum value in the domain of the
    /// hypergeometric distribution representable by a 64-bit
    /// integer
    ///
    /// # Formula
    ///
    /// ```text
    /// min(K, n)
    /// ```
    ///
    /// where `K` is successes and `n` is draws
    fn max(&self) -> u64 {
        cmp::min(self.successes, self.draws)
    }

    /// Calculates the survival function for the hypergeometric
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// 1 - ((n choose k+1) * (N-n choose K-k-1)) / (N choose K) * 3_F_2(1,
    /// k+1-K, k+1-n; k+2, N+k+2-K-n; 1)
    /// ```
    ///
    /// where `N` is population, `K` is successes, `n` is draws,
    /// and `p_F_q` is the [generalized hypergeometric function](https://en.wikipedia.org/wiki/Generalized_hypergeometric_function)
    ///
    /// Calculated as a discrete integral over the probability mass
    /// function evaluated from (k+1)..max
    pub fn sf(&self, x: u64) -> f64 {
        if x < self.min() {
            1.0
        } else if x >= self.max() {
            0.0
        } else {
            let k = x;
            let ln_denom = ln_binomial(self.population, self.draws);
            ((k + 1)..=self.max()).fold(0.0, |acc, i| {
                acc + (ln_binomial(self.successes, i)
                    + ln_binomial(self.population - self.successes, self.draws - i)
                    - ln_denom)
                    .exp()
            })
        }
    }
}

/// Computes the logarithm of the gamma function
/// with an accuracy of 16 floating point digits.
/// The implementation is derived from
/// "An Analysis of the Lanczos Gamma Approximation",
/// Glendon Ralph Pugh, 2004 p. 116
fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (f64_from_usize(t.0) - x));

        LN_PI
            - (std::f64::consts::PI * x).sin().ln()
            - s.ln()
            - LN_2_SQRT_E_OVER_PI
            - (0.5 - x) * ((0.5 - x + GAMMA_R) / std::f64::consts::E).ln()
    } else {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| {
                s + t.1 / (x + f64_from_usize(t.0) - 1.0)
            });

        s.ln() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / std::f64::consts::E).ln()
    }
}

/// Computes the logarithmic factorial function `x -> ln(x!)`
/// for `x >= 0`.
///
/// # Remarks
///
/// Returns `0.0` if `x <= 1`
fn ln_factorial(x: u64) -> f64 {
    let x = usize::try_from(x).expect("x must be castable to usize");
    FCACHE
        .get(x)
        .map_or_else(|| ln_gamma(f64_from_usize(x) + 1.0), |&fac| fac.ln())
}

/// Computes the natural logarithm of the binomial coefficient
/// `ln(n choose k)` where `k` and `n` are non-negative values
///
/// # Remarks
///
/// Returns `f64::NEG_INFINITY` if `k > n`
pub fn ln_binomial(n: u64, k: u64) -> f64 {
    if k > n {
        f64::NEG_INFINITY
    } else {
        ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_fcache() {
        assert!((FCACHE[0] - 1.0).abs() < f64::EPSILON);
        assert!((FCACHE[1] - 1.0).abs() < f64::EPSILON);
        assert!((FCACHE[2] - 2.0).abs() < f64::EPSILON);
        assert!((FCACHE[3] - 6.0).abs() < f64::EPSILON);
        assert!((FCACHE[4] - 24.0).abs() < f64::EPSILON);
        assert!(
            (
                FCACHE[70] -
                11978571669969890000000000000000000000000000000000000000000000000000000000000000000000000000000000000.0
            )
            .abs() < f64::EPSILON
        );
        assert!(
            (
                FCACHE[170] -
                7257415615307994000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.0
            )
            .abs() < f64::EPSILON
        );
    }

    #[test]
    fn test_hypergeom_build() {
        let mut result = Hypergeometric::new(2, 2, 2);
        assert!(result.is_ok());

        result = Hypergeometric::new(2, 3, 2);
        assert!(result.is_err());
    }

    #[test]
    fn test_hypergeom_max() {
        let hyper = Hypergeometric::new(50, 25, 13).unwrap();
        assert_eq!(hyper.max(), 13);

        let hyper = Hypergeometric::new(50, 10, 13).unwrap();
        assert_eq!(hyper.max(), 10);
    }

    #[test]
    fn min() {
        let hyper = Hypergeometric::new(50, 25, 30).unwrap();
        assert_eq!(hyper.min(), 5);

        let hyper = Hypergeometric::new(50, 40, 30).unwrap();
        assert_eq!(hyper.min(), 20);

        let hyper = Hypergeometric::new(50, 10, 13).unwrap();
        assert_eq!(hyper.min(), 0);
    }

    #[test]
    fn test_hypergeom_cdf() {
        // Numbers calculated here https://statisticsbyjim.com/probability/hypergeometric-distribution/
        let hyper = Hypergeometric::new(50, 25, 13).unwrap();

        // more than 1 == 2 or more
        assert!((hyper.sf(1) - 0.9996189832542451).abs() < f64::EPSILON);
        // more than 3 == 4 or more
        assert!((hyper.sf(3) - 0.9746644799047702).abs() < f64::EPSILON);
        // more than 7 == 8 or more
        assert!((hyper.sf(7) - 0.26009737477738537).abs() < f64::EPSILON);
        // more than 12 == 13 or more
        assert!((hyper.sf(12) - 0.000014654490222007184).abs() < f64::EPSILON);

        assert!(hyper.sf(13) < f64::EPSILON);
    }
}
