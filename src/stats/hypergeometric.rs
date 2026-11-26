//! Hypergeometric upper-tail probabilities.

/// A double-double number (`hi + lo`), giving about 31 significant digits.
///
/// The tail probability is a ratio of two sums of rapidly varying terms;
/// carrying it out in extended precision makes the final `f64` correctly
/// rounded in practice.
#[derive(Debug, Clone, Copy)]
struct DoubleDouble {
    hi: f64,
    lo: f64,
}

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn quick_two_sum(a: f64, b: f64) -> DoubleDouble {
    let s = a + b;
    DoubleDouble {
        hi: s,
        lo: b - (s - a),
    }
}

impl DoubleDouble {
    const ZERO: DoubleDouble = DoubleDouble { hi: 0.0, lo: 0.0 };
    const ONE: DoubleDouble = DoubleDouble { hi: 1.0, lo: 0.0 };

    fn add(self, other: DoubleDouble) -> DoubleDouble {
        let (s, e) = two_sum(self.hi, other.hi);
        quick_two_sum(s, e + self.lo + other.lo)
    }

    fn mul_f64(self, b: f64) -> DoubleDouble {
        let p = self.hi * b;
        let e = self.hi.mul_add(b, -p) + self.lo * b;
        quick_two_sum(p, e)
    }

    fn sub_scaled(self, other: DoubleDouble, factor: f64) -> DoubleDouble {
        let scaled = other.mul_f64(factor);
        self.add(DoubleDouble {
            hi: -scaled.hi,
            lo: -scaled.lo,
        })
    }

    fn div_f64(self, b: f64) -> DoubleDouble {
        self.div(DoubleDouble { hi: b, lo: 0.0 })
    }

    fn div(self, other: DoubleDouble) -> DoubleDouble {
        let q1 = self.hi / other.hi;
        let r = self.sub_scaled(other, q1);
        let q2 = r.hi / other.hi;
        let r = r.sub_scaled(other, q2);
        let q3 = r.hi / other.hi;
        quick_two_sum(q1, q2).add(DoubleDouble { hi: q3, lo: 0.0 })
    }

    fn is_zero(self) -> bool {
        self.hi == 0.0
    }

    fn to_f64(self) -> f64 {
        self.hi + self.lo
    }
}

/// `P(X >= k)` for `X ~ Hypergeometric(successes, failures, draws)`:
/// the number of successes among `draws` items drawn without replacement
/// from a population of `successes + failures` items.
///
/// The tail is computed from pmf *ratios* anchored at the mode, so no
/// factorials or log-gammas are needed and relative accuracy stays at
/// machine precision even for very small tails.
pub fn upper_tail(successes: u64, failures: u64, draws: u64, k: i64) -> f64 {
    let population = successes + failures;
    if draws > population {
        return f64::NAN;
    }
    let lo = draws.saturating_sub(failures);
    let hi = successes.min(draws);
    if k <= lo as i64 {
        return 1.0;
    }
    if k > hi as i64 {
        return 0.0;
    }
    let k = k as u64;

    // pmf(x+1)/pmf(x) = (s-x)(d-x) / ((x+1)(f-d+x+1)); both products are
    // exact in f64 for any realistic population size.
    let numerator = |x: u64| (successes - x) as f64 * (draws - x) as f64;
    let denominator = |x: u64| (x + 1) as f64 * (failures + x + 1 - draws) as f64;
    let mode =
        (((draws + 1) as f64 * (successes + 1) as f64) / (population + 2) as f64).floor() as u64;
    let mode = mode.clamp(lo, hi);

    // Unnormalised pmf relative to the mode.
    let mut total = DoubleDouble::ONE;
    let mut tail = if mode >= k {
        DoubleDouble::ONE
    } else {
        DoubleDouble::ZERO
    };
    let mut term = DoubleDouble::ONE;
    for x in mode..hi {
        term = term.mul_f64(numerator(x)).div_f64(denominator(x));
        total = total.add(term);
        if x + 1 >= k {
            tail = tail.add(term);
        }
        if term.is_zero() {
            break;
        }
    }
    term = DoubleDouble::ONE;
    for x in (lo..mode).rev() {
        term = term.mul_f64(denominator(x)).div_f64(numerator(x));
        total = total.add(term);
        if x >= k {
            tail = tail.add(term);
        }
        if term.is_zero() {
            break;
        }
    }
    tail.div(total).to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values: Distributions.jl `ccdf(Hypergeometric(s, f, n), k - 1)`
    // (Rmath `phyper`). The reference is itself only accurate to roughly
    // twelve digits in the tails, hence the tolerance.
    #[test]
    fn agrees_with_distributions_jl() {
        let cases = [
            (8, 3, 4, 2, 0.9757575757575758),
            (8, 3, 4, 4, 0.2121212121212121),
            (50, 20, 30, 28, 0.0003298444166752913),
            (200, 150, 120, 100, 1.6230405434290412e-13),
        ];
        for (s, f, n, k, expected) in cases {
            let p = upper_tail(s, f, n, k);
            assert!(
                (p - expected).abs() <= 1e-11 * expected,
                "{s} {f} {n} {k}: {p} vs {expected}"
            );
        }
    }

    // Exact rational values, correctly rounded.
    #[test]
    fn simple_cases_are_correctly_rounded() {
        assert_eq!(upper_tail(8, 3, 4, 2), 322.0 / 330.0);
        assert_eq!(upper_tail(8, 3, 4, 4), 70.0 / 330.0);
        assert_eq!(upper_tail(2, 2, 2, 1), 5.0 / 6.0);
    }

    #[test]
    fn boundaries() {
        assert_eq!(upper_tail(8, 3, 4, 0), 1.0);
        assert_eq!(upper_tail(8, 3, 4, 1), 1.0); // 4 draws, 3 failures: one success is forced
        assert_eq!(upper_tail(8, 3, 4, 5), 0.0);
        assert!(upper_tail(2, 2, 9, 1).is_nan());
    }

    #[test]
    fn tail_is_monotone_and_complements_to_one() {
        let (s, f, n) = (30u64, 45u64, 20u64);
        let mut previous = 1.0;
        for k in 0..=21 {
            let p = upper_tail(s, f, n, k);
            assert!(p <= previous + 1e-15);
            previous = p;
        }
        assert_eq!(upper_tail(s, f, n, 21), 0.0);
    }
}
