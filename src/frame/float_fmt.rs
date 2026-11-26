//! Julia-compatible `Float64` text formatting.
//!
//! Julia prints floats with the shortest digit string that round-trips
//! (Ryu) and switches to exponent notation outside `1e-4 <= |x| < 1e6`.
//! Reproducing this exactly is what keeps numeric TSV columns byte-identical
//! with the original implementation.

/// Format `x` exactly as Julia's `string(x::Float64)` does.
pub fn format_float(x: f64) -> String {
    if x.is_nan() {
        return "NaN".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    let sign = if x.is_sign_negative() { "-" } else { "" };
    if x == 0.0 {
        return format!("{sign}0.0");
    }

    // `{:e}` yields the shortest round-trip mantissa, e.g. "1.2345e-5".
    let scientific = format!("{:e}", x.abs());
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("LowerExp output always contains an exponent");
    let exponent: i32 = exponent.parse().expect("exponent is an integer");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let n_digits = digits.len() as i32;
    // Number of digits in front of the decimal point.
    let point = exponent + 1;

    let body = if -4 < point && point <= 6 {
        if point <= 0 {
            format!("0.{}{}", "0".repeat((-point) as usize), digits)
        } else if point >= n_digits {
            format!("{}{}.0", digits, "0".repeat((point - n_digits) as usize))
        } else {
            let (int_part, frac_part) = digits.split_at(point as usize);
            format!("{int_part}.{frac_part}")
        }
    } else {
        let (head, tail) = digits.split_at(1);
        let tail = if tail.is_empty() { "0" } else { tail };
        format!("{head}.{tail}e{exponent}")
    };
    format!("{sign}{body}")
}

/// Julia's `round(x, digits=d)` (round-half-even on the scaled value).
pub fn round_digits(x: f64, digits: i32) -> f64 {
    let scale = 10f64.powi(digits);
    let rounded = (x * scale).round_ties_even() / scale;
    if rounded.is_finite() {
        rounded
    } else {
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_julia_reference_output() {
        // Reference strings produced by Julia 1.13 `string(x)`.
        let cases: &[(f64, &str)] = &[
            (1e-5, "1.0e-5"),
            (1e-4, "0.0001"),
            (0.00012345, "0.00012345"),
            (123456.7, "123456.7"),
            (999999.0, "999999.0"),
            (1e6, "1.0e6"),
            (1.5e6, "1.5e6"),
            (1e21, "1.0e21"),
            (0.1 + 0.2, "0.30000000000000004"),
            (1.0 / 3.0, "0.3333333333333333"),
            (f64::NAN, "NaN"),
            (f64::INFINITY, "Inf"),
            (f64::NEG_INFINITY, "-Inf"),
            (-0.0, "-0.0"),
            (100.0, "100.0"),
            (1e-10, "1.0e-10"),
            (5.0e-324, "5.0e-324"),
            (1.0e15, "1.0e15"),
            (123456789012345680.0, "1.2345678901234568e17"),
            (0.5, "0.5"),
            (2.0, "2.0"),
            (1e5, "100000.0"),
            (12345678.9, "1.23456789e7"),
            (0.001, "0.001"),
            (-2.5, "-2.5"),
        ];
        for (value, expected) in cases {
            assert_eq!(&format_float(*value), expected, "value {value:?}");
        }
    }

    #[test]
    fn rounding_is_half_even_on_scaled_value() {
        assert_eq!(round_digits(0.8, 3), 0.8);
        assert_eq!(round_digits(2.0 / 3.0, 3), 0.667);
        assert_eq!(round_digits(0.0005, 3), 0.0); // 0.5 -> even
        assert!(round_digits(f64::NAN, 1).is_nan());
    }
}
