//! Value parsers enforcing the argument ranges of the original CLI.

fn parse_f64(text: &str) -> Result<f64, String> {
    text.parse::<f64>()
        .map_err(|_| format!("'{text}' is not a number"))
}

/// A fraction in `[0, 1]`.
pub fn unit_interval(text: &str) -> Result<f64, String> {
    let value = parse_f64(text)?;
    if (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(format!("{value} is not in the range 0..1"))
    }
}

/// A fraction in `(0, 1]`.
pub fn positive_unit_interval(text: &str) -> Result<f64, String> {
    let value = unit_interval(text)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err("value must be greater than 0".to_string())
    }
}

/// A number `>= 0`.
pub fn non_negative(text: &str) -> Result<f64, String> {
    let value = parse_f64(text)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(format!("{value} is negative"))
    }
}

/// An integer `>= 1`.
pub fn at_least_one(text: &str) -> Result<i64, String> {
    match text.parse::<i64>() {
        Ok(value) if value >= 1 => Ok(value),
        Ok(value) => Err(format!("{value} is smaller than 1")),
        Err(_) => Err(format!("'{text}' is not an integer")),
    }
}

/// An integer `>= 1`, as a size.
pub fn positive_size(text: &str) -> Result<usize, String> {
    at_least_one(text).map(|v| v as usize)
}
