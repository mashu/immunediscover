//! Cell values and their comparison/equality semantics.
//!
//! Ordering follows Julia's `isless`, equality follows `isequal`:
//! `missing` sorts last, `NaN` sorts after every other number and equals
//! itself, integers and floats compare numerically.

use std::cmp::Ordering;

use super::float_fmt::format_float;

/// A borrowed view of one table cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value<'a> {
    Missing,
    Int(i64),
    Float(f64),
    Bool(bool),
    Str(&'a str),
}

impl<'a> Value<'a> {
    pub fn is_missing(&self) -> bool {
        matches!(self, Value::Missing)
    }

    /// Numeric view (`Bool` counts as 0/1, like Julia's `Bool <: Integer`).
    pub fn as_f64(&self) -> Option<f64> {
        match *self {
            Value::Int(v) => Some(v as f64),
            Value::Float(v) => Some(v),
            Value::Bool(v) => Some(if v { 1.0 } else { 0.0 }),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&'a str> {
        match *self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    /// Text exactly as Julia's `string(x)` renders it.
    ///
    /// Note that `missing` becomes the literal `"missing"`; callers that
    /// rely on this are reproducing upstream behaviour on purpose.
    pub fn julia_string(&self) -> String {
        match *self {
            Value::Missing => "missing".to_string(),
            Value::Int(v) => v.to_string(),
            Value::Float(v) => format_float(v),
            Value::Bool(v) => v.to_string(),
            Value::Str(s) => s.to_string(),
        }
    }

    /// Hashable identity used for grouping, de-duplication and joins.
    pub fn key(&self) -> Key<'a> {
        match *self {
            Value::Missing => Key::Missing,
            Value::Int(v) => Key::Int(v),
            Value::Bool(v) => Key::Int(v as i64),
            Value::Float(v) => float_key(v),
            Value::Str(s) => Key::Str(s),
        }
    }

    /// Total order equivalent to Julia's `isless`.
    pub fn total_cmp(&self, other: &Value<'_>) -> Ordering {
        match (self.rank(), other.rank()) {
            (a, b) if a != b => a.cmp(&b),
            _ => match (self, other) {
                (Value::Str(a), Value::Str(b)) => a.as_bytes().cmp(b.as_bytes()),
                (Value::Int(a), Value::Int(b)) => a.cmp(b),
                (Value::Bool(a), Value::Bool(b)) => a.cmp(b),
                (Value::Missing, Value::Missing) => Ordering::Equal,
                _ => cmp_numbers(
                    self.as_f64().unwrap_or(f64::NAN),
                    other.as_f64().unwrap_or(f64::NAN),
                ),
            },
        }
    }

    /// Coarse type rank: numbers < strings < missing.
    fn rank(&self) -> u8 {
        match self {
            Value::Int(_) | Value::Float(_) | Value::Bool(_) => 0,
            Value::Str(_) => 1,
            Value::Missing => 2,
        }
    }
}

/// `isless` for floats: `NaN` is larger than everything, `-0.0 < 0.0`.
fn cmp_numbers(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => a.total_cmp(&b),
    }
}

/// Hashable cell identity (`isequal` semantics).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Key<'a> {
    Missing,
    Int(i64),
    /// Bit pattern of a non-integral (or out-of-range) float; NaNs unified.
    Float(u64),
    Str(&'a str),
}

fn float_key<'a>(v: f64) -> Key<'a> {
    if v.is_nan() {
        return Key::Float(f64::NAN.to_bits());
    }
    // isequal(1, 1.0) is true, isequal(0.0, -0.0) is false.
    let integral = v.fract() == 0.0 && v.abs() < 9.0e15 && !(v == 0.0 && v.is_sign_negative());
    if integral {
        Key::Int(v as i64)
    } else {
        Key::Float(v.to_bits())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordering_follows_isless() {
        assert_eq!(Value::Int(1).total_cmp(&Value::Float(1.5)), Ordering::Less);
        assert_eq!(
            Value::Float(f64::NAN).total_cmp(&Value::Float(1e300)),
            Ordering::Greater
        );
        assert_eq!(
            Value::Missing.total_cmp(&Value::Str("zzz")),
            Ordering::Greater
        );
        assert_eq!(Value::Str("10").total_cmp(&Value::Str("9")), Ordering::Less);
    }

    #[test]
    fn keys_follow_isequal() {
        assert_eq!(Value::Int(1).key(), Value::Float(1.0).key());
        assert_eq!(Value::Bool(true).key(), Value::Int(1).key());
        assert_ne!(Value::Float(0.0).key(), Value::Float(-0.0).key());
        assert_eq!(Value::Float(f64::NAN).key(), Value::Float(-f64::NAN).key());
    }
}
