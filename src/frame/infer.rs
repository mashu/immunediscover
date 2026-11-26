//! Column type inference for delimited text, following CSV.jl 0.10.
//!
//! CSV.jl detects a column's type from its first value and promotes it as
//! later values disagree, so the result can depend on value order:
//!
//! * the first value decides between integer, float, boolean and text;
//! * an integer column becomes float when a float appears;
//! * any other disagreement turns the column into text (cells verbatim).
//!
//! Numeric cells may be padded with spaces and carry a sign or leading
//! zeros (`007` reads as 7). Booleans are `true/True/TRUE/T` and
//! `false/False/FALSE/F` (and `1`/`0` once a column is boolean).

use super::column::Column;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Int,
    Float,
    Bool,
    Text,
}

fn parse_int(cell: &str) -> Option<i64> {
    let trimmed = cell.trim_matches(' ');
    let digits = trimmed.strip_prefix(['+', '-']).unwrap_or(trimmed);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    trimmed.strip_prefix('+').unwrap_or(trimmed).parse().ok()
}

fn parse_float(cell: &str) -> Option<f64> {
    let trimmed = cell.trim_matches(' ');
    let unsigned = trimmed.strip_prefix(['+', '-']).unwrap_or(trimmed);
    let special = ["nan", "inf", "infinity"]
        .iter()
        .any(|s| unsigned.eq_ignore_ascii_case(s));
    let looks_numeric = unsigned
        .bytes()
        .next()
        .is_some_and(|b| b.is_ascii_digit() || b == b'.');
    if !(special || looks_numeric) {
        return None;
    }
    trimmed.parse().ok()
}

fn parse_bool(cell: &str) -> Option<bool> {
    match cell.trim_matches(' ') {
        "true" | "True" | "TRUE" | "T" | "1" => Some(true),
        "false" | "False" | "FALSE" | "F" | "0" => Some(false),
        _ => None,
    }
}

/// Type suggested by a single value.
fn detect(cell: &str) -> Kind {
    if parse_int(cell).is_some() {
        Kind::Int
    } else if parse_float(cell).is_some() {
        Kind::Float
    } else if parse_bool(cell).is_some() {
        Kind::Bool
    } else {
        Kind::Text
    }
}

/// Column type after seeing `cell` in a column currently typed `kind`.
fn promote(kind: Kind, cell: &str) -> Kind {
    match kind {
        Kind::Int if parse_int(cell).is_some() => Kind::Int,
        Kind::Int | Kind::Float if parse_float(cell).is_some() => Kind::Float,
        Kind::Bool if parse_bool(cell).is_some() => Kind::Bool,
        _ => Kind::Text,
    }
}

/// Turn raw cells (`None` = empty, unquoted field) into a typed column.
///
/// With `force_text` the cells are kept verbatim (CSV.jl `types=String`).
pub fn infer_column(cells: Vec<Option<String>>, force_text: bool) -> Column {
    if force_text {
        return Column::Str(cells);
    }
    let mut kind: Option<Kind> = None;
    for cell in cells.iter().flatten() {
        kind = Some(match kind {
            None => detect(cell),
            Some(kind) => promote(kind, cell),
        });
        if kind == Some(Kind::Text) {
            break;
        }
    }
    match kind {
        Some(Kind::Int) => Column::Int(
            cells
                .iter()
                .map(|c| c.as_deref().and_then(parse_int))
                .collect(),
        ),
        Some(Kind::Float) => Column::Float(
            cells
                .iter()
                .map(|c| c.as_deref().and_then(parse_float))
                .collect(),
        ),
        Some(Kind::Bool) => Column::Bool(
            cells
                .iter()
                .map(|c| c.as_deref().and_then(parse_bool))
                .collect(),
        ),
        // Text, or a column with no values at all.
        Some(Kind::Text) | None => Column::Str(cells),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn infer(cells: &[&str]) -> Column {
        infer_column(
            cells
                .iter()
                .map(|c| (!c.is_empty()).then(|| c.to_string()))
                .collect(),
            false,
        )
    }

    // Expectations mirror CSV.jl 0.10.15 on the same inputs.
    #[test]
    fn integers_allow_padding_sign_and_leading_zeros() {
        assert_eq!(infer(&["007", "12", "1"]), Column::from_ints([7, 12, 1]));
        assert_eq!(infer(&["+5", "-3", "3"]), Column::from_ints([5, -3, 3]));
        assert_eq!(infer(&[" 5", "6 ", "7"]), Column::from_ints([5, 6, 7]));
    }

    #[test]
    fn floats_include_special_values_and_promote_ints() {
        assert_eq!(
            infer(&["1.0", "2", "3.50"]),
            Column::from_floats([1.0, 2.0, 3.5])
        );
        assert_eq!(
            infer(&["1e5", ".5", "5."]),
            Column::from_floats([1e5, 0.5, 5.0])
        );
        assert_eq!(infer(&["1", "2.5"]), Column::from_floats([1.0, 2.5]));
        assert_eq!(
            infer(&["99999999999999999999", "1"]),
            Column::from_floats([1e20, 1.0])
        );
        let special = infer(&["NaN", "Inf", "-Inf"]);
        let values = special.to_f64().unwrap();
        assert!(values[0].unwrap().is_nan());
        assert_eq!(values[1], Some(f64::INFINITY));
        assert_eq!(values[2], Some(f64::NEG_INFINITY));
    }

    #[test]
    fn booleans_accept_letter_forms() {
        assert_eq!(infer(&["true", "false"]), Column::from_bools([true, false]));
        assert_eq!(
            infer(&["T", "F", "T"]),
            Column::from_bools([true, false, true])
        );
        // Once boolean, 1/0 are accepted too.
        assert_eq!(
            infer(&["T", "1", "0"]),
            Column::from_bools([true, true, false])
        );
    }

    #[test]
    fn detection_depends_on_value_order() {
        // An integer column does not become boolean.
        assert!(matches!(infer(&["1", "T"]), Column::Str(_)));
        assert!(matches!(infer(&["1", "abc"]), Column::Str(_)));
        assert!(matches!(infer(&["T", "A"]), Column::Str(_)));
        assert!(matches!(infer(&["1.5", "true"]), Column::Str(_)));
    }

    #[test]
    fn missing_cells_are_preserved() {
        assert_eq!(
            infer(&["1", "", "3"]),
            Column::Int(vec![Some(1), None, Some(3)])
        );
        assert_eq!(infer(&["", ""]), Column::Str(vec![None, None]));
    }

    #[test]
    fn forced_text_keeps_cells_verbatim() {
        let c = infer_column(vec![Some("007".into()), Some("8".into())], true);
        assert_eq!(c, Column::from_strings(["007", "8"]));
    }
}
