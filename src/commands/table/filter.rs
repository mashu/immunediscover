//! `table filter`: keep rows by regex match or numeric comparison.

use std::path::Path;
use std::str::FromStr;

use super::common::{is_match, load, regex, save_gz};
use crate::error::{Error, Result};
use crate::frame::{Column, Frame, Value};
use crate::info;

/// Numeric comparison operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparison {
    Less,
    LessEqual,
    GreaterEqual,
    Greater,
}

impl Comparison {
    pub const SYMBOLS: [&'static str; 4] = ["<", "<=", ">=", ">"];

    fn holds(self, value: f64, threshold: f64) -> bool {
        match self {
            Comparison::Less => value < threshold,
            Comparison::LessEqual => value <= threshold,
            Comparison::GreaterEqual => value >= threshold,
            Comparison::Greater => value > threshold,
        }
    }
}

impl FromStr for Comparison {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "<" => Ok(Comparison::Less),
            "<=" => Ok(Comparison::LessEqual),
            ">=" => Ok(Comparison::GreaterEqual),
            ">" => Ok(Comparison::Greater),
            other => Err(Error::invalid(format!(
                "Invalid operator: {other}. Must be one of: <, <=, >=, >"
            ))),
        }
    }
}

/// What rows to keep.
#[derive(Debug, Clone)]
pub enum RowCondition {
    /// The cell text matches the regular expression.
    Pattern(String),
    /// The numeric cell compares to the threshold.
    Numeric(Comparison, f64),
}

impl RowCondition {
    /// Build from the optional command-line values, enforcing that exactly
    /// one filtering method is given.
    pub fn from_options(
        pattern: Option<String>,
        operator: Option<Comparison>,
        threshold: Option<f64>,
    ) -> Result<Self> {
        match (pattern, operator, threshold) {
            (Some(_), Some(_), _) | (Some(_), _, Some(_)) => Err(Error::invalid(
                "Cannot use both regex pattern and numeric operator/threshold. Choose one filtering method.",
            )),
            (Some(pattern), None, None) => Ok(RowCondition::Pattern(pattern)),
            (None, Some(operator), Some(threshold)) => Ok(RowCondition::Numeric(operator, threshold)),
            _ => Err(Error::invalid(
                "Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering",
            )),
        }
    }
}

fn numeric_mask(
    column: &Column,
    name: &str,
    comparison: Comparison,
    threshold: f64,
) -> Result<Vec<bool>> {
    if column.is_numeric() && column.iter().any(|v| !v.is_missing()) {
        // Missing values never pass.
        return Ok(column
            .iter()
            .map(|v| v.as_f64().is_some_and(|x| comparison.holds(x, threshold)))
            .collect());
    }
    // Text column: every cell must parse as a number.
    column
        .iter()
        .map(|value| match value {
            Value::Str(text) => text.trim().parse::<f64>().ok(),
            _ => None,
        })
        .map(|parsed| {
            parsed
                .map(|x| comparison.holds(x, threshold))
                .ok_or_else(|| {
                    Error::invalid(format!(
                        "Column '{name}' cannot be converted to numeric values for filtering"
                    ))
                })
        })
        .collect()
}

/// Filter rows of `input` on `column` and write the result.
pub fn run(input: &Path, output: &Path, column: &str, condition: &RowCondition) -> Result<Frame> {
    let table = load(input)?;
    let cells = table.column(column)?;
    let mask: Vec<bool> = match condition {
        RowCondition::Pattern(pattern) => {
            info!("Filtering column '{column}' using regex pattern: {pattern}");
            let regex = regex(pattern)?;
            cells
                .iter()
                .map(|v| !v.is_missing() && is_match(&regex, &v.julia_string()))
                .collect()
        }
        RowCondition::Numeric(comparison, threshold) => {
            info!(
                "Filtering column '{column}' using numeric operation: {comparison:?} {threshold}"
            );
            numeric_mask(cells, column, *comparison, *threshold)?
        }
    };
    let filtered = table.filter(&mask);
    info!(
        "Filter results: {} rows kept, {} rows removed",
        filtered.nrows(),
        table.nrows() - filtered.nrows()
    );
    let saved = save_gz(output, &filtered)?;
    info!("Filtered data saved to: {}", saved.display());
    Ok(filtered)
}
