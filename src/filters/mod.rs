//! Composable row filters for discovery pipelines.
//!
//! A [`GermlineFilter`] applies a sequence of [`RowFilter`] criteria to a
//! table and logs how many rows each one kept, so a run's log documents
//! where candidates were lost.

use crate::error::{Error, Result};
use crate::frame::{Frame, Value};
use crate::info;

/// One filtering criterion: decides, per row, whether the row is kept.
pub trait RowFilter {
    /// Human-readable name used in the log.
    fn label(&self) -> &str;

    /// `true` for every row that passes.
    fn mask(&self, frame: &Frame) -> Result<Vec<bool>>;
}

fn numeric(frame: &Frame, column: &str) -> Result<Vec<f64>> {
    frame
        .column(column)?
        .iter()
        .map(|value| {
            value.as_f64().ok_or_else(|| {
                Error::invalid(format!(
                    "Column '{column}' holds a non-numeric or missing value and cannot be filtered on"
                ))
            })
        })
        .collect()
}

/// Keep rows with `column >= value`.
pub struct MinThreshold {
    pub column: String,
    pub value: f64,
    pub label: String,
}

impl MinThreshold {
    pub fn new(column: &str, value: f64, label: &str) -> Self {
        MinThreshold {
            column: column.to_string(),
            value,
            label: label.to_string(),
        }
    }
}

impl RowFilter for MinThreshold {
    fn label(&self) -> &str {
        &self.label
    }

    fn mask(&self, frame: &Frame) -> Result<Vec<bool>> {
        Ok(numeric(frame, &self.column)?
            .iter()
            .map(|v| *v >= self.value)
            .collect())
    }
}

/// Keep rows with `column <= value`.
pub struct MaxThreshold {
    pub column: String,
    pub value: f64,
    pub label: String,
}

impl MaxThreshold {
    pub fn new(column: &str, value: f64, label: &str) -> Self {
        MaxThreshold {
            column: column.to_string(),
            value,
            label: label.to_string(),
        }
    }
}

impl RowFilter for MaxThreshold {
    fn label(&self) -> &str {
        &self.label
    }

    fn mask(&self, frame: &Frame) -> Result<Vec<bool>> {
        Ok(numeric(frame, &self.column)?
            .iter()
            .map(|v| *v <= self.value)
            .collect())
    }
}

/// Keep rows whose text in `column` is at least `value` characters long.
pub struct MinStringLength {
    pub column: String,
    pub value: usize,
    pub label: String,
}

impl MinStringLength {
    pub fn new(column: &str, value: usize, label: &str) -> Self {
        MinStringLength {
            column: column.to_string(),
            value,
            label: label.to_string(),
        }
    }
}

impl RowFilter for MinStringLength {
    fn label(&self) -> &str {
        &self.label
    }

    fn mask(&self, frame: &Frame) -> Result<Vec<bool>> {
        frame
            .column(&self.column)?
            .iter()
            .map(|value| match value {
                Value::Str(s) => Ok(s.chars().count() >= self.value),
                _ => Err(Error::invalid(format!(
                    "Column '{}' holds a non-text or missing value",
                    self.column
                ))),
            })
            .collect()
    }
}

/// Keep rows with `column >= 0`.
pub struct NonNegative {
    pub column: String,
    pub label: String,
}

impl NonNegative {
    pub fn new(column: &str, label: &str) -> Self {
        NonNegative {
            column: column.to_string(),
            label: label.to_string(),
        }
    }
}

impl RowFilter for NonNegative {
    fn label(&self) -> &str {
        &self.label
    }

    fn mask(&self, frame: &Frame) -> Result<Vec<bool>> {
        Ok(numeric(frame, &self.column)?
            .iter()
            .map(|v| *v >= 0.0)
            .collect())
    }
}

/// Keep rows accepted by an arbitrary predicate over `(table, row index)`.
pub struct CustomFilter<F> {
    pub predicate: F,
    pub label: String,
}

impl<F> CustomFilter<F>
where
    F: Fn(&Frame, usize) -> Result<bool>,
{
    pub fn new(label: &str, predicate: F) -> Self {
        CustomFilter {
            predicate,
            label: label.to_string(),
        }
    }
}

impl<F> RowFilter for CustomFilter<F>
where
    F: Fn(&Frame, usize) -> Result<bool>,
{
    fn label(&self) -> &str {
        &self.label
    }

    fn mask(&self, frame: &Frame) -> Result<Vec<bool>> {
        (0..frame.nrows())
            .map(|row| (self.predicate)(frame, row))
            .collect()
    }
}

/// An ordered chain of criteria applied one after another.
#[derive(Default)]
pub struct GermlineFilter<'a> {
    criteria: Vec<Box<dyn RowFilter + 'a>>,
}

impl<'a> GermlineFilter<'a> {
    pub fn new() -> Self {
        GermlineFilter {
            criteria: Vec::new(),
        }
    }

    /// Append a criterion (builder style).
    pub fn with(mut self, criterion: impl RowFilter + 'a) -> Self {
        self.criteria.push(Box::new(criterion));
        self
    }

    pub fn push(&mut self, criterion: impl RowFilter + 'a) {
        self.criteria.push(Box::new(criterion));
    }

    pub fn is_empty(&self) -> bool {
        self.criteria.is_empty()
    }

    /// Filter `frame` in place, logging kept/total after every criterion.
    pub fn apply(&self, frame: &mut Frame) -> Result<()> {
        for criterion in &self.criteria {
            let before = frame.nrows();
            let mask = criterion.mask(frame)?;
            *frame = frame.filter(&mask);
            info!(
                "{}: {}/{} rows kept",
                criterion.label(),
                frame.nrows(),
                before
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    fn sample() -> Frame {
        Frame::from_columns([
            ("count", Column::from_ints([10, 5, 2, 8, 1])),
            ("ratio", Column::from_floats([1.0, 0.5, 0.2, 0.8, 0.1])),
            (
                "qseq",
                Column::from_strings([
                    "ATCGATCGATCG",
                    "ATCG",
                    "ATCGATCGATCGATCG",
                    "ATCGATCG",
                    "AT",
                ]),
            ),
            ("mismatch", Column::from_ints([0, 3, 5, 1, -1])),
        ])
        .unwrap()
    }

    fn kept(filter: GermlineFilter<'_>) -> usize {
        let mut frame = sample();
        filter.apply(&mut frame).unwrap();
        frame.nrows()
    }

    // Same expectations as the upstream `filters.jl` test set.
    #[test]
    fn single_criteria() {
        assert_eq!(
            kept(GermlineFilter::new().with(MinThreshold::new("count", 5.0, "Min count"))),
            3
        );
        assert_eq!(
            kept(GermlineFilter::new().with(MaxThreshold::new("mismatch", 3.0, "Max mismatch"))),
            4
        );
        assert_eq!(
            kept(GermlineFilter::new().with(MinStringLength::new("qseq", 8, "Min read length"))),
            3
        );
        assert_eq!(
            kept(GermlineFilter::new().with(NonNegative::new("mismatch", "Non-negative"))),
            4
        );
        let custom = CustomFilter::new("High ratio", |f: &Frame, row| {
            Ok(f.value(row, "ratio")?.as_f64().unwrap() > 0.3)
        });
        assert_eq!(kept(GermlineFilter::new().with(custom)), 3);
    }

    #[test]
    fn criteria_compose_in_order() {
        let chain = GermlineFilter::new()
            .with(MinThreshold::new("count", 5.0, "Min count"))
            .with(MaxThreshold::new("mismatch", 3.0, "Max mismatch"));
        assert_eq!(kept(chain), 3);
    }
}
