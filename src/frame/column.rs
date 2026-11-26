//! Typed, nullable table columns.

use super::value::Value;

/// One homogeneous column. `None` is Julia's `missing`.
#[derive(Debug, Clone, PartialEq)]
pub enum Column {
    Int(Vec<Option<i64>>),
    Float(Vec<Option<f64>>),
    Bool(Vec<Option<bool>>),
    Str(Vec<Option<String>>),
}

impl Column {
    pub fn from_strings<I, S>(values: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Column::Str(values.into_iter().map(|s| Some(s.into())).collect())
    }

    pub fn from_ints<I: IntoIterator<Item = i64>>(values: I) -> Self {
        Column::Int(values.into_iter().map(Some).collect())
    }

    pub fn from_floats<I: IntoIterator<Item = f64>>(values: I) -> Self {
        Column::Float(values.into_iter().map(Some).collect())
    }

    pub fn from_bools<I: IntoIterator<Item = bool>>(values: I) -> Self {
        Column::Bool(values.into_iter().map(Some).collect())
    }

    pub fn len(&self) -> usize {
        match self {
            Column::Int(v) => v.len(),
            Column::Float(v) => v.len(),
            Column::Bool(v) => v.len(),
            Column::Str(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn type_name(&self) -> &'static str {
        match self {
            Column::Int(_) => "Int",
            Column::Float(_) => "Float",
            Column::Bool(_) => "Bool",
            Column::Str(_) => "String",
        }
    }

    pub fn is_numeric(&self) -> bool {
        !matches!(self, Column::Str(_))
    }

    pub fn has_missing(&self) -> bool {
        match self {
            Column::Int(v) => v.iter().any(Option::is_none),
            Column::Float(v) => v.iter().any(Option::is_none),
            Column::Bool(v) => v.iter().any(Option::is_none),
            Column::Str(v) => v.iter().any(Option::is_none),
        }
    }

    /// Cell at `row`.
    pub fn get(&self, row: usize) -> Value<'_> {
        match self {
            Column::Int(v) => v[row].map_or(Value::Missing, Value::Int),
            Column::Float(v) => v[row].map_or(Value::Missing, Value::Float),
            Column::Bool(v) => v[row].map_or(Value::Missing, Value::Bool),
            Column::Str(v) => v[row].as_deref().map_or(Value::Missing, Value::Str),
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = Value<'_>> + '_ {
        (0..self.len()).map(move |row| self.get(row))
    }

    /// New column holding the given rows, in the given order.
    pub fn take(&self, rows: &[usize]) -> Column {
        match self {
            Column::Int(v) => Column::Int(rows.iter().map(|&r| v[r]).collect()),
            Column::Float(v) => Column::Float(rows.iter().map(|&r| v[r]).collect()),
            Column::Bool(v) => Column::Bool(rows.iter().map(|&r| v[r]).collect()),
            Column::Str(v) => Column::Str(rows.iter().map(|&r| v[r].clone()).collect()),
        }
    }

    /// Like [`Column::take`], with `None` producing a missing cell (joins).
    pub fn take_optional(&self, rows: &[Option<usize>]) -> Column {
        match self {
            Column::Int(v) => Column::Int(rows.iter().map(|r| r.and_then(|r| v[r])).collect()),
            Column::Float(v) => Column::Float(rows.iter().map(|r| r.and_then(|r| v[r])).collect()),
            Column::Bool(v) => Column::Bool(rows.iter().map(|r| r.and_then(|r| v[r])).collect()),
            Column::Str(v) => {
                Column::Str(rows.iter().map(|r| r.and_then(|r| v[r].clone())).collect())
            }
        }
    }

    /// Keep rows whose mask entry is `true`.
    pub fn filter(&self, mask: &[bool]) -> Column {
        let rows: Vec<usize> = mask
            .iter()
            .enumerate()
            .filter_map(|(i, keep)| keep.then_some(i))
            .collect();
        self.take(&rows)
    }

    /// Append `other`, promoting the element type like `vcat` does:
    /// Int+Float becomes Float, anything else mixed becomes text.
    pub fn append(self, other: Column) -> Column {
        match (self, other) {
            (Column::Int(mut a), Column::Int(b)) => {
                a.extend(b);
                Column::Int(a)
            }
            (Column::Float(mut a), Column::Float(b)) => {
                a.extend(b);
                Column::Float(a)
            }
            (Column::Bool(mut a), Column::Bool(b)) => {
                a.extend(b);
                Column::Bool(a)
            }
            (Column::Str(mut a), Column::Str(b)) => {
                a.extend(b);
                Column::Str(a)
            }
            (Column::Int(a), Column::Float(mut b)) => {
                let mut out: Vec<Option<f64>> =
                    a.into_iter().map(|v| v.map(|v| v as f64)).collect();
                out.append(&mut b);
                Column::Float(out)
            }
            (Column::Float(mut a), Column::Int(b)) => {
                a.extend(b.into_iter().map(|v| v.map(|v| v as f64)));
                Column::Float(a)
            }
            (a, b) => {
                let mut out = a.into_text();
                out.extend(b.into_text());
                Column::Str(out)
            }
        }
    }

    /// Render every cell as text (missing stays missing).
    pub fn into_text(self) -> Vec<Option<String>> {
        match self {
            Column::Str(v) => v,
            other => other
                .iter()
                .map(|v| (!v.is_missing()).then(|| v.julia_string()))
                .collect(),
        }
    }

    /// Text cells, or `None` when the column is not textual.
    pub fn as_strs(&self) -> Option<&[Option<String>]> {
        match self {
            Column::Str(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_ints(&self) -> Option<&[Option<i64>]> {
        match self {
            Column::Int(v) => Some(v),
            _ => None,
        }
    }

    /// Numeric cells widened to `f64` (`None` for text columns).
    pub fn to_f64(&self) -> Option<Vec<Option<f64>>> {
        self.is_numeric()
            .then(|| self.iter().map(|v| v.as_f64()).collect())
    }
}
