//! `table transform`: rewrite a column from regex capture groups.

use std::path::Path;

use super::common::{load, regex, save_gz};
use crate::error::Result;
use crate::frame::{Column, Frame};
use crate::info;

/// Rewrite one cell: the first match of `regex` is replaced by
/// `replacement`, in which `\\1`, `\\2`, ... stand for the capture groups;
/// text around the match is kept. Returns `(new value, concatenated
/// captures)`; a cell without a match is returned unchanged.
fn rewrite(regex: &fancy_regex::Regex, value: &str, replacement: &str) -> (String, String) {
    let Ok(Some(captures)) = regex.captures(value) else {
        return (value.to_string(), String::new());
    };
    let matched = captures.get(0).expect("group 0 is the whole match");
    let mut substitute = replacement.to_string();
    let mut extracted = String::new();
    // Highest group first, so that `\\1` does not eat the start of `\\10`.
    for group in (1..captures.len()).rev() {
        let text = captures.get(group).map_or("", |c| c.as_str());
        substitute = substitute.replace(&format!("\\{group}"), text);
    }
    for group in 1..captures.len() {
        if let Some(capture) = captures.get(group) {
            extracted.push_str(capture.as_str());
        }
    }
    let result = format!(
        "{}{}{}",
        &value[..matched.start()],
        substitute,
        &value[matched.end()..]
    );
    (result, extracted)
}

/// Transform `column` of `input`; optionally store the captured text in
/// `new_column`. The column becomes text; missing cells stay missing.
pub fn run(
    input: &Path,
    output: &Path,
    column: &str,
    pattern: &str,
    replacement: &str,
    new_column: Option<&str>,
) -> Result<Frame> {
    let mut table = load(input)?;
    info!("Transforming column '{column}' using pattern: {pattern}");
    let regex = regex(pattern)?;
    let original: Vec<Option<String>> = table
        .column(column)?
        .iter()
        .map(|v| (!v.is_missing()).then(|| v.julia_string()))
        .collect();
    let rewritten: Vec<Option<(String, String)>> = original
        .iter()
        .map(|cell| {
            cell.as_deref()
                .map(|value| rewrite(&regex, value, replacement))
        })
        .collect();

    let changed = original
        .iter()
        .zip(&rewritten)
        .filter(|(before, after)| before.as_deref() != after.as_ref().map(|(new, _)| new.as_str()))
        .count();
    info!("Transformed {changed} values in column '{column}'");
    let (transformed, extracted): (Vec<Option<String>>, Vec<String>) = rewritten
        .into_iter()
        .map(|cell| match cell {
            Some((new, captured)) => (Some(new), captured),
            None => (None, String::new()),
        })
        .unzip();
    table.set_column(column, Column::Str(transformed))?;
    if let Some(name) = new_column {
        table.set_column(name, Column::from_strings(extracted))?;
        info!("Added new column '{name}'");
    }
    let saved = save_gz(output, &table)?;
    info!("Transformed data saved to: {}", saved.display());
    Ok(table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_the_match_is_replaced() {
        let regex = regex("(.*)_(.*)").unwrap();
        assert_eq!(
            rewrite(&regex, "ABC_123", "\\1-\\2"),
            ("ABC-123".into(), "ABC123".into())
        );
        let digits = super::regex("(\\d+)").unwrap();
        assert_eq!(rewrite(&digits, "ID_123_x", "<\\1>").0, "ID_<123>_x");
        assert_eq!(
            rewrite(&digits, "none", "\\1"),
            ("none".into(), String::new())
        );
        // A group that did not take part in the match expands to nothing.
        let optional = super::regex("a(b)?c").unwrap();
        assert_eq!(rewrite(&optional, "xacx", "[\\1]").0, "x[]x");
    }
}
