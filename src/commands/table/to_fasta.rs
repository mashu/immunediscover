//! `table fasta`: export sequences of a table as FASTA.

use std::io::Write;
use std::path::PathBuf;

use super::common::{is_match, regex, remove_matches};
use crate::error::{Error, Result};
use crate::frame::{semi_join_rows, Column, Frame};
use crate::info;
use crate::io::compress;
use crate::io::tsv::read_tsv;

#[derive(Debug, Clone)]
pub struct ToFastaParams {
    pub input: PathBuf,
    pub output: PathBuf,
    /// Column with sequence names.
    pub name_column: String,
    /// Sequence column, or comma-separated columns to concatenate.
    pub sequence_columns: String,
    /// Column appended to the FASTA header as description.
    pub description_column: Option<String>,
    /// Keep rows whose name matches. Names are exported unchanged; use
    /// `cleanup` to remove text from them.
    pub filter: Option<String>,
    /// Regex removed from names.
    pub cleanup: Option<String>,
    /// Keep rows whose description matches; capture group 1, if any,
    /// becomes the exported description.
    pub description_filter: Option<String>,
    pub sort_by_name: bool,
    /// Minimum number of donors carrying the allele.
    pub min_cases: usize,
    pub case_column: String,
    /// Export each distinct sequence once (first name wins).
    pub unique_sequences: bool,
}

fn require(table: &Frame, column: &str, context: &str) -> Result<()> {
    if table.has_column(column) {
        Ok(())
    } else {
        Err(Error::invalid(format!(
            "Input file must have {column} column{context}"
        )))
    }
}

/// Rows whose text in `column` matches `pattern`; an empty cell matches
/// nothing.
fn rows_matching(table: &Frame, column: &str, pattern: &str) -> Result<Frame> {
    let regex = regex(pattern)?;
    let mask: Vec<bool> = table
        .column(column)?
        .iter()
        .map(|value| !value.is_missing() && is_match(&regex, &value.julia_string()))
        .collect();
    Ok(table.filter(&mask))
}

/// Keep alleles (name + sequence) observed in at least `min_cases` donors.
fn filter_by_donor_count(
    table: &Frame,
    key: &[String],
    case_column: &str,
    min_cases: usize,
) -> Result<Frame> {
    let groups = table.group_by(key)?;
    let qualifying: Vec<usize> = groups
        .iter()
        .filter(|rows| {
            table
                .take(rows)
                .unique_rows(&[case_column])
                .is_ok_and(|u| u.len() >= min_cases)
        })
        .map(|rows| rows[0])
        .collect();
    let qualifying = table.take(&qualifying).select(key)?;
    if key
        .iter()
        .any(|c| table.column(c).is_ok_and(Column::has_missing))
    {
        return Err(Error::invalid(
            "Missing values in name or sequence columns are not allowed with --mincase",
        ));
    }
    Ok(table.take(&semi_join_rows(table, &qualifying, key)?))
}

/// Export the table as single-line FASTA.
pub fn run(params: &ToFastaParams) -> Result<()> {
    info!("Extracting sequences to FASTA format");
    let mut table = read_tsv(&params.input)?;
    let sequence_columns: Vec<String> = params
        .sequence_columns
        .split(',')
        .map(|c| c.trim().to_string())
        .collect();
    require(&table, &params.name_column, "")?;
    for column in &sequence_columns {
        require(&table, column, "")?;
    }
    if sequence_columns.len() > 1 {
        info!(
            "Using multiple sequence columns for concatenation: {}",
            sequence_columns.join(", ")
        );
    }
    if let Some(column) = &params.description_column {
        require(&table, column, " when coldesc is specified")?;
        info!("Using {column} column for FASTA descriptions");
    }
    info!("Loaded {} rows from input file", table.nrows());

    if let Some(pattern) = &params.filter {
        info!(
            "Filtering rows where {} matches regex pattern '{pattern}'",
            params.name_column
        );
        table = rows_matching(&table, &params.name_column, pattern)?;
        info!("After filtering: {} rows remaining", table.nrows());
    }
    if let Some(pattern) = &params.description_filter {
        let column = params.description_column.as_ref().ok_or_else(|| {
            Error::invalid("Description filter pattern specified but no coldesc column provided")
        })?;
        info!("Filtering rows where {column} matches regex pattern '{pattern}'");
        table = rows_matching(&table, column, pattern)?;
        info!(
            "After description filtering: {} rows remaining",
            table.nrows()
        );
    }

    let mut key = vec![params.name_column.clone()];
    key.extend(sequence_columns.iter().cloned());
    if params.min_cases > 1 {
        require(&table, &params.case_column, " when mincase > 1")?;
        info!(
            "Filtering alleles present in at least {} donors",
            params.min_cases
        );
        table = filter_by_donor_count(&table, &key, &params.case_column, params.min_cases)?;
        info!(
            "After mincase filtering (≥{} donors): {} rows remaining",
            params.min_cases,
            table.nrows()
        );
    }

    let before = table.nrows();
    table = table.unique_by(&key)?;
    if table.nrows() != before {
        info!(
            "Removed {} duplicate records (by {} and sequence columns)",
            before - table.nrows(),
            params.name_column
        );
    }
    if params.unique_sequences {
        info!("Extracting unique sequences (ignoring sequence names)");
        let before = table.nrows();
        table = table.unique_by(&sequence_columns)?;
        if table.nrows() != before {
            info!(
                "Removed {} duplicate sequences with different names",
                before - table.nrows()
            );
        }
    }
    if params.sort_by_name {
        info!("Sorting alleles by name");
        table = table.sorted_asc(&[params.name_column.as_str()])?;
    }

    info!(
        "Writing {} sequences to {}",
        table.nrows(),
        params.output.display()
    );
    let cleanup = params.cleanup.as_deref().map(regex).transpose()?;
    let description_filter = params
        .description_filter
        .as_deref()
        .map(regex)
        .transpose()?;
    let sequences: Vec<&Column> = sequence_columns
        .iter()
        .map(|c| table.column(c))
        .collect::<Result<_>>()?;
    let names = table.column(&params.name_column)?;
    let descriptions = params
        .description_column
        .as_ref()
        .map(|c| table.column(c))
        .transpose()?;

    let mut out = compress::create(&params.output, false)?;
    for row in 0..table.nrows() {
        let mut name = names.get(row).julia_string();
        let sequence: String = sequences
            .iter()
            .map(|c| c.get(row))
            .filter(|cell| !cell.is_missing())
            .map(|cell| cell.julia_string())
            .collect();
        if let Some(regex) = cleanup.as_ref().filter(|r| is_match(r, &name)) {
            name = remove_matches(regex, &name).trim().to_string();
        }
        if let Some(column) = descriptions {
            let cell = column.get(row);
            let description = if cell.is_missing() {
                String::new()
            } else {
                cell.julia_string()
            };
            if !description.is_empty() {
                match &description_filter {
                    None => name = format!("{name} {description}"),
                    Some(regex) => match regex.captures(&description) {
                        Ok(Some(captures)) => {
                            if let Some(group) = captures.get(1) {
                                name = format!("{name} {}", group.as_str());
                            }
                        }
                        _ => continue,
                    },
                }
            }
        }
        writeln!(out, ">{name}\n{sequence}").map_err(|e| Error::io(&params.output, e))?;
    }
    out.flush().map_err(|e| Error::io(&params.output, e))?;
    info!(
        "Successfully extracted {} unique sequences to {}",
        table.nrows(),
        params.output.display()
    );
    Ok(())
}
