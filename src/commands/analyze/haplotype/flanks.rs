//! Reassembling an allele with its reported flanks.
//!
//! Two rows with the same allele name but different flanks are different
//! variants; the "extended sequence" (flank + allele + flank) tells them
//! apart. Which flank columns exist depends on how the search was run.

use crate::error::Result;
use crate::frame::{Frame, Value};
use crate::seq::GeneType;

/// Flank columns present in a search result table.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FlankColumns {
    prefix: bool,
    suffix: bool,
    heptamer: bool,
    spacer: bool,
    nonamer: bool,
    pre_heptamer: bool,
    pre_spacer: bool,
    pre_nonamer: bool,
    post_heptamer: bool,
    post_spacer: bool,
    post_nonamer: bool,
}

impl FlankColumns {
    pub fn detect(table: &Frame) -> Self {
        let has = |name: &str| table.has_column(name);
        FlankColumns {
            prefix: has("prefix"),
            suffix: has("suffix"),
            heptamer: has("heptamer"),
            spacer: has("spacer"),
            nonamer: has("nonamer"),
            pre_heptamer: has("pre_heptamer"),
            pre_spacer: has("pre_spacer"),
            pre_nonamer: has("pre_nonamer"),
            post_heptamer: has("post_heptamer"),
            post_spacer: has("post_spacer"),
            post_nonamer: has("post_nonamer"),
        }
    }
}

/// Text of a cell; missing reads as empty.
fn cell(table: &Frame, row: usize, column: &str) -> Result<String> {
    Ok(match table.value(row, column)? {
        Value::Missing => String::new(),
        other => other.julia_string(),
    })
}

/// Concatenate the listed columns that are present.
fn concat(table: &Frame, row: usize, columns: &[(bool, &str)]) -> Result<String> {
    let mut out = String::new();
    for (present, name) in columns {
        if *present {
            out.push_str(&cell(table, row, name)?);
        }
    }
    Ok(out)
}

fn upstream(table: &Frame, row: usize, present: &FlankColumns, gene: GeneType) -> Result<String> {
    if present.prefix {
        return cell(table, row, "prefix");
    }
    match gene {
        GeneType::V => Ok(String::new()),
        GeneType::J => concat(
            table,
            row,
            &[
                (present.nonamer, "nonamer"),
                (present.spacer, "spacer"),
                (present.heptamer, "heptamer"),
            ],
        ),
        GeneType::D => concat(
            table,
            row,
            &[
                (present.pre_nonamer, "pre_nonamer"),
                (present.pre_spacer, "pre_spacer"),
                (present.pre_heptamer, "pre_heptamer"),
            ],
        ),
    }
}

fn downstream(table: &Frame, row: usize, present: &FlankColumns, gene: GeneType) -> Result<String> {
    if present.suffix {
        return cell(table, row, "suffix");
    }
    match gene {
        GeneType::V => concat(
            table,
            row,
            &[
                (present.heptamer, "heptamer"),
                (present.spacer, "spacer"),
                (present.nonamer, "nonamer"),
            ],
        ),
        GeneType::J => Ok(String::new()),
        GeneType::D => concat(
            table,
            row,
            &[
                (present.post_heptamer, "post_heptamer"),
                (present.post_spacer, "post_spacer"),
                (present.post_nonamer, "post_nonamer"),
            ],
        ),
    }
}

/// Upstream flank + allele sequence + downstream flank for one row.
/// Flanks are only added when the gene name reveals the segment type.
pub fn extended_sequence(
    table: &Frame,
    row: usize,
    gene_column: &str,
    present: &FlankColumns,
) -> Result<String> {
    let gene_name = cell(table, row, gene_column)?;
    let sequence = cell(table, row, "sequence")?;
    match GeneType::from_name(&gene_name) {
        Some(gene) => Ok(format!(
            "{}{}{}",
            upstream(table, row, present, gene)?,
            sequence,
            downstream(table, row, present, gene)?
        )),
        None => Ok(sequence),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frame::Column;

    #[test]
    fn v_rows_append_the_signal_and_prefer_explicit_suffix() {
        let rss = Frame::from_columns([
            ("gene", Column::from_strings(["IGHV1"])),
            ("sequence", Column::from_strings(["SEQ"])),
            ("prefix", Column::from_strings(["pp"])),
            ("heptamer", Column::from_strings(["H"])),
            ("nonamer", Column::from_strings(["N"])),
        ])
        .unwrap();
        let present = FlankColumns::detect(&rss);
        assert_eq!(
            extended_sequence(&rss, 0, "gene", &present).unwrap(),
            "ppSEQHN"
        );

        let extension = Frame::from_columns([
            ("gene", Column::from_strings(["IGHJ4", "TRAC"])),
            ("sequence", Column::from_strings(["SEQ", "SEQ"])),
            ("prefix", Column::Str(vec![None, Some("x".into())])),
            ("suffix", Column::from_strings(["ss", "y"])),
        ])
        .unwrap();
        let present = FlankColumns::detect(&extension);
        assert_eq!(
            extended_sequence(&extension, 0, "gene", &present).unwrap(),
            "SEQss"
        );
        // Unknown segment type: no flanks.
        assert_eq!(
            extended_sequence(&extension, 1, "gene", &present).unwrap(),
            "SEQ"
        );
    }
}
