//! The demultiplexed read table shared by all search and discovery commands.

use std::path::Path;

use crate::error::{Error, Result};
use crate::frame::{Column, Frame};
use crate::info;
use crate::io::tsv::{read_table, ReadOptions};

/// Columns every demultiplexed table must provide.
pub const REQUIRED_COLUMNS: [&str; 4] = ["well", "case", "name", "genomic_sequence"];

/// Reads with their plate well and donor (case) assignment.
#[derive(Debug, Clone)]
pub struct DemuxTable {
    /// Well identifiers rendered as text (`string(well)`).
    pub well: Vec<String>,
    /// Donor identifiers; always read as text, never as numbers.
    pub case: Vec<String>,
    pub name: Vec<String>,
    pub sequence: Vec<String>,
    /// The table as read, with inferred column types.
    pub frame: Frame,
}

impl DemuxTable {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Validate and wrap an already loaded table.
    pub fn from_frame(frame: Frame) -> Result<Self> {
        if REQUIRED_COLUMNS.iter().any(|c| !frame.has_column(c)) {
            return Err(Error::invalid(
                "File must contain following columns: well, case, name, genomic_sequence",
            ));
        }
        let text = |name: &str| -> Result<Vec<String>> {
            Ok(frame
                .column(name)?
                .iter()
                .map(|v| v.julia_string())
                .collect())
        };
        let sequence = match frame.column("genomic_sequence")? {
            Column::Str(cells) => cells
                .iter()
                .enumerate()
                .map(|(row, cell)| {
                    cell.clone().ok_or_else(|| {
                        Error::invalid(format!("Row {}: genomic_sequence is empty", row + 1))
                    })
                })
                .collect::<Result<Vec<String>>>()?,
            other => {
                return Err(Error::ColumnType {
                    column: "genomic_sequence".to_string(),
                    expected: "String",
                    found: other.type_name(),
                })
            }
        };
        Ok(DemuxTable {
            well: text("well")?,
            case: text("case")?,
            name: text("name")?,
            sequence,
            frame,
        })
    }

    /// First `n` reads.
    pub fn head(&self, n: usize) -> Result<Self> {
        Self::from_frame(self.frame.head(n))
    }
}

/// Load a demultiplexed TSV (`well`, `case`, `name`, `genomic_sequence`).
pub fn load_demultiplex(path: &Path, limit: Option<usize>) -> Result<DemuxTable> {
    let options = ReadOptions::tsv().text_column("case").limit(limit);
    let table = DemuxTable::from_frame(read_table(path, &options)?)?;
    info!("Reading {}.", path.display());
    Ok(table)
}
