//! FASTA reading and writing.
//!
//! Header handling follows FASTX.jl: the *description* is the whole header
//! line (trailing blanks included), the *identifier* is its first
//! whitespace-delimited token.

use std::io::{BufRead, Write};
use std::path::Path;

use super::compress;
use crate::error::{Error, Result};

/// Line width used by FASTX.jl's `FASTA.Writer`.
pub const DEFAULT_LINE_WIDTH: usize = 70;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// Full header line without the leading `>`.
    pub description: String,
    pub sequence: String,
}

impl FastaRecord {
    pub fn new(description: impl Into<String>, sequence: impl Into<String>) -> Self {
        FastaRecord {
            description: description.into(),
            sequence: sequence.into(),
        }
    }

    /// Header up to the first whitespace character.
    pub fn identifier(&self) -> &str {
        identifier_of(&self.description)
    }

    /// Header text after the identifier and one separating character.
    pub fn comment(&self) -> &str {
        let id_len = self.identifier().len();
        self.description.get(id_len + 1..).unwrap_or("")
    }
}

/// First whitespace-delimited token of a header.
pub fn identifier_of(description: &str) -> &str {
    description
        .split(|c: char| c.is_ascii_whitespace())
        .next()
        .unwrap_or("")
}

/// Read every record of a (possibly gzip-compressed) FASTA file.
pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let reader = compress::open_maybe_gzip(path)?;
    let mut records: Vec<FastaRecord> = Vec::new();
    for (number, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| Error::io(path, e))?;
        let line = line.strip_suffix('\r').unwrap_or(&line);
        if let Some(header) = line.strip_prefix('>') {
            records.push(FastaRecord::new(header, String::new()));
        } else if line.is_empty() {
            continue;
        } else {
            let record = records.last_mut().ok_or_else(|| {
                Error::parse(
                    path,
                    format!("line {}: sequence data before first '>' header", number + 1),
                )
            })?;
            record.sequence.push_str(line);
        }
    }
    Ok(records)
}

/// Streaming FASTA writer.
pub struct FastaWriter<W: Write> {
    out: W,
    line_width: Option<usize>,
}

impl<W: Write> FastaWriter<W> {
    /// `line_width = None` writes each sequence on a single line.
    pub fn new(out: W, line_width: Option<usize>) -> Self {
        FastaWriter { out, line_width }
    }

    pub fn write(&mut self, description: &str, sequence: &str) -> std::io::Result<()> {
        writeln!(self.out, ">{description}")?;
        match self.line_width {
            Some(width) if width > 0 => {
                for chunk in sequence.as_bytes().chunks(width) {
                    self.out.write_all(chunk)?;
                    self.out.write_all(b"\n")?;
                }
            }
            _ => writeln!(self.out, "{sequence}")?,
        }
        Ok(())
    }

    pub fn finish(mut self) -> std::io::Result<()> {
        self.out.flush()
    }
}

/// Write records to `path`, wrapping sequence lines at `line_width`.
pub fn write_fasta<'a, I>(path: &Path, records: I, line_width: Option<usize>) -> Result<()>
where
    I: IntoIterator<Item = (&'a str, &'a str)>,
{
    let mut writer = FastaWriter::new(compress::create(path, false)?, line_width);
    for (description, sequence) in records {
        writer
            .write(description, sequence)
            .map_err(|e| Error::io(path, e))?;
    }
    writer.finish().map_err(|e| Error::io(path, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn headers_keep_trailing_blanks_like_fastx() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("f.fa");
        std::fs::write(
            &path,
            b">A*01 desc here  \nACGT\nAC\n>B*01\t\nGG\r\n>C*01 \nTT\n",
        )
        .unwrap();
        let records = read_fasta(&path).unwrap();
        assert_eq!(records[0].description, "A*01 desc here  ");
        assert_eq!(records[0].identifier(), "A*01");
        assert_eq!(records[0].sequence, "ACGTAC");
        assert_eq!(records[1].description, "B*01\t");
        assert_eq!(records[1].sequence, "GG");
        assert_eq!(records[2].description, "C*01 ");
    }

    #[test]
    fn writer_wraps_at_width() {
        let mut buffer = Vec::new();
        let mut writer = FastaWriter::new(&mut buffer, Some(DEFAULT_LINE_WIDTH));
        writer.write("name here", &"A".repeat(150)).unwrap();
        writer.write("n2", "").unwrap();
        let text = String::from_utf8(buffer).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], ">name here");
        assert_eq!(
            (lines[1].len(), lines[2].len(), lines[3].len()),
            (70, 70, 10)
        );
        assert_eq!(lines[4], ">n2");
        assert_eq!(lines.len(), 5);
    }
}
