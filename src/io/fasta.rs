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

/// Reads the records of a FASTA text one at a time, for inputs too large
/// to hold as a whole (a genome), from a file or from a stream.
pub struct FastaStream<R: BufRead> {
    reader: R,
    /// Where the text comes from, for error messages.
    origin: String,
    /// Header of the record that follows, read with the record before.
    next_header: Option<String>,
    line: Vec<u8>,
    started: bool,
}

impl<R: BufRead> FastaStream<R> {
    pub fn new(reader: R, origin: impl Into<String>) -> Self {
        FastaStream {
            reader,
            origin: origin.into(),
            next_header: None,
            line: Vec::new(),
            started: false,
        }
    }

    /// Read the next record: its header (without `>`) is returned and its
    /// sequence replaces the content of `sequence`, in upper case and
    /// without line breaks. `None` at the end of the text.
    pub fn next_record(&mut self, sequence: &mut Vec<u8>) -> Result<Option<String>> {
        sequence.clear();
        let mut header = self.next_header.take();
        loop {
            self.line.clear();
            let read = self
                .reader
                .read_until(b'\n', &mut self.line)
                .map_err(|e| Error::io(&self.origin, e))?;
            if read == 0 {
                self.started = true;
                return Ok(header);
            }
            let line = self.line.trim_ascii_end();
            if let Some(text) = line.strip_prefix(b">") {
                let text = String::from_utf8_lossy(text).into_owned();
                if header.is_some() {
                    self.next_header = Some(text);
                    return Ok(header);
                }
                header = Some(text);
                self.started = true;
            } else if !line.is_empty() {
                if header.is_none() {
                    let what = if self.started {
                        "sequence data before a '>' header"
                    } else {
                        "not a FASTA text: it does not start with '>'"
                    };
                    return Err(Error::parse(&self.origin, what));
                }
                let from = sequence.len();
                sequence.extend_from_slice(line);
                sequence[from..].make_ascii_uppercase();
            }
        }
    }
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
    fn a_stream_gives_one_record_at_a_time_in_upper_case() {
        let text = b">chr1 first\nacgtNN\nACgt\r\n\n>chr2\n>chr3\nTT\n";
        let mut stream = FastaStream::new(&text[..], "test");
        let mut sequence = Vec::new();
        let mut records = Vec::new();
        while let Some(header) = stream.next_record(&mut sequence).unwrap() {
            records.push((header, String::from_utf8(sequence.clone()).unwrap()));
        }
        let expected = [("chr1 first", "ACGTNNACGT"), ("chr2", ""), ("chr3", "TT")];
        let found: Vec<(&str, &str)> = records
            .iter()
            .map(|(header, sequence)| (header.as_str(), sequence.as_str()))
            .collect();
        assert_eq!(found, expected);
        // The end stays the end.
        assert!(stream.next_record(&mut sequence).unwrap().is_none());

        let mut empty = FastaStream::new(&b""[..], "test");
        assert!(empty.next_record(&mut sequence).unwrap().is_none());
        let mut other = FastaStream::new(&b"<html>not found</html>\n"[..], "test");
        let error = other.next_record(&mut sequence).unwrap_err().to_string();
        assert!(error.contains("not a FASTA text"), "{error}");
    }

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
