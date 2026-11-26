//! FASTQ reading and writing (four-line records).

use std::io::{BufRead, Write};
use std::path::Path;

use super::compress;
use super::fasta::identifier_of;
use crate::error::{Error, Result};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastqRecord {
    /// Full header line without the leading `@`.
    pub description: String,
    pub sequence: String,
    pub quality: String,
}

impl FastqRecord {
    /// Header up to the first whitespace character.
    pub fn identifier(&self) -> &str {
        identifier_of(&self.description)
    }

    pub fn write_to(&self, out: &mut dyn Write) -> std::io::Result<()> {
        writeln!(
            out,
            "@{}\n{}\n+\n{}",
            self.description, self.sequence, self.quality
        )
    }
}

fn read_line(reader: &mut dyn BufRead, buffer: &mut String, path: &Path) -> Result<bool> {
    buffer.clear();
    let n = reader.read_line(buffer).map_err(|e| Error::io(path, e))?;
    while buffer.ends_with(['\n', '\r']) {
        buffer.pop();
    }
    Ok(n > 0)
}

/// Call `visit` for every record of a (possibly gzip-compressed) FASTQ file.
pub fn for_each_record(
    path: &Path,
    mut visit: impl FnMut(&FastqRecord) -> Result<()>,
) -> Result<()> {
    let mut reader = compress::open_maybe_gzip(path)?;
    let mut record = FastqRecord::default();
    let mut line = String::new();
    let mut index = 0usize;
    loop {
        if !read_line(&mut *reader, &mut line, path)? {
            return Ok(());
        }
        if line.is_empty() {
            continue;
        }
        index += 1;
        let malformed = |what: &str| Error::parse(path, format!("record {index}: {what}"));
        record.description = line
            .strip_prefix('@')
            .ok_or_else(|| malformed("header does not start with '@'"))?
            .to_string();
        if !read_line(&mut *reader, &mut record.sequence, path)? {
            return Err(malformed("missing sequence line"));
        }
        if !read_line(&mut *reader, &mut line, path)? || !line.starts_with('+') {
            return Err(malformed("missing '+' separator line"));
        }
        if !read_line(&mut *reader, &mut record.quality, path)? {
            return Err(malformed("missing quality line"));
        }
        if record.quality.len() != record.sequence.len() {
            return Err(malformed("sequence and quality lengths differ"));
        }
        visit(&record)?;
    }
}

/// Write records to an uncompressed FASTQ file.
pub fn write_fastq<'a>(
    path: &Path,
    records: impl IntoIterator<Item = &'a FastqRecord>,
) -> Result<()> {
    let mut out = compress::create(path, false)?;
    for record in records {
        record.write_to(&mut *out).map_err(|e| Error::io(path, e))?;
    }
    out.flush().map_err(|e| Error::io(path, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.fastq");
        let records = vec![
            FastqRecord {
                description: "r1 extra".into(),
                sequence: "ACGT".into(),
                quality: "IIII".into(),
            },
            FastqRecord {
                description: "r2".into(),
                sequence: "GG".into(),
                quality: "II".into(),
            },
        ];
        write_fastq(&path, &records).unwrap();
        let mut seen = Vec::new();
        for_each_record(&path, |r| {
            seen.push(r.clone());
            Ok(())
        })
        .unwrap();
        assert_eq!(seen, records);
        assert_eq!(seen[0].identifier(), "r1");
    }
}
