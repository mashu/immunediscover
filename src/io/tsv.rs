//! Delimited text tables, read and written the way CSV.jl 0.10 does it.
//!
//! Reading: gzip is detected from the content, empty unquoted fields are
//! missing, `""` is an empty string, column types are inferred, short rows
//! are padded with missing values, extra fields open new `ColumnN` columns
//! and blank lines are skipped. Low-cardinality text columns are marked as
//! pooled (see [`crate::frame::Pool`]).
//!
//! Writing: missing values and empty strings are both empty fields, fields
//! containing the delimiter, a quote or a line break are quoted, floats use
//! Julia's formatting.

use std::io::{Read, Write};
use std::path::Path;

use memchr::{memchr, memchr2};

use super::compress;
use crate::error::{Error, Result};
use crate::frame::{format_float, infer_column, Column, Frame, Pool, Value};

/// CSV.jl's default pooling rule: at most this many distinct values ...
const POOL_MAX_LEVELS: usize = 500;
/// ... making up at most this fraction of the rows.
const POOL_MAX_FRACTION: f64 = 0.2;

/// Options for [`read_table`].
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// Field separator; `None` picks tab or comma from the header line.
    pub delimiter: Option<u8>,
    /// Columns that must stay text regardless of their content.
    pub text_columns: Vec<String>,
    /// Maximum number of data rows to read.
    pub limit: Option<usize>,
    /// Column names for files without a header line.
    pub header: Option<Vec<String>>,
}

impl ReadOptions {
    /// Tab-separated input.
    pub fn tsv() -> Self {
        ReadOptions {
            delimiter: Some(b'\t'),
            ..Self::default()
        }
    }

    pub fn text_column(mut self, name: &str) -> Self {
        self.text_columns.push(name.to_string());
        self
    }

    pub fn limit(mut self, limit: Option<usize>) -> Self {
        self.limit = limit;
        self
    }

    /// The file has no header line; use these column names.
    pub fn header<S: AsRef<str>>(mut self, names: &[S]) -> Self {
        self.header = Some(names.iter().map(|n| n.as_ref().to_string()).collect());
        self
    }
}

/// Read a tab-separated table with inferred column types.
pub fn read_tsv(path: &Path) -> Result<Frame> {
    read_table(path, &ReadOptions::tsv())
}

/// Read a delimited table.
pub fn read_table(path: &Path, options: &ReadOptions) -> Result<Frame> {
    let mut bytes = Vec::new();
    compress::open_maybe_gzip(path)?
        .read_to_end(&mut bytes)
        .map_err(|e| Error::io(path, e))?;
    parse_table(&bytes, options).map_err(|message| Error::parse(path, message))
}

fn parse_table(bytes: &[u8], options: &ReadOptions) -> std::result::Result<Frame, String> {
    let delimiter = options.delimiter.unwrap_or_else(|| detect_delimiter(bytes));
    let mut parser = RecordParser::new(bytes, delimiter);

    let names = match &options.header {
        Some(names) => names.clone(),
        None => match parser.next_record()? {
            Some(header) => unique_names(header),
            None => return Ok(Frame::new()),
        },
    };
    let mut names = names;
    let mut cells: Vec<Vec<Option<String>>> = vec![Vec::new(); names.len()];
    let mut ragged = 0usize;
    let mut nrows = 0usize;

    while options.limit.is_none_or(|limit| nrows < limit) {
        let Some(record) = parser.next_record()? else {
            break;
        };
        if record.len() != names.len() {
            ragged += 1;
        }
        // Extra fields widen the table, as CSV.jl does.
        while record.len() > names.len() {
            names.push(format!("Column{}", names.len() + 1));
            cells.push(vec![None; nrows]);
        }
        let mut record = record.into_iter();
        for column in cells.iter_mut() {
            column.push(record.next().flatten());
        }
        nrows += 1;
    }
    if ragged > 0 {
        crate::warn!(
            "{ragged} rows did not have the expected number of fields; missing fields were set to missing, extra fields opened new columns"
        );
    }

    let mut frame = Frame::new();
    for (name, column) in names.into_iter().zip(cells) {
        let force_text = options.text_columns.contains(&name);
        let column = infer_column(column, force_text);
        let pool = pool_of(&column);
        frame
            .push_column(name.clone(), column)
            .map_err(|e| e.to_string())?;
        if let Some(pool) = pool {
            frame.set_pool(&name, pool).map_err(|e| e.to_string())?;
        }
    }
    Ok(frame)
}

/// Levels of a text column if CSV.jl would store it as a pooled array:
/// no more than 500 distinct values (missing counts as one) and
/// `(levels - 1) / rows <= 0.2`.
fn pool_of(column: &Column) -> Option<Pool> {
    let cells = column.as_strs()?;
    let mut seen: std::collections::HashSet<Option<&str>> = std::collections::HashSet::new();
    let mut levels: Vec<Option<String>> = Vec::new();
    for cell in cells {
        if seen.insert(cell.as_deref()) {
            levels.push(cell.clone());
            if levels.len() > POOL_MAX_LEVELS {
                return None;
            }
        }
    }
    let fraction = (levels.len() as f64 - 1.0) / cells.len() as f64;
    (fraction <= POOL_MAX_FRACTION).then(|| std::sync::Arc::new(levels))
}

fn detect_delimiter(bytes: &[u8]) -> u8 {
    let first_line = &bytes[..memchr(b'\n', bytes).unwrap_or(bytes.len())];
    let tabs = first_line.iter().filter(|&&b| b == b'\t').count();
    let commas = first_line.iter().filter(|&&b| b == b',').count();
    if tabs >= commas && tabs > 0 {
        b'\t'
    } else {
        b','
    }
}

/// Header cells to column names: blanks become `ColumnN`, repeats get `_1`, `_2`.
fn unique_names(header: Vec<Option<String>>) -> Vec<String> {
    let mut names: Vec<String> = Vec::with_capacity(header.len());
    for (i, cell) in header.into_iter().enumerate() {
        let base = cell
            .filter(|c| !c.is_empty())
            .unwrap_or_else(|| format!("Column{}", i + 1));
        let mut name = base.clone();
        let mut suffix = 0;
        while names.contains(&name) {
            suffix += 1;
            name = format!("{base}_{suffix}");
        }
        names.push(name);
    }
    names
}

/// Splits a byte buffer into records of optional fields.
struct RecordParser<'a> {
    bytes: &'a [u8],
    pos: usize,
    delimiter: u8,
}

impl<'a> RecordParser<'a> {
    fn new(bytes: &'a [u8], delimiter: u8) -> Self {
        // Skip a UTF-8 byte order mark.
        let pos = if bytes.starts_with(&[0xEF, 0xBB, 0xBF]) {
            3
        } else {
            0
        };
        RecordParser {
            bytes,
            pos,
            delimiter,
        }
    }

    fn skip_blank_lines(&mut self) {
        loop {
            match self.bytes.get(self.pos) {
                Some(b'\n') => self.pos += 1,
                Some(b'\r') if self.bytes.get(self.pos + 1) == Some(&b'\n') => self.pos += 2,
                _ => break,
            }
        }
    }

    fn next_record(&mut self) -> std::result::Result<Option<Vec<Option<String>>>, String> {
        self.skip_blank_lines();
        if self.pos >= self.bytes.len() {
            return Ok(None);
        }
        let mut fields = Vec::new();
        loop {
            let (field, end_of_record) = if self.bytes.get(self.pos) == Some(&b'"') {
                self.quoted_field()?
            } else {
                self.bare_field()?
            };
            fields.push(field);
            if end_of_record {
                return Ok(Some(fields));
            }
        }
    }

    /// Unquoted field: runs to the next delimiter or line break.
    fn bare_field(&mut self) -> std::result::Result<(Option<String>, bool), String> {
        let rest = &self.bytes[self.pos..];
        let stop = memchr2(self.delimiter, b'\n', rest).unwrap_or(rest.len());
        let at_delimiter = rest.get(stop) == Some(&self.delimiter);
        let mut raw = &rest[..stop];
        if !at_delimiter && raw.last() == Some(&b'\r') {
            raw = &raw[..raw.len() - 1];
        }
        self.pos += (stop + 1).min(rest.len());
        let field = if raw.is_empty() {
            None
        } else {
            Some(to_text(raw)?)
        };
        Ok((field, !at_delimiter))
    }

    /// Quoted field: `""` is an escaped quote; may span lines.
    fn quoted_field(&mut self) -> std::result::Result<(Option<String>, bool), String> {
        let mut value = Vec::new();
        let mut i = self.pos + 1;
        loop {
            let Some(offset) = memchr(b'"', &self.bytes[i..]) else {
                return Err("unterminated quoted field".to_string());
            };
            value.extend_from_slice(&self.bytes[i..i + offset]);
            i += offset + 1;
            if self.bytes.get(i) == Some(&b'"') {
                value.push(b'"');
                i += 1;
            } else {
                break;
            }
        }
        // Skip to the field terminator.
        let rest = &self.bytes[i..];
        let stop = memchr2(self.delimiter, b'\n', rest).unwrap_or(rest.len());
        let at_delimiter = rest.get(stop) == Some(&self.delimiter);
        self.pos = i + (stop + 1).min(rest.len());
        Ok((Some(to_text(&value)?), !at_delimiter))
    }
}

fn to_text(raw: &[u8]) -> std::result::Result<String, String> {
    std::str::from_utf8(raw)
        .map(str::to_owned)
        .map_err(|_| "invalid UTF-8 in table".to_string())
}

/// Options for [`write_table`].
#[derive(Debug, Clone, Copy)]
pub struct WriteOptions {
    pub delimiter: u8,
    pub gzip: bool,
}

impl WriteOptions {
    pub fn tsv() -> Self {
        WriteOptions {
            delimiter: b'\t',
            gzip: false,
        }
    }

    pub fn tsv_gz() -> Self {
        WriteOptions {
            delimiter: b'\t',
            gzip: true,
        }
    }

    pub fn csv() -> Self {
        WriteOptions {
            delimiter: b',',
            gzip: false,
        }
    }
}

/// Write `frame` to `path`.
pub fn write_table(path: &Path, frame: &Frame, options: WriteOptions) -> Result<()> {
    let mut out = compress::create(path, options.gzip)?;
    write_frame(&mut out, frame, options.delimiter).map_err(|e| Error::io(path, e))?;
    out.flush().map_err(|e| Error::io(path, e))
}

/// Plain tab-separated output.
pub fn write_tsv(path: &Path, frame: &Frame) -> Result<()> {
    write_table(path, frame, WriteOptions::tsv())
}

/// Gzip-compressed tab-separated output.
pub fn write_tsv_gz(path: &Path, frame: &Frame) -> Result<()> {
    write_table(path, frame, WriteOptions::tsv_gz())
}

fn write_frame(out: &mut dyn Write, frame: &Frame, delimiter: u8) -> std::io::Result<()> {
    let mut line: Vec<u8> = Vec::with_capacity(256);
    for (i, name) in frame.names().iter().enumerate() {
        if i > 0 {
            line.push(delimiter);
        }
        push_text(&mut line, name, delimiter);
    }
    line.push(b'\n');
    out.write_all(&line)?;

    let columns: &[Column] = frame.columns();
    for row in 0..frame.nrows() {
        line.clear();
        for (i, column) in columns.iter().enumerate() {
            if i > 0 {
                line.push(delimiter);
            }
            match column.get(row) {
                Value::Missing => {}
                Value::Int(v) => line.extend_from_slice(v.to_string().as_bytes()),
                Value::Float(v) => line.extend_from_slice(format_float(v).as_bytes()),
                Value::Bool(v) => line.extend_from_slice(if v { b"true" } else { b"false" }),
                Value::Str(s) => push_text(&mut line, s, delimiter),
            }
        }
        line.push(b'\n');
        out.write_all(&line)?;
    }
    Ok(())
}

fn push_text(line: &mut Vec<u8>, text: &str, delimiter: u8) {
    let bytes = text.as_bytes();
    let needs_quotes = bytes
        .iter()
        .any(|&b| b == delimiter || b == b'"' || b == b'\n' || b == b'\r');
    if !needs_quotes {
        line.extend_from_slice(bytes);
        return;
    }
    line.push(b'"');
    for &b in bytes {
        if b == b'"' {
            line.push(b'"');
        }
        line.push(b);
    }
    line.push(b'"');
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> Frame {
        parse_table(text.as_bytes(), &ReadOptions::tsv()).unwrap()
    }

    fn render(frame: &Frame) -> String {
        let mut out = Vec::new();
        write_frame(&mut out, frame, b'\t').unwrap();
        String::from_utf8(out).unwrap()
    }

    // Reference behaviour: CSV.jl 0.10.15 on the same text.
    #[test]
    fn quoting_and_empty_fields_match_csv_jl() {
        let input = "a\tb\tc\td\te\n\"\"\t\t1.0\ttrue\t\"a\"\"b\"\nx\t\"\"\tNaN\tfalse\t\"c\td\"\ny z\tq\t2.5\ttrue\t\"e\nf\"\n";
        let frame = parse(input);
        assert_eq!(frame.strs("a").unwrap()[0].as_deref(), Some(""));
        assert_eq!(frame.strs("b").unwrap()[0], None);
        assert_eq!(frame.strs("e").unwrap()[2].as_deref(), Some("e\nf"));
        // Empty strings and missing values are both written as empty fields.
        let expected = "a\tb\tc\td\te\n\t\t1.0\ttrue\t\"a\"\"b\"\nx\t\tNaN\tfalse\t\"c\td\"\ny z\tq\t2.5\ttrue\t\"e\nf\"\n";
        assert_eq!(render(&frame), expected);
    }

    #[test]
    fn ragged_rows_blank_lines_and_duplicate_headers() {
        let frame = parse("a\ta\tb\n1\t2\n3\t4\t5\t6\n\n7\t8\t9\r\n");
        assert_eq!(frame.names(), ["a", "a_1", "b", "Column4"]);
        assert_eq!(frame.nrows(), 3);
        assert_eq!(
            frame.column("b").unwrap(),
            &Column::Int(vec![None, Some(5), Some(9)])
        );
        assert_eq!(
            frame.column("Column4").unwrap(),
            &Column::Int(vec![None, Some(6), None])
        );
    }

    #[test]
    fn low_cardinality_text_columns_are_pooled() {
        let mut text = String::from("s\tn\tu\n");
        for i in 1..=30 {
            text.push_str(&format!(
                "{}\t{}\tu{}\n",
                if i % 3 == 0 { "zz" } else { "aa" },
                i % 4,
                i
            ));
        }
        let frame = parse(&text);
        let pool = frame
            .pool_at(0)
            .expect("two levels in thirty rows are pooled");
        assert_eq!(
            pool.as_slice(),
            [Some("aa".to_string()), Some("zz".to_string())]
        );
        assert!(
            frame.pool_at(1).is_none(),
            "numeric columns are never pooled"
        );
        assert!(
            frame.pool_at(2).is_none(),
            "all-distinct text is not pooled"
        );
        // Two levels in six rows: (2 - 1) / 6 <= 0.2, pooled; in two rows not.
        assert!(parse("s\na\nb\na\nb\na\nb\n").pool_at(0).is_some());
        assert!(parse("s\na\nb\n").pool_at(0).is_none());
    }

    #[test]
    fn delimiter_is_detected_and_limit_applies() {
        let frame =
            parse_table(b"x,y\n1,2\n3,4\n", &ReadOptions::default().limit(Some(1))).unwrap();
        assert_eq!(frame.names(), ["x", "y"]);
        assert_eq!(frame.nrows(), 1);
    }

    #[test]
    fn text_columns_are_not_inferred() {
        let options = ReadOptions::tsv().text_column("case");
        let frame = parse_table(b"well\tcase\n1\t007\n", &options).unwrap();
        assert_eq!(frame.strs("case").unwrap()[0].as_deref(), Some("007"));
        assert_eq!(frame.column("well").unwrap(), &Column::from_ints([1]));
    }
}
