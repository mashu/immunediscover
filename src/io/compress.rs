//! Transparent gzip handling for readers and writers.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::error::{Error, Result};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const BUFFER: usize = 1 << 16;

/// Open `path` for reading, decompressing when the content is gzip
/// (detected from the magic bytes, not the file name).
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).map_err(|e| Error::io(path, e))?;
    let mut reader = BufReader::with_capacity(BUFFER, file);
    let is_gzip = reader
        .fill_buf()
        .map_err(|e| Error::io(path, e))?
        .starts_with(&GZIP_MAGIC);
    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(
            BUFFER,
            MultiGzDecoder::new(reader),
        )))
    } else {
        Ok(Box::new(reader))
    }
}

/// Read the whole (possibly gzip-compressed) file as text.
pub fn read_to_string(path: &Path) -> Result<String> {
    let mut text = String::new();
    open_maybe_gzip(path)?
        .read_to_string(&mut text)
        .map_err(|e| Error::io(path, e))?;
    Ok(text)
}

/// Create `path` for writing, gzip-compressed when `gzip` is set.
pub fn create(path: &Path, gzip: bool) -> Result<Box<dyn Write + Send>> {
    let file = File::create(path).map_err(|e| Error::io(path, e))?;
    let writer = BufWriter::with_capacity(BUFFER, file);
    if gzip {
        Ok(Box::new(GzEncoder::new(writer, Compression::default())))
    } else {
        Ok(Box::new(writer))
    }
}
