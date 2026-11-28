//! Where a genome is read from: a file, the standard input, or the
//! internet.
//!
//! A genome is read as a stream and never stored, so one that lives on a
//! public server can be searched without keeping a copy: the download is
//! piped through the search. `curl` does the downloading.

use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::process::{Child, ChildStdout, Command, Stdio};

use flate2::read::MultiGzDecoder;

use crate::error::{Error, Result};

/// Human genome assemblies that can be named instead of a file: the names
/// they answer to, what they are, and where UCSC serves them.
pub const KNOWN_GENOMES: [(&[&str], &str, &str); 3] = [
    (
        &["grch37", "hg19"],
        "GRCh37 (hg19)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    ),
    (
        &["grch38", "hg38"],
        "GRCh38 (hg38)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    ),
    (
        &["t2t", "chm13", "hs1"],
        "T2T-CHM13 v2.0 (hs1)",
        "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz",
    ),
];

const BUFFER: usize = 1 << 20;
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// What a genome argument refers to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Origin {
    File(String),
    StandardInput,
    /// An address, with the name of the assembly if it is a known one.
    Internet {
        url: String,
        assembly: Option<&'static str>,
    },
}

impl Origin {
    /// Decide what `argument` means: an existing file, `-` for the
    /// standard input, the name of a known assembly, or an address.
    pub fn of(argument: &str) -> Result<Origin> {
        if argument == "-" {
            return Ok(Origin::StandardInput);
        }
        if Path::new(argument).is_file() {
            return Ok(Origin::File(argument.to_string()));
        }
        let lower = argument.to_ascii_lowercase();
        let known = KNOWN_GENOMES
            .iter()
            .find(|(names, _, _)| names.contains(&lower.as_str()));
        if let Some((_, assembly, url)) = known {
            return Ok(Origin::Internet {
                url: url.to_string(),
                assembly: Some(assembly),
            });
        }
        let schemes = ["https://", "http://", "ftp://", "file://"];
        if schemes.iter().any(|scheme| lower.starts_with(scheme)) {
            return Ok(Origin::Internet {
                url: argument.to_string(),
                assembly: None,
            });
        }
        let names: Vec<&str> = KNOWN_GENOMES.iter().map(|(names, _, _)| names[0]).collect();
        Err(Error::invalid(format!(
            "Genome '{argument}' is not a file, not an address and not one of the known genomes ({})",
            names.join(", ")
        )))
    }

    /// How the genome is called in messages and tables.
    pub fn label(&self) -> String {
        match self {
            Origin::File(path) => path.clone(),
            Origin::StandardInput => "standard input".to_string(),
            Origin::Internet {
                assembly: Some(assembly),
                ..
            } => assembly.to_string(),
            Origin::Internet { url, .. } => url.clone(),
        }
    }
}

/// An opened genome: FASTA text, decompressed if it was gzip.
pub struct Source {
    pub reader: Box<dyn BufRead + Send>,
    download: Option<(Child, String)>,
}

/// Decompress when the content is gzip (told by its first bytes, not by a
/// name).
fn decompressed<R: Read + Send + 'static>(input: R) -> std::io::Result<Box<dyn BufRead + Send>> {
    let mut reader = BufReader::with_capacity(BUFFER, input);
    let is_gzip = reader.fill_buf()?.starts_with(&GZIP_MAGIC);
    Ok(if is_gzip {
        Box::new(BufReader::with_capacity(
            BUFFER,
            MultiGzDecoder::new(reader),
        ))
    } else {
        Box::new(reader)
    })
}

fn download(url: &str) -> Result<(Child, ChildStdout)> {
    let mut child = Command::new("curl")
        .args(["--fail", "--silent", "--show-error", "--location"])
        .args(["--retry", "3", "--connect-timeout", "30"])
        .arg(url)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|e| {
            Error::tool(
                "curl",
                format!(
                    "could not be started ({e}). Reading a genome from the internet needs curl; without it, download the file and give its path"
                ),
            )
        })?;
    let output = child
        .stdout
        .take()
        .ok_or_else(|| Error::tool("curl", "gave no output to read"))?;
    Ok((child, output))
}

impl Source {
    pub fn open(origin: &Origin) -> Result<Source> {
        let label = origin.label();
        let failed = |e: std::io::Error| Error::io(&label, e);
        match origin {
            Origin::File(path) => {
                let file = std::fs::File::open(path).map_err(failed)?;
                Ok(Source {
                    reader: decompressed(file).map_err(failed)?,
                    download: None,
                })
            }
            Origin::StandardInput => Ok(Source {
                reader: decompressed(std::io::stdin()).map_err(failed)?,
                download: None,
            }),
            Origin::Internet { url, .. } => {
                let (child, output) = download(url)?;
                // Nothing arriving is a failed download; `finish` tells why.
                let reader = decompressed(output).map_err(failed)?;
                Ok(Source {
                    reader,
                    download: Some((child, url.clone())),
                })
            }
        }
    }

    /// Close the genome after it was read to its end. A download that
    /// broke off is an error: a genome read in part would give wrong
    /// answers.
    pub fn finish(self) -> Result<()> {
        drop(self.reader);
        let Some((mut child, url)) = self.download else {
            return Ok(());
        };
        let status = child
            .wait()
            .map_err(|e| Error::tool("curl", e.to_string()))?;
        if status.success() {
            Ok(())
        } else {
            Err(Error::tool(
                "curl",
                format!("the download of {url} failed ({status})"),
            ))
        }
    }

    /// Close a genome that could not be read to its end. If its download
    /// had already failed, that failure is returned: it is the reason.
    pub fn abandon(self) -> Option<Error> {
        let (mut child, url) = self.download?;
        let ended = child.try_wait().ok().flatten();
        let _ = child.kill();
        let _ = child.wait();
        drop(self.reader);
        ended
            .filter(|status| !status.success())
            .map(|status| Error::tool("curl", format!("the download of {url} failed ({status})")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn names_files_and_addresses_are_told_apart() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("genome.fa");
        std::fs::write(&file, ">chr1\nACGT\n").unwrap();
        let path = file.to_string_lossy().into_owned();
        assert_eq!(Origin::of(&path).unwrap(), Origin::File(path));
        assert_eq!(Origin::of("-").unwrap(), Origin::StandardInput);
        for (name, assembly) in [("GRCh38", "GRCh38 (hg38)"), ("hg19", "GRCh37 (hg19)")] {
            let origin = Origin::of(name).unwrap();
            assert_eq!(origin.label(), assembly);
            assert!(matches!(origin, Origin::Internet { url, .. } if url.ends_with(".fa.gz")));
        }
        assert_eq!(Origin::of("t2t").unwrap().label(), "T2T-CHM13 v2.0 (hs1)");
        let address = "https://example.org/genome.fa.gz";
        assert_eq!(Origin::of(address).unwrap().label(), address);
        let error = Origin::of("no-such-genome.fa").unwrap_err().to_string();
        assert!(error.contains("grch37, grch38, t2t"), "{error}");
    }

    #[test]
    fn a_file_is_read_plain_or_gzipped_whatever_its_name() {
        let dir = tempfile::tempdir().unwrap();
        let plain = dir.path().join("plain.txt");
        std::fs::write(&plain, ">chr1\nACGT\n").unwrap();
        let zipped = dir.path().join("zipped.txt");
        let mut encoder = flate2::write::GzEncoder::new(
            std::fs::File::create(&zipped).unwrap(),
            flate2::Compression::fast(),
        );
        encoder.write_all(b">chr1\nACGT\n").unwrap();
        encoder.finish().unwrap();
        for path in [plain, zipped] {
            let origin = Origin::of(&path.to_string_lossy()).unwrap();
            let mut source = Source::open(&origin).unwrap();
            let mut text = String::new();
            source.reader.read_to_string(&mut text).unwrap();
            assert_eq!(text, ">chr1\nACGT\n");
            source.finish().unwrap();
        }
    }
}
