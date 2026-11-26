//! Work-directory layout for caches and intermediate files.

use std::path::{Component, Path, PathBuf};

use md5::{Digest, Md5};

use crate::error::{Error, Result};
use crate::io::path::stem_before_first_dot;

/// Default work directory, relative to the current directory.
pub const DEFAULT_WORK_DIR: &str = ".immunediscover";

/// Absolute, lexically normalised form of `path` (no symlink resolution).
pub fn absolute(path: &Path) -> Result<PathBuf> {
    let joined = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let mut normalised = PathBuf::new();
    for component in joined.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalised.pop();
            }
            other => normalised.push(other),
        }
    }
    Ok(normalised)
}

/// Resolve the `--work-dir` argument; an empty value means the default.
pub fn resolve_work_dir(argument: &str) -> Result<PathBuf> {
    let trimmed = argument.trim();
    if trimmed.is_empty() {
        return absolute(Path::new(DEFAULT_WORK_DIR));
    }
    let expanded = match trimmed.strip_prefix("~/") {
        Some(rest) => match std::env::var_os("HOME") {
            Some(home) => PathBuf::from(home).join(rest),
            None => PathBuf::from(trimmed),
        },
        None => PathBuf::from(trimmed),
    };
    absolute(&expanded)
}

/// Files and directories of one `discover blast` run.
#[derive(Debug, Clone)]
pub struct WorkLayout {
    pub work_dir: PathBuf,
    /// Per-input directory holding the BLAST cache and query FASTA.
    pub run_dir: PathBuf,
    pub query_fasta: PathBuf,
    pub hits: PathBuf,
    pub affixes: PathBuf,
    pub combined_fasta: PathBuf,
    pub extended_fasta: PathBuf,
}

impl WorkLayout {
    /// Layout for a demultiplexed input table and a database FASTA.
    ///
    /// The run directory is keyed by the MD5 of the input's absolute path;
    /// database-derived files are named after the FASTA file name up to its
    /// first dot. The paths only say where a cache lives; whether it is
    /// still valid is decided by its [`Fingerprint`].
    pub fn new(work_dir: &Path, input_tsv: &Path, database_fasta: &Path) -> Result<Self> {
        let absolute_input = absolute(input_tsv)?;
        let digest = Md5::digest(absolute_input.to_string_lossy().as_bytes());
        let tag: String = digest.iter().take(8).map(|b| format!("{b:02x}")).collect();
        let run_dir = work_dir.join("blast").join(tag);
        let stem = stem_before_first_dot(database_fasta);
        Ok(WorkLayout {
            work_dir: work_dir.to_path_buf(),
            query_fasta: run_dir.join("query.fasta"),
            hits: run_dir.join("hits.blast.gz"),
            run_dir,
            affixes: work_dir.join(format!("{stem}.affixes")),
            combined_fasta: work_dir.join(format!("{stem}-combined.fasta")),
            extended_fasta: work_dir.join(format!("{stem}-combined-extended.fasta")),
        })
    }

    pub fn create_directories(&self) -> Result<()> {
        std::fs::create_dir_all(&self.run_dir).map_err(|e| Error::io(&self.run_dir, e))
    }
}

/// Digest of everything a cached file depends on.
///
/// A cached file is reused only while the fingerprint stored next to it
/// (`<file>.key`) equals the fingerprint of the current inputs and
/// parameters; otherwise it is rebuilt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint(String);

impl Fingerprint {
    /// Fingerprint of a list of parts (order matters).
    pub fn of<S: AsRef<str>>(parts: &[S]) -> Self {
        let mut hasher = Md5::new();
        for part in parts {
            let part = part.as_ref();
            hasher.update((part.len() as u64).to_le_bytes());
            hasher.update(part.as_bytes());
        }
        Fingerprint(
            hasher
                .finalize()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect(),
        )
    }

    fn key_path(cached: &Path) -> PathBuf {
        let mut name = cached.as_os_str().to_os_string();
        name.push(".key");
        PathBuf::from(name)
    }

    /// Does `cached` exist and was it built from these inputs?
    pub fn matches(&self, cached: &Path) -> bool {
        cached.is_file()
            && std::fs::read_to_string(Self::key_path(cached)).is_ok_and(|key| key.trim() == self.0)
    }

    /// Record that `cached` was built from these inputs.
    pub fn record(&self, cached: &Path) -> Result<()> {
        let path = Self::key_path(cached);
        std::fs::write(&path, &self.0).map_err(|e| Error::io(&path, e))
    }
}

/// Identity of an input file that is too large to hash on every run:
/// absolute path, size and modification time.
pub fn file_identity(path: &Path) -> Result<String> {
    let metadata = std::fs::metadata(path).map_err(|e| Error::io(path, e))?;
    let modified = metadata
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_nanos());
    Ok(format!(
        "{}|{}|{modified}",
        absolute(path)?.display(),
        metadata.len()
    ))
}

/// MD5 of a file's content, as hex text.
pub fn content_digest(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path).map_err(|e| Error::io(path, e))?;
    Ok(Md5::digest(&bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect())
}

/// `dir/<name up to first dot>.<extension>` next to `path`.
pub fn replace_extension(path: &Path, extension: &str) -> PathBuf {
    let file = format!("{}.{extension}", stem_before_first_dot(path));
    match path.parent() {
        Some(parent) => parent.join(file),
        None => PathBuf::from(file),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_is_keyed_by_input_path_and_database_stem() {
        let layout = WorkLayout::new(
            Path::new("/w"),
            Path::new("/data/a/../reads.tsv.gz"),
            Path::new("/db/IGHV.v2.fasta"),
        )
        .unwrap();
        // md5("/data/reads.tsv.gz") = 0d6b4d2b0aa5a0ba...
        assert!(layout.run_dir.starts_with("/w/blast"));
        assert_eq!(layout.run_dir.file_name().unwrap().len(), 16);
        assert_eq!(layout.affixes, PathBuf::from("/w/IGHV.affixes"));
        assert_eq!(
            layout.extended_fasta,
            PathBuf::from("/w/IGHV-combined-extended.fasta")
        );
        let same = WorkLayout::new(
            Path::new("/w"),
            Path::new("/data/reads.tsv.gz"),
            Path::new("x.fa"),
        )
        .unwrap();
        assert_eq!(layout.run_dir, same.run_dir);
    }

    #[test]
    fn cache_is_valid_only_for_the_inputs_it_was_built_from() {
        let dir = tempfile::tempdir().unwrap();
        let cached = dir.path().join("hits.blast.gz");
        let key = Fingerprint::of(&["reads|10|1", "db-digest", "-task megablast"]);
        assert!(!key.matches(&cached), "nothing cached yet");
        std::fs::write(&cached, b"x").unwrap();
        assert!(!key.matches(&cached), "a file without a key is not trusted");
        key.record(&cached).unwrap();
        assert!(key.matches(&cached));
        let other = Fingerprint::of(&["reads|10|1", "db-digest", "-task blastn"]);
        assert!(!other.matches(&cached));
        // Parts are length-prefixed: moving a boundary changes the key.
        assert_ne!(Fingerprint::of(&["ab", "c"]), Fingerprint::of(&["a", "bc"]));
    }

    #[test]
    fn file_identity_changes_with_content_size() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("reads.tsv");
        std::fs::write(&path, b"abc").unwrap();
        let before = file_identity(&path).unwrap();
        std::fs::write(&path, b"abcd").unwrap();
        assert_ne!(before, file_identity(&path).unwrap());
        assert_eq!(content_digest(&path).unwrap().len(), 32);
    }

    #[test]
    fn replace_extension_cuts_at_first_dot() {
        assert_eq!(
            replace_extension(Path::new("/d/reads.tsv.gz"), "blast"),
            PathBuf::from("/d/reads.blast")
        );
    }
}
