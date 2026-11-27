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
    /// Directory of this combination of reads, database and pseudogenes.
    /// Every intermediate file lives here, so that runs on other reads or
    /// with another database never touch each other's files.
    pub run_dir: PathBuf,
    pub query_fasta: PathBuf,
    pub hits: PathBuf,
    pub affixes: PathBuf,
    pub combined_fasta: PathBuf,
    pub extended_fasta: PathBuf,
}

impl WorkLayout {
    /// Layout for a demultiplexed input table, a database FASTA and an
    /// optional pseudogene FASTA.
    ///
    /// The run directory is keyed by the MD5 of the absolute paths of all
    /// three. The paths only say where a cache lives; whether it is still
    /// valid is decided by its [`Fingerprint`].
    pub fn new(
        work_dir: &Path,
        input_tsv: &Path,
        database_fasta: &Path,
        pseudo_fasta: Option<&Path>,
    ) -> Result<Self> {
        let mut hasher = Md5::new();
        for path in [Some(input_tsv), Some(database_fasta), pseudo_fasta] {
            let text = match path {
                Some(path) => absolute(path)?.to_string_lossy().into_owned(),
                None => String::new(),
            };
            hasher.update((text.len() as u64).to_le_bytes());
            hasher.update(text.as_bytes());
        }
        let tag: String = hasher
            .finalize()
            .iter()
            .take(8)
            .map(|b| format!("{b:02x}"))
            .collect();
        let run_dir = work_dir.join("blast").join(tag);
        let stem = stem_before_first_dot(database_fasta);
        Ok(WorkLayout {
            work_dir: work_dir.to_path_buf(),
            query_fasta: run_dir.join("query.fasta"),
            hits: run_dir.join("hits.blast.gz"),
            affixes: run_dir.join(format!("{stem}.affixes")),
            combined_fasta: run_dir.join(format!("{stem}-combined.fasta")),
            extended_fasta: run_dir.join(format!("{stem}-combined-extended.fasta")),
            run_dir,
        })
    }

    pub fn create_directories(&self) -> Result<()> {
        std::fs::create_dir_all(&self.run_dir).map_err(|e| Error::io(&self.run_dir, e))
    }
}

/// Exclusive use of a run directory, for as long as the value lives.
///
/// Two runs on the same reads and database would write the same files.
/// The second one is refused with a message instead of being left to
/// corrupt the first. The lock is a file holding the process number; one
/// left behind by a run that was killed is recognised and replaced.
#[derive(Debug)]
pub struct RunLock {
    path: PathBuf,
}

impl RunLock {
    pub fn acquire(run_dir: &Path) -> Result<Self> {
        let path = run_dir.join("lock");
        for _ in 0..2 {
            let created = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path);
            match created {
                Ok(mut file) => {
                    use std::io::Write as _;
                    write!(file, "{}", std::process::id()).map_err(|e| Error::io(&path, e))?;
                    return Ok(RunLock { path });
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                    let owner = std::fs::read_to_string(&path)
                        .ok()
                        .and_then(|text| text.trim().parse::<u32>().ok());
                    if owner.is_some_and(|pid| !process_is_running(pid)) {
                        // Left behind by a run that did not end normally.
                        let _ = std::fs::remove_file(&path);
                        continue;
                    }
                    return Err(Error::invalid(format!(
                        "Another discover blast run{} is working on the same reads and database \
                         in {}. Wait for it to finish or give this run another --work-dir. \
                         If no such run exists, delete {}.",
                        owner.map_or(String::new(), |pid| format!(" (process {pid})")),
                        run_dir.display(),
                        path.display()
                    )));
                }
                Err(e) => return Err(Error::io(&path, e)),
            }
        }
        Err(Error::invalid(format!(
            "Could not lock {}",
            run_dir.display()
        )))
    }
}

impl Drop for RunLock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// Whether a process with this number exists. Where that cannot be found
/// out, it is taken to exist.
fn process_is_running(pid: u32) -> bool {
    if pid == std::process::id() {
        return true;
    }
    if cfg!(unix) {
        return std::process::Command::new("kill")
            .args(["-0", &pid.to_string()])
            .stderr(std::process::Stdio::null())
            .status()
            .map_or(true, |status| status.success());
    }
    true
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

    fn layout(input: &str, database: &str, pseudo: Option<&str>) -> WorkLayout {
        WorkLayout::new(
            Path::new("/w"),
            Path::new(input),
            Path::new(database),
            pseudo.map(Path::new),
        )
        .unwrap()
    }

    #[test]
    fn every_combination_of_reads_database_and_pseudogenes_has_its_own_directory() {
        let v = layout("/data/a/../reads.tsv.gz", "/db/IGHV.v2.fasta", None);
        assert!(v.run_dir.starts_with("/w/blast"));
        assert_eq!(v.run_dir.file_name().unwrap().len(), 16);
        // Everything derived from the run lives in its directory.
        for file in [&v.hits, &v.query_fasta, &v.affixes, &v.combined_fasta] {
            assert_eq!(file.parent(), Some(v.run_dir.as_path()));
        }
        assert_eq!(
            v.extended_fasta,
            v.run_dir.join("IGHV-combined-extended.fasta")
        );
        // The same files, written differently: the same directory.
        let again = layout("/data/reads.tsv.gz", "/db/x/../IGHV.v2.fasta", None);
        assert_eq!(v.run_dir, again.run_dir);
        // Other database, same-named database elsewhere, other reads,
        // pseudogenes: other directories.
        let others = [
            layout("/data/reads.tsv.gz", "/db/IGHD.fasta", None),
            layout("/data/reads.tsv.gz", "/other/IGHV.v2.fasta", None),
            layout("/data/plate2.tsv.gz", "/db/IGHV.v2.fasta", None),
            layout(
                "/data/reads.tsv.gz",
                "/db/IGHV.v2.fasta",
                Some("/db/pseudo.fasta"),
            ),
        ];
        let mut directories: Vec<&PathBuf> = others.iter().map(|l| &l.run_dir).collect();
        directories.push(&v.run_dir);
        directories.sort();
        directories.dedup();
        assert_eq!(directories.len(), 5);
    }

    #[test]
    fn a_run_directory_is_used_by_one_run_at_a_time() {
        let dir = tempfile::tempdir().unwrap();
        let lock = RunLock::acquire(dir.path()).unwrap();
        let refused = RunLock::acquire(dir.path()).unwrap_err().to_string();
        assert!(refused.contains("Another discover blast run"), "{refused}");
        drop(lock);
        assert!(!dir.path().join("lock").exists());
        // A lock left by a process that no longer exists is replaced.
        std::fs::write(dir.path().join("lock"), "4294967294").unwrap();
        let _lock = RunLock::acquire(dir.path()).unwrap();
        assert_eq!(
            std::fs::read_to_string(dir.path().join("lock")).unwrap(),
            std::process::id().to_string()
        );
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
