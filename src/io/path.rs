//! Output path conventions.

use std::path::{Path, PathBuf};

/// Return `path` with a `.gz` suffix, appending it when absent.
pub fn always_gz(path: &Path) -> PathBuf {
    if path.as_os_str().to_string_lossy().ends_with(".gz") {
        path.to_path_buf()
    } else {
        let mut with_suffix = path.as_os_str().to_os_string();
        with_suffix.push(".gz");
        PathBuf::from(with_suffix)
    }
}

/// File name up to (not including) its first dot, e.g. `IGHV.v2.fasta` -> `IGHV`.
pub fn stem_before_first_dot(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    name.split('.').next().unwrap_or_default().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gz_suffix_is_added_once() {
        assert_eq!(always_gz(Path::new("a.tsv")), PathBuf::from("a.tsv.gz"));
        assert_eq!(always_gz(Path::new("a.tsv.gz")), PathBuf::from("a.tsv.gz"));
    }

    #[test]
    fn stem_stops_at_first_dot() {
        assert_eq!(stem_before_first_dot(Path::new("/x/IGHV.v2.fasta")), "IGHV");
    }
}
