//! Crate-wide error type.

use std::path::PathBuf;

/// Convenience alias used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Every failure the library can report.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("{path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("I/O error: {0}")]
    BareIo(#[from] std::io::Error),

    #[error("{path}: {message}")]
    Parse { path: PathBuf, message: String },

    #[error("Column '{column}' not found. Available: {available}")]
    MissingColumn { column: String, available: String },

    #[error("Column '{column}' has type {found}, expected {expected}")]
    ColumnType {
        column: String,
        expected: &'static str,
        found: &'static str,
    },

    #[error("Invalid regular expression '{pattern}': {message}")]
    Regex { pattern: String, message: String },

    #[error("{tool}: {message}")]
    ExternalTool { tool: &'static str, message: String },

    #[error("{0}")]
    Invalid(String),
}

impl Error {
    /// Free-form validation failure.
    pub fn invalid(message: impl Into<String>) -> Self {
        Error::Invalid(message.into())
    }

    /// Wrap an I/O error with the path it happened on.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Error::Io {
            path: path.into(),
            source,
        }
    }

    /// Malformed file content.
    pub fn parse(path: impl Into<PathBuf>, message: impl Into<String>) -> Self {
        Error::Parse {
            path: path.into(),
            message: message.into(),
        }
    }

    /// Failure of an external program (blastn, makeblastdb, bwa).
    pub fn tool(tool: &'static str, message: impl Into<String>) -> Self {
        Error::ExternalTool {
            tool,
            message: message.into(),
        }
    }
}
