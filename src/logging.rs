//! Minimal stderr logging in the style of Julia's `@info` / `@warn`.
//!
//! Messages go to stderr so that commands printing results to stdout
//! (`fasta hash`, `fasta diff`) stay pipe-friendly.

use std::sync::atomic::{AtomicBool, Ordering};

static QUIET: AtomicBool = AtomicBool::new(false);

/// Silence informational messages (warnings are still printed).
pub fn set_quiet(quiet: bool) {
    QUIET.store(quiet, Ordering::Relaxed);
}

#[doc(hidden)]
pub fn emit_info(message: std::fmt::Arguments<'_>) {
    if !QUIET.load(Ordering::Relaxed) {
        eprintln!("[ Info: {message}");
    }
}

#[doc(hidden)]
pub fn emit_warn(message: std::fmt::Arguments<'_>) {
    eprintln!("[ Warning: {message}");
}

/// Informational message, mirrors Julia's `@info`.
#[macro_export]
macro_rules! info {
    ($($arg:tt)*) => { $crate::logging::emit_info(format_args!($($arg)*)) };
}

/// Warning message, mirrors Julia's `@warn`.
#[macro_export]
macro_rules! warn {
    ($($arg:tt)*) => { $crate::logging::emit_warn(format_args!($($arg)*)) };
}
