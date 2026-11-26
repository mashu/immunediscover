//! Appending the command line to `immunediscover.log`.

use std::fs::OpenOptions;
use std::io::Write;

/// File in the working directory that records every invocation.
pub const LOG_FILE: &str = "immunediscover.log";

/// Record the invocation; failures are ignored (the log is a convenience).
pub fn record(version: &str, args: &[String]) {
    let Ok(mut file) = OpenOptions::new().create(true).append(true).open(LOG_FILE) else {
        return;
    };
    let timestamp = chrono::Local::now().format("%Y-%m-%dT%H:%M:%S%.3f");
    let quoted: Vec<String> = args.iter().map(|a| format!("{a:?}")).collect();
    let _ = writeln!(
        file,
        "[ Info: {version} {timestamp} - Parsing command line arguments: [{}]",
        quoted.join(", ")
    );
}
