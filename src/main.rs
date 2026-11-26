//! `immunediscover` command-line entry point.

use std::process::ExitCode;

fn main() -> ExitCode {
    match immunediscover::cli::main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ERROR: {error}");
            ExitCode::FAILURE
        }
    }
}
