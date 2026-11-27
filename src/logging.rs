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

/// `1234567` as `1,234,567`.
pub fn with_thousands(number: usize) -> String {
    let digits = number.to_string();
    let mut text = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, digit) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            text.push(',');
        }
        text.push(digit);
    }
    text
}

/// A span of time as people say it: `8.4 s`, `2 min 41 s`, `1 h 05 min`.
pub fn spoken_duration(seconds: f64) -> String {
    let whole = seconds.round() as u64;
    if seconds < 60.0 {
        format!("{seconds:.1} s")
    } else if whole < 3600 {
        format!("{} min {:02} s", whole / 60, whole % 60)
    } else {
        format!("{} h {:02} min", whole / 3600, whole % 3600 / 60)
    }
}

/// Progress of a long step, shown on stderr.
///
/// On a terminal one line is redrawn a few times per second, with a bar,
/// the speed and the time left. When stderr goes to a file, a plain line
/// is written every ten seconds instead. [`Progress::finish`] always
/// writes the total time and the average speed.
///
/// [`Progress::add`] may be called from many threads; it costs one atomic
/// addition unless a line is due.
pub struct Progress {
    doing: &'static str,
    done_as: &'static str,
    unit: &'static str,
    total: usize,
    done: std::sync::atomic::AtomicUsize,
    /// Milliseconds since `started` at which the last line was written.
    last_line: std::sync::atomic::AtomicU64,
    started: std::time::Instant,
    on_terminal: bool,
    silent: bool,
}

const BAR_WIDTH: usize = 24;

impl Progress {
    /// Start counting `total` items of `unit` ("reads"). `doing` and
    /// `done_as` name the step while it runs and afterwards ("Aligning",
    /// "Aligned").
    pub fn new(
        doing: &'static str,
        done_as: &'static str,
        unit: &'static str,
        total: usize,
    ) -> Self {
        use std::io::IsTerminal;
        Progress {
            doing,
            done_as,
            unit,
            total,
            done: std::sync::atomic::AtomicUsize::new(0),
            last_line: std::sync::atomic::AtomicU64::new(0),
            started: std::time::Instant::now(),
            on_terminal: std::io::stderr().is_terminal(),
            silent: QUIET.load(Ordering::Relaxed),
        }
    }

    /// Milliseconds between two lines.
    fn interval(&self) -> u64 {
        if self.on_terminal {
            200
        } else {
            10_000
        }
    }

    /// Count `count` more items as done.
    pub fn add(&self, count: usize) {
        let done = self.done.fetch_add(count, Ordering::Relaxed) + count;
        if self.silent {
            return;
        }
        let now = self.started.elapsed().as_millis() as u64;
        let last = self.last_line.load(Ordering::Relaxed);
        if now < last + self.interval() {
            return;
        }
        // Only the thread that moves the clock forward writes.
        let moved =
            self.last_line
                .compare_exchange(last, now, Ordering::Relaxed, Ordering::Relaxed);
        if moved.is_ok() {
            self.write_line(done, now as f64 / 1000.0);
        }
    }

    fn line(&self, done: usize, seconds: f64) -> String {
        let done = done.min(self.total);
        let share = if self.total == 0 {
            1.0
        } else {
            done as f64 / self.total as f64
        };
        let rate = if seconds > 0.0 {
            done as f64 / seconds
        } else {
            0.0
        };
        let left = if rate > 0.0 {
            format!(
                ", about {} left",
                spoken_duration((self.total - done) as f64 / rate)
            )
        } else {
            String::new()
        };
        let counts = format!(
            "{:>3.0}%  {} of {} {}, {} {} per second{left}",
            100.0 * share,
            with_thousands(done),
            with_thousands(self.total),
            self.unit,
            with_thousands(rate.round() as usize),
            self.unit,
        );
        if self.on_terminal {
            let filled = (share * BAR_WIDTH as f64).round() as usize;
            format!(
                "{} [{}{}] {counts}",
                self.doing,
                "#".repeat(filled),
                "-".repeat(BAR_WIDTH - filled)
            )
        } else {
            format!("{}: {counts}", self.doing)
        }
    }

    fn write_line(&self, done: usize, seconds: f64) {
        let line = self.line(done, seconds);
        if self.on_terminal {
            // Back to the start of the line, write, erase what is left of
            // the line before.
            eprint!("\r[ Info: {line}\x1b[K");
        } else {
            eprintln!("[ Info: {line}");
        }
    }

    /// End the step: the total time and the average speed are written.
    pub fn finish(self) {
        let seconds = self.started.elapsed().as_secs_f64();
        if self.silent {
            return;
        }
        if self.on_terminal {
            eprint!("\r\x1b[K");
        }
        let done = self.done.load(Ordering::Relaxed);
        let rate = if seconds > 0.0 {
            done as f64 / seconds
        } else {
            done as f64
        };
        eprintln!(
            "[ Info: {} {} {} in {} ({} {} per second)",
            self.done_as,
            with_thousands(done),
            self.unit,
            spoken_duration(seconds),
            with_thousands(rate.round() as usize),
            self.unit
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_and_durations_read_naturally() {
        assert_eq!(with_thousands(0), "0");
        assert_eq!(with_thousands(999), "999");
        assert_eq!(with_thousands(1000), "1,000");
        assert_eq!(with_thousands(1234567), "1,234,567");
        assert_eq!(spoken_duration(8.44), "8.4 s");
        assert_eq!(spoken_duration(161.2), "2 min 41 s");
        assert_eq!(spoken_duration(3900.0), "1 h 05 min");
    }

    #[test]
    fn a_progress_line_tells_share_speed_and_time_left() {
        let mut progress = Progress::new("Aligning", "Aligned", "reads", 200_000);
        progress.on_terminal = false;
        assert_eq!(
            progress.line(50_000, 10.0),
            "Aligning:  25%  50,000 of 200,000 reads, 5,000 reads per second, about 30.0 s left"
        );
        progress.on_terminal = true;
        assert_eq!(
            progress.line(50_000, 10.0),
            "Aligning [######------------------]  25%  50,000 of 200,000 reads, 5,000 reads per second, about 30.0 s left"
        );
        // Nothing done yet: no speed to tell the time left from.
        assert!(!progress.line(0, 0.0).contains("left"));
        progress.add(123);
        assert_eq!(progress.done.load(Ordering::Relaxed), 123);
    }
}
