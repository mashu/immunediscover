//! Gene-specific parameter presets for `discover blast`.

use std::fmt;
use std::str::FromStr;

use crate::info;

/// Defaults of the options a preset can change.
pub mod defaults {
    pub const FORWARD: usize = 20;
    pub const REVERSE: usize = 20;
    pub const MIN_FULL_RATIO: f64 = 0.1;
    pub const LENGTH: usize = 290;
    pub const MAX_DISTANCE: i64 = 20;
    pub const MIN_FULL_COUNT: i64 = 5;
    pub const EDGE: i64 = 0;
    pub const SUBJECT_COVERAGE: f64 = 0.1;
    pub const MIN_QUALITY: f64 = 0.75;
    pub const BLAST_ARGS: &str =
        "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50";
    pub const MIN_CORE_COVERAGE: f64 = 0.6;
}

/// Gene segment a preset is tuned for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    V,
    D,
    J,
}

impl FromStr for Preset {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V" => Ok(Preset::V),
            "D" => Ok(Preset::D),
            "J" => Ok(Preset::J),
            other => Err(format!("'{other}' is not one of V, D, J")),
        }
    }
}

impl fmt::Display for Preset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Preset::V => "V",
            Preset::D => "D",
            Preset::J => "J",
        })
    }
}

/// Values for the options a preset can change; `None` leaves an option
/// alone. Used both for what a preset sets (`PresetValues`) and for what
/// was given explicitly on the command line (`PresetValues<String>`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PresetValues<S = &'static str> {
    pub forward: Option<usize>,
    pub reverse: Option<usize>,
    pub min_full_ratio: Option<f64>,
    pub length: Option<usize>,
    pub max_distance: Option<i64>,
    pub min_full_count: Option<i64>,
    pub edge: Option<i64>,
    pub subject_coverage: Option<f64>,
    pub min_quality: Option<f64>,
    pub min_core_coverage: Option<f64>,
    pub blast_args: Option<S>,
}

impl Preset {
    pub fn values(self) -> PresetValues {
        match self {
            // Tuned for genomic IGHV novel recovery: relaxed trimming and
            // core coverage, lower allelic ratio.
            Preset::V => PresetValues {
                forward: Some(20),
                reverse: Some(20),
                min_full_ratio: Some(0.035),
                length: Some(283),
                max_distance: Some(14),
                min_full_count: Some(5),
                min_quality: Some(0.62),
                min_core_coverage: Some(0.50),
                blast_args: Some("-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"),
                ..PresetValues::default()
            },
            // 40 nt of context on each side splits the reads of a D allele
            // over many extension variants, and the count threshold applies
            // to the reported variant as well as to the allele: hence 5
            // (upstream: 10).
            Preset::D => PresetValues {
                forward: Some(40),
                reverse: Some(40),
                min_full_ratio: Some(0.2),
                length: Some(5),
                max_distance: Some(20),
                min_full_count: Some(5),
                edge: Some(10),
                subject_coverage: Some(0.25),
                min_quality: Some(0.5),
                blast_args: Some(
                    "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5",
                ),
                ..PresetValues::default()
            },
            Preset::J => PresetValues {
                forward: Some(12),
                reverse: Some(12),
                min_full_ratio: Some(0.1),
                length: Some(10),
                max_distance: Some(10),
                min_full_count: Some(10),
                blast_args: Some("-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"),
                ..PresetValues::default()
            },
        }
    }
}

/// Final values of the options a preset can change.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolved {
    pub forward: usize,
    pub reverse: usize,
    pub min_full_ratio: f64,
    pub length: usize,
    pub max_distance: i64,
    pub min_full_count: i64,
    pub edge: i64,
    pub subject_coverage: f64,
    pub min_quality: f64,
    pub min_core_coverage: f64,
    pub blast_args: String,
}

/// Precedence for one option: explicit value, then preset, then default.
fn choose<T: fmt::Display>(
    (preset, log): (Option<Preset>, bool),
    name: &str,
    explicit: Option<T>,
    preset_value: Option<T>,
    default: T,
) -> T {
    match (explicit, preset, preset_value) {
        (Some(value), Some(preset), Some(preset_value)) => {
            if log {
                info!("Keeping explicit {name}: {value} ({preset} preset: {preset_value})");
            }
            value
        }
        (Some(value), _, _) => value,
        (None, Some(preset), Some(preset_value)) => {
            if log {
                info!("Applying {preset} preset {name}: {default} → {preset_value}");
            }
            preset_value
        }
        (None, _, _) => default,
    }
}

/// Combine the explicit command-line values with a preset and the defaults.
///
/// An option given on the command line always wins, whatever its value;
/// a preset only fills in options that were not given.
pub fn resolve(preset: Option<Preset>, explicit: PresetValues<String>) -> Resolved {
    resolve_with(preset, explicit, true)
}

fn resolve_with(preset: Option<Preset>, explicit: PresetValues<String>, log: bool) -> Resolved {
    let v = preset.map(Preset::values).unwrap_or_default();
    let preset = (preset, log);
    Resolved {
        forward: choose(
            preset,
            "forward",
            explicit.forward,
            v.forward,
            defaults::FORWARD,
        ),
        reverse: choose(
            preset,
            "reverse",
            explicit.reverse,
            v.reverse,
            defaults::REVERSE,
        ),
        min_full_ratio: choose(
            preset,
            "minfullratio",
            explicit.min_full_ratio,
            v.min_full_ratio,
            defaults::MIN_FULL_RATIO,
        ),
        length: choose(
            preset,
            "length",
            explicit.length,
            v.length,
            defaults::LENGTH,
        ),
        max_distance: choose(
            preset,
            "maxdist",
            explicit.max_distance,
            v.max_distance,
            defaults::MAX_DISTANCE,
        ),
        min_full_count: choose(
            preset,
            "minfullcount",
            explicit.min_full_count,
            v.min_full_count,
            defaults::MIN_FULL_COUNT,
        ),
        edge: choose(preset, "edge", explicit.edge, v.edge, defaults::EDGE),
        subject_coverage: choose(
            preset,
            "subjectcov",
            explicit.subject_coverage,
            v.subject_coverage,
            defaults::SUBJECT_COVERAGE,
        ),
        min_quality: choose(
            preset,
            "minquality",
            explicit.min_quality,
            v.min_quality,
            defaults::MIN_QUALITY,
        ),
        min_core_coverage: choose(
            preset,
            "min-corecov",
            explicit.min_core_coverage,
            v.min_core_coverage,
            defaults::MIN_CORE_COVERAGE,
        ),
        blast_args: choose(
            preset,
            "args",
            explicit.blast_args,
            v.blast_args.map(str::to_string),
            defaults::BLAST_ARGS.to_string(),
        ),
    }
}

/// The values every preset resolves to, as text (`--show-presets`).
pub fn describe_presets() -> String {
    let mut text = String::new();
    let resolved: Vec<Resolved> = [None, Some(Preset::V), Some(Preset::D), Some(Preset::J)]
        .into_iter()
        .map(|preset| resolve_with(preset, PresetValues::default(), false))
        .collect();
    type Field = fn(&Resolved) -> String;
    let rows: [(&str, Field); 10] = [
        ("--forward", |r| r.forward.to_string()),
        ("--reverse", |r| r.reverse.to_string()),
        ("--minfullratio", |r| r.min_full_ratio.to_string()),
        ("--length", |r| r.length.to_string()),
        ("--maxdist", |r| r.max_distance.to_string()),
        ("--minfullcount", |r| r.min_full_count.to_string()),
        ("--edge", |r| r.edge.to_string()),
        ("--subjectcov", |r| r.subject_coverage.to_string()),
        ("--minquality", |r| r.min_quality.to_string()),
        ("--min-corecov", |r| r.min_core_coverage.to_string()),
    ];
    text.push_str(&format!(
        "{:<15}{:>9}{:>9}{:>9}{:>9}\n",
        "option", "default", "-g V", "-g D", "-g J"
    ));
    for (name, value) in &rows {
        text.push_str(&format!("{name:<15}"));
        for column in &resolved {
            text.push_str(&format!("{:>9}", value(column)));
        }
        text.push('\n');
    }
    for (label, column) in ["default", "-g V", "-g D", "-g J"].iter().zip(&resolved) {
        text.push_str(&format!("--args ({label}): {}\n", column.blast_args));
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_fills_in_what_was_not_given() {
        let explicit = PresetValues {
            max_distance: Some(7),
            ..PresetValues::default()
        };
        let r = resolve(Some(Preset::D), explicit);
        assert_eq!(r.max_distance, 7);
        assert_eq!((r.forward, r.reverse), (40, 40));
        assert_eq!(r.edge, 10);
        assert_eq!(r.length, 5);
        assert!(r.blast_args.contains("-word_size 7"));
        assert_eq!(r.min_core_coverage, defaults::MIN_CORE_COVERAGE);
    }

    #[test]
    fn explicit_values_win_even_when_they_equal_the_default() {
        let explicit = PresetValues {
            forward: Some(defaults::FORWARD),
            min_core_coverage: Some(0.9),
            ..PresetValues::default()
        };
        let r = resolve(Some(Preset::D), explicit);
        assert_eq!(r.forward, defaults::FORWARD);
        assert_eq!(r.reverse, 40);
        let v = resolve(
            Some(Preset::V),
            PresetValues {
                min_core_coverage: Some(0.9),
                ..PresetValues::default()
            },
        );
        assert_eq!(v.min_core_coverage, 0.9);
        assert_eq!(
            resolve(Some(Preset::V), PresetValues::default()).min_core_coverage,
            0.5
        );
    }

    #[test]
    fn without_a_preset_defaults_apply() {
        let r = resolve(None, PresetValues::default());
        assert_eq!(r.length, defaults::LENGTH);
        assert_eq!(r.blast_args, defaults::BLAST_ARGS);
        assert!(describe_presets().contains("--min-corecov"));
    }
}
