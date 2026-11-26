//! Small text reports printed to stdout.

use std::fmt::Display;

/// Print a two-column label/count table.
///
/// This is the plain-text fallback the original uses when UnicodePlots is
/// not loaded, which is the case in its released binaries.
pub fn print_counts_table<L: Display, C: Display>(labels: &[L], counts: &[C]) {
    print!("{}", counts_table(labels, counts));
}

fn counts_table<L: Display, C: Display>(labels: &[L], counts: &[C]) -> String {
    if labels.is_empty() {
        return String::new();
    }
    let labels: Vec<String> = labels.iter().map(ToString::to_string).collect();
    let counts: Vec<String> = counts.iter().map(ToString::to_string).collect();
    let label_width = labels
        .iter()
        .map(|l| l.chars().count())
        .max()
        .unwrap_or(0)
        .max(12);
    let count_width = counts.iter().map(|c| c.chars().count()).max().unwrap_or(0);
    let mut out = String::new();
    for (label, count) in labels.iter().zip(&counts) {
        out.push_str(&format!("  {label:<label_width$}  {count:>count_width$}\n"));
    }
    out.push('\n');
    out
}

#[cfg(test)]
mod tests {
    use super::counts_table;

    #[test]
    fn pads_labels_to_at_least_twelve_columns() {
        let text = counts_table(&["Donor1", "D2"], &[1200, 7]);
        assert_eq!(text, "  Donor1        1200\n  D2               7\n\n");
        assert_eq!(counts_table::<&str, i32>(&[], &[]), "");
    }
}
