//! The evidence columns of the output: how much of every reported
//! sequence a donor holds, and what the rest of the plate says about it.
//!
//! The evidence is computed on all clusters of the run, before the
//! abundance cuts: the sequences that are cut are what the artifacts of
//! the plate look like.

use std::collections::HashMap;

use crate::error::Result;
use crate::evidence::{self, Kind, Observation};
use crate::frame::{Column, Frame};

/// Names of the columns added, in output order.
pub const COLUMNS: [&str; 8] = [
    "db_distance",
    "copies",
    "evidence",
    "relation",
    "donors",
    "carriers",
    "segregates",
    "usual_parent_ratio",
];

fn text(frame: &Frame, name: &str) -> Result<Vec<String>> {
    Ok(frame
        .column(name)?
        .iter()
        .map(|v| {
            if v.is_missing() {
                String::new()
            } else {
                v.julia_string()
            }
        })
        .collect())
}

/// Three significant digits: the estimates are not more precise.
fn rounded(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let scale = 10f64.powi(2 - x.abs().log10().floor() as i32);
    (x * scale).round() / scale
}

fn relation(o: &Observation) -> &'static str {
    match Kind::of(o) {
        Kind::NoParent => "top of gene",
        kind => kind.label(),
    }
}

/// Add the [`COLUMNS`] to the reported clusters. `all_clusters` are the
/// clusters of the whole run before the abundance cuts.
pub fn add_evidence_columns(reported: &mut Frame, all_clusters: &Frame) -> Result<()> {
    let (mut observations, donors) =
        evidence::observations::from_clusters_with_donors(all_clusters)?;
    evidence::assess(&mut observations, donors.len());

    let mut seen_in: HashMap<&str, usize> = HashMap::new();
    let mut by_key: HashMap<(&str, &str, &str), &Observation> = HashMap::new();
    for o in &observations {
        *seen_in.entry(o.sequence.as_str()).or_insert(0) += 1;
        let (well, case) = &donors[o.donor];
        by_key.insert((well.as_str(), case.as_str(), o.sequence.as_str()), o);
    }

    let (wells, cases, sequences) = (
        text(reported, "well")?,
        text(reported, "case")?,
        text(reported, "aln_qseq")?,
    );
    let rows: Vec<Option<&Observation>> = (0..reported.nrows())
        .map(|row| {
            by_key
                .get(&(
                    wells[row].as_str(),
                    cases[row].as_str(),
                    sequences[row].as_str(),
                ))
                .copied()
        })
        .collect();
    let floats = |value: &dyn Fn(&Observation) -> Option<f64>| {
        Column::Float(
            rows.iter()
                .map(|o| o.and_then(value).map(rounded))
                .collect(),
        )
    };
    let ints = |value: &dyn Fn(&Observation) -> i64| {
        Column::Int(rows.iter().map(|o| o.map(value)).collect())
    };

    reported.set_column("db_distance", ints(&|o| o.db_distance as i64))?;
    reported.set_column("copies", floats(&|o| o.evidence.dosage))?;
    reported.set_column("evidence", floats(&|o| o.evidence.score))?;
    reported.set_column(
        "relation",
        Column::Str(
            rows.iter()
                .map(|o| o.map(|o| relation(o).to_string()))
                .collect(),
        ),
    )?;
    reported.set_column("donors", ints(&|o| seen_in[o.sequence.as_str()] as i64))?;
    reported.set_column("carriers", ints(&|o| o.evidence.carriers as i64))?;
    reported.set_column(
        "segregates",
        Column::Bool(
            rows.iter()
                .map(|o| o.map(|o| !o.evidence.ubiquitous))
                .collect(),
        ),
    )?;
    reported.set_column(
        "usual_parent_ratio",
        floats(&|o| o.evidence.recurrent_ratio.filter(|r| *r > 0.0)),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Six donors with 200 reads of the allele of gene A; donor 1 also has
    /// 100 reads of a variant, donor 2 three reads of it.
    fn clusters() -> Frame {
        let mut wells = Vec::new();
        let mut sequences = Vec::new();
        let mut counts = Vec::new();
        let mut distances = Vec::new();
        for well in 1..=6 {
            wells.push(well);
            sequences.push("ACGTACGT");
            counts.push(200);
            distances.push(0);
            if well <= 2 {
                wells.push(well);
                sequences.push("ACGTACGA");
                counts.push(if well == 1 { 100 } else { 3 });
                distances.push(1);
            }
        }
        let n = wells.len();
        Frame::from_columns([
            ("well", Column::from_ints(wells)),
            ("case", Column::from_strings(vec!["D"; n])),
            ("sseqid", Column::from_strings(vec!["A*01"; n])),
            ("full_count", Column::from_ints(counts.clone())),
            (
                "full_ratio",
                Column::from_floats(counts.iter().map(|c| *c as f64 / 200.0)),
            ),
            (
                "ratio",
                Column::from_floats(counts.iter().map(|c| *c as f64 / 200.0)),
            ),
            ("aln_qseq", Column::from_strings(sequences)),
            ("aln_mismatch", Column::from_ints(distances)),
        ])
        .unwrap()
    }

    #[test]
    fn reported_rows_get_the_evidence_of_their_sequence_in_their_donor() {
        let all = clusters();
        // Reported: the allele and the variant in donor 1.
        let mut reported = all.filter(&[true, true, false, false, false, false, false, false]);
        add_evidence_columns(&mut reported, &all).unwrap();
        for name in COLUMNS {
            assert!(reported.column(name).is_ok(), "{name}");
        }
        let floats = |name: &str| reported.column(name).unwrap().to_f64().unwrap();
        // 200 reads are two copies; the variant holds one.
        assert_eq!(floats("copies"), [Some(2.0), Some(1.0)]);
        assert_eq!(reported.ints("db_distance").unwrap(), [Some(0), Some(1)]);
        assert_eq!(reported.ints("donors").unwrap(), [Some(6), Some(2)]);
        assert_eq!(reported.ints("carriers").unwrap(), [Some(6), Some(1)]);
        let relation: Vec<String> = text(&reported, "relation").unwrap();
        assert_eq!(relation, ["top of gene", "1 difference"]);
        assert_eq!(
            *reported.column("segregates").unwrap(),
            Column::Bool(vec![Some(false), Some(true)])
        );
        // No artifacts were seen on this tiny plate, so their tail is the
        // default one and the evidence modest.
        assert!(floats("evidence")[1].unwrap() > 0.5);
    }

    #[test]
    fn three_significant_digits() {
        assert_eq!(rounded(1.23456), 1.23);
        assert_eq!(rounded(0.0123456), 0.0123);
        assert_eq!(rounded(1234.56), 1230.0);
        assert_eq!(rounded(0.0), 0.0);
    }
}
