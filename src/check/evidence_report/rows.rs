//! The rows `discover blast` reports today, seen through the evidence.

use std::fmt::Write as _;

use super::{number, Class, Plate};
use crate::evidence::Observation;
use crate::stats::median;

/// One line per class for the sequence-in-donor rows passing the current
/// thresholds: this is what a user of the count thresholds gets.
pub(super) fn write_passing_rows(out: &mut String, plate: &Plate) {
    let _ = writeln!(
        out,
        "\n== Sequence-in-donor rows passing {} (what is reported today)",
        plate.thresholds.label()
    );
    let _ = writeln!(
        out,
        "  class            n   count  dosage  evidence: none    <1   1-2   >=2  not assessed  carriers  does not segregate  chimera"
    );
    for class in Class::ALL {
        let items: Vec<&Observation> = plate
            .of_class(class)
            .filter(|o| plate.thresholds.passes(o))
            .collect();
        if items.is_empty() {
            continue;
        }
        let med = |f: &dyn Fn(&Observation) -> Option<f64>| -> f64 {
            median(&items.iter().filter_map(|o| f(o)).collect::<Vec<f64>>())
        };
        let count = |f: &dyn Fn(&Observation) -> bool| items.iter().filter(|o| f(o)).count();
        let within = |low: f64, high: f64| {
            count(&|o| o.evidence.score.is_some_and(|s| (low..high).contains(&s)))
        };
        let _ = writeln!(
            out,
            "  {:<12} {:5}  {:6.0}  {}           {:5} {:5} {:5} {:5}  {:12}    {:6.0}  {:18}  {:7}",
            class.label(),
            items.len(),
            med(&|o| Some(o.count as f64)),
            number(med(&|o| o.evidence.dosage)),
            count(&|o| o.evidence.score == Some(0.0)),
            within(f64::MIN_POSITIVE, 1.0),
            within(1.0, 2.0),
            within(2.0, f64::INFINITY),
            count(&|o| o.evidence.score.is_none()),
            med(&|o| Some(o.evidence.carriers as f64)),
            count(&|o| o.evidence.ubiquitous),
            count(&|o| o.evidence.chimera),
        );
    }
}
