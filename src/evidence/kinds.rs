//! Which artifact could have produced a sequence.

use super::Observation;

/// A sequence in a donor, by its relation to the more abundant sequences
/// of the same gene in that donor. Every kind has its own artifacts, and
/// its own rate of them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Kind {
    /// One, two or three differences from a more abundant sequence.
    Near(usize),
    /// Four or more differences from the nearest more abundant sequence.
    Far,
    /// Starts like one more abundant sequence and ends like another.
    Chimera,
    /// The most abundant sequence of its gene in the donor.
    NoParent,
}

impl Kind {
    pub const ALL: [Kind; 6] = [
        Kind::Near(1),
        Kind::Near(2),
        Kind::Near(3),
        Kind::Far,
        Kind::Chimera,
        Kind::NoParent,
    ];

    pub fn of(o: &Observation) -> Kind {
        match o.evidence.parent {
            _ if o.evidence.chimera => Kind::Chimera,
            None => Kind::NoParent,
            Some((_, d @ 1..=3)) => Kind::Near(d),
            Some(_) => Kind::Far,
        }
    }

    pub fn index(self) -> usize {
        Kind::ALL
            .iter()
            .position(|k| *k == self)
            .expect("listed in ALL")
    }

    pub fn label(self) -> &'static str {
        match self {
            Kind::Near(1) => "1 difference",
            Kind::Near(2) => "2 differences",
            Kind::Near(_) => "3 differences",
            Kind::Far => "4+ differences",
            Kind::Chimera => "chimera",
            Kind::NoParent => "top of gene",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::testing::observation;

    #[test]
    fn kinds() {
        let mut o = observation(0, "A", "x", 10);
        o.evidence.parent = Some((0, 2));
        assert_eq!(Kind::of(&o), Kind::Near(2));
        o.evidence.parent = Some((0, 7));
        assert_eq!(Kind::of(&o), Kind::Far);
        o.evidence.chimera = true;
        assert_eq!(Kind::of(&o), Kind::Chimera);
        o.evidence.chimera = false;
        o.evidence.parent = None;
        assert_eq!(Kind::of(&o), Kind::NoParent);
    }
}
