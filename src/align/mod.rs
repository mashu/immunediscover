//! Pairwise alignment kernels.

pub mod banded;
pub mod distance;
pub mod edits;
pub mod genome;
pub mod gotoh;
pub mod needleman_wunsch;
pub mod search;
pub mod seeds;
pub mod submat;

pub use banded::{place_near_diagonal, Placement};
pub use distance::{hamming, levenshtein};
pub use needleman_wunsch::{align, AlignedPair, Alignment, GapScores};
pub use submat::{ednafull, EncodedDna};

/// Semi-global alignment as BioAlignments.jl defines it: all of `a` is
/// aligned, leading and trailing parts of `b` may stay unaligned for free.
/// The returned columns still cover the whole of `b`.
pub fn semiglobal(a: &EncodedDna, b: &EncodedDna, gap_open: i32, gap_extend: i32) -> Alignment {
    align(
        a.codes(),
        b.codes(),
        ednafull,
        GapScores::free_ends(gap_open, gap_extend),
        GapScores::uniform(gap_open, gap_extend),
    )
}
