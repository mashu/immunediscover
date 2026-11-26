//! Statistical building blocks.

pub mod descriptive;
pub mod hclust;
pub mod hypergeometric;
pub mod multiple_testing;
pub mod robust;
pub mod tails;

pub use descriptive::{mean, median, sample_std};
pub use hclust::{hclust, Dendrogram, Linkage};
pub use multiple_testing::benjamini_hochberg;
