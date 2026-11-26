//! Sequence-level domain primitives shared by all commands.

pub mod database;
pub mod demux;
pub mod dna;
pub mod gene;
pub mod naming;
pub mod search;
pub mod slice;

pub use database::{load_database, Allele};
pub use demux::{load_demultiplex, DemuxTable};
pub use gene::{gene_of, GeneType};
pub use naming::unique_name;
