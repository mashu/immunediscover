//! immunediscover: germline immune receptor allele discovery from NGS data.
//!
//! This crate is a Rust port of `immunediscover.jl`. The layering is:
//!
//! * [`frame`], [`io`] — typed tables and file formats;
//! * [`seq`], [`align`], [`stats`], [`filters`] — domain primitives and algorithms;
//! * [`commands`] — one module per CLI command, usable as a library;
//! * [`cli`] — argument parsing and dispatch.

pub mod align;
pub mod check;
pub mod cli;
pub mod commands;
pub mod error;
pub mod filters;
pub mod frame;
pub mod io;
pub mod logging;
pub mod seq;
pub mod stats;

pub use error::{Error, Result};
