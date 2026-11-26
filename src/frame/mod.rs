//! A compact typed table (`Frame`) with the relational operations the
//! pipelines need, mirroring the semantics of CSV.jl and DataFrames.jl.

mod column;
mod float_fmt;
mod group;
mod infer;
mod join;
mod sort;
mod table;
mod value;

pub use column::Column;
pub use float_fmt::{format_float, round_digits};
pub use group::Groups;
pub use infer::infer_column;
pub use join::{join, semi_join_rows, JoinKind};
pub use sort::SortKey;
pub use table::{Frame, Pool};
pub use value::{Key, Value};
