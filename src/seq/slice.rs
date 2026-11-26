//! 1-based, inclusive string ranges with Julia semantics.

/// `s[first:last]` as Julia evaluates it for ASCII text: an empty string
/// when `last < first`, otherwise the inclusive 1-based byte range.
///
/// Bounds are clamped rather than panicking; callers pass ranges that the
/// original code already clamps, so clamping never changes a result.
pub fn range(s: &str, first: i64, last: i64) -> &str {
    if last < first {
        return "";
    }
    let len = s.len() as i64;
    let start = (first.max(1) - 1).min(len) as usize;
    let end = last.min(len).max(0) as usize;
    if end <= start {
        ""
    } else {
        &s[start..end]
    }
}

#[cfg(test)]
mod tests {
    use super::range;

    #[test]
    fn behaves_like_julia_ranges() {
        let s = "ABCDEFG";
        assert_eq!(range(s, 2, 4), "BCD");
        assert_eq!(range(s, 5, 4), "");
        assert_eq!(range(s, 1, 0), "");
        assert_eq!(range(s, 8, 7), "");
        assert_eq!(range(s, 1, 1), "A");
        assert_eq!(range(s, 7, 7), "G");
    }
}
