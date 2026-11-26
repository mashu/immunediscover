//! A set of keys with one associated value each (sequence -> allele name).

use std::collections::HashMap;

use crate::{info, warn};

/// Insertion-ordered key/value set. The first value seen for a key wins.
#[derive(Debug, Clone, Default)]
pub struct KeyedSet {
    order: Vec<String>,
    values: HashMap<String, String>,
}

impl KeyedSet {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build from `(key, value)` pairs.
    pub fn from_pairs<I, K, V>(pairs: I) -> Self
    where
        I: IntoIterator<Item = (K, V)>,
        K: Into<String>,
        V: Into<String>,
    {
        let mut set = KeyedSet::new();
        for (key, value) in pairs {
            set.insert(key.into(), value.into());
        }
        set
    }

    /// Add a pair; an existing key keeps its value (a differing value is
    /// reported as a warning).
    pub fn insert(&mut self, key: String, value: String) {
        match self.values.get(&key) {
            None => {
                self.order.push(key.clone());
                self.values.insert(key, value);
            }
            Some(existing) if *existing != value => {
                warn!(
                    "Key {key} with value {value} already exists in KeyedSet with value {existing}"
                );
            }
            Some(_) => info!("Duplicate key {key} with value {value} already exists in KeyedSet"),
        }
    }

    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    pub fn contains(&self, key: &str) -> bool {
        self.values.contains_key(key)
    }

    pub fn get(&self, key: &str) -> Option<&str> {
        self.values.get(key).map(String::as_str)
    }

    /// `(key, value)` pairs in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &str)> + '_ {
        self.order
            .iter()
            .map(move |k| (k.as_str(), self.values[k].as_str()))
    }

    /// Keys of either set; values of `self` take precedence.
    pub fn union(&self, other: &KeyedSet) -> KeyedSet {
        let mut out = self.clone();
        for (key, value) in other.iter() {
            match out.get(key) {
                Some(existing) if existing != value => info!(
                    "Key {key} exists in both sets with different names. Using name '{existing}' from the first set instead of '{value}' from the second set."
                ),
                Some(_) => {}
                None => out.insert(key.to_string(), value.to_string()),
            }
        }
        out
    }

    /// Keys present in both sets, with the values of `self`.
    pub fn intersection(&self, other: &KeyedSet) -> KeyedSet {
        let mut out = KeyedSet::new();
        for (key, value) in self.iter() {
            if let Some(theirs) = other.get(key) {
                if theirs != value {
                    info!("Sequence with key {key} is identical but has different names: {value} in set 1, {theirs} in set 2");
                }
                out.insert(key.to_string(), value.to_string());
            }
        }
        out
    }

    /// Keys of `self` that are absent from `other`.
    pub fn difference(&self, other: &KeyedSet) -> KeyedSet {
        KeyedSet::from_pairs(self.iter().filter(|(key, _)| !other.contains(key)))
    }
}

impl PartialEq for KeyedSet {
    fn eq(&self, other: &Self) -> bool {
        self.values == other.values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Mirrors the upstream `keyedsets.jl` test set.
    #[test]
    fn set_operations() {
        let a = KeyedSet::from_pairs([("key1", "value1"), ("key2", "value2")]);
        let b = KeyedSet::from_pairs([("key2", "value2"), ("key3", "value3")]);
        assert!(a.contains("key1") && !a.contains("key3"));
        assert_eq!(a.union(&b).len(), 3);
        let both = a.intersection(&b);
        assert_eq!(both.iter().collect::<Vec<_>>(), [("key2", "value2")]);
        let only_a = a.difference(&b);
        assert_eq!(only_a.iter().collect::<Vec<_>>(), [("key1", "value1")]);
    }

    #[test]
    fn first_value_wins_and_equality_ignores_order() {
        let mut set = KeyedSet::new();
        set.insert("k".into(), "v1".into());
        set.insert("k".into(), "v2".into());
        assert_eq!(set.get("k"), Some("v1"));
        assert_eq!(set.len(), 1);
        let x = KeyedSet::from_pairs([("a", "1"), ("b", "2")]);
        let y = KeyedSet::from_pairs([("b", "2"), ("a", "1")]);
        assert_eq!(x, y);
        assert_ne!(x, KeyedSet::from_pairs([("a", "1"), ("c", "3")]));
    }
}
