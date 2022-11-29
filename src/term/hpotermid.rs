use core::fmt::Debug;
use std::fmt::Display;


#[derive(Copy, Clone, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct HpoTermId {
    inner: u32,
}

impl HpoTermId {
    fn new(s: &str) -> HpoTermId {
        HpoTermId {
            inner: s[3..].parse::<u32>().unwrap(),
        }
    }

    pub fn to_usize(&self) -> usize {
        self.inner.try_into().unwrap()
    }
}

impl From<&str> for HpoTermId {
    fn from(s: &str) -> Self {
        HpoTermId::new(s)
    }
}

impl From<usize> for HpoTermId {
    fn from(n: usize) -> Self {
        Self { inner: n.try_into().unwrap() }
    }
}

impl From<[char; 10]> for HpoTermId {
    fn from(s: [char; 10]) -> Self {
        let mut num = String::with_capacity(7);
        for c in &s[3..] {
            num.push(*c);
        }
        HpoTermId { inner: num.parse::<u32>().unwrap() }
    }
}

impl Debug for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HpoTermId({})", self)
    }
}

impl Display for HpoTermId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HP:{:07}", self.inner)
    }
}

impl PartialEq<str> for HpoTermId {
    fn eq(&self, other: &str) -> bool {
        self == &HpoTermId::new(other)
    }
}