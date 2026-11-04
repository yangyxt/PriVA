//! CPython 3.12's 64-bit little-endian string/text-pair set iteration order.
//!
//! The gene-statistics producer constructs dictionaries by iterating newly built
//! sets of transcript/domain strings. Preserve that order for an explicit
//! `PYTHONHASHSEED`; a different process's unspecified random secret is unknown.
//! This small table supports construction from an iterable and the string-set
//! copy/update paths, without deletion or general set algebra. Unicode hashes
//! use Python's one-, two- or four-byte character kind.
//!
//! Algorithm references (CPython v3.12.3): Objects/setobject.c, Objects/tupleobject.c,
//! Python/pyhash.c and Python/bootstrap_hash.c in https://github.com/python/cpython.

#[derive(Clone, Copy)]
pub(crate) struct StringHash {
    keys: [u64; 2],
}

impl StringHash {
    pub(crate) fn seeded(mut seed: u32) -> Self {
        let mut bytes = [0u8; 16];
        if seed != 0 {
            for byte in &mut bytes {
                seed = seed.wrapping_mul(214_013).wrapping_add(2_531_011);
                *byte = (seed >> 16) as u8;
            }
        }
        Self {
            keys: [
                u64::from_le_bytes(bytes[..8].try_into().unwrap()),
                u64::from_le_bytes(bytes[8..].try_into().unwrap()),
            ],
        }
    }

    pub(crate) fn hash(&self, text: &str) -> i64 {
        if text.is_empty() {
            return 0;
        }
        let largest = text.chars().map(u32::from).max().unwrap();
        let width = if largest <= 255 {
            1
        } else if largest <= 65535 {
            2
        } else {
            4
        };
        let mut data = Vec::with_capacity(text.chars().count() * width);
        for character in text.chars() {
            data.extend_from_slice(&u32::from(character).to_le_bytes()[..width]);
        }
        let value = siphash13(self.keys, &data) as i64;
        if value == -1 { -2 } else { value }
    }

    /// `tuplehash` for exactly two Unicode strings, including empty strings.
    /// Hash randomization remains an explicit caller/test profile parameter.
    pub(crate) fn pair_hash(&self, first: &str, second: &str) -> i64 {
        const PRIME_1: u64 = 11_400_714_785_074_694_791;
        const PRIME_2: u64 = 14_029_467_366_897_019_727;
        const PRIME_5: u64 = 2_870_177_450_012_600_261;
        let mut value = PRIME_5;
        for text in [first, second] {
            value = value.wrapping_add((self.hash(text) as u64).wrapping_mul(PRIME_2));
            value = value.rotate_left(31).wrapping_mul(PRIME_1);
        }
        value = value.wrapping_add(2 ^ (PRIME_5 ^ 3_527_539));
        if value == u64::MAX {
            1_546_275_796
        } else {
            value as i64
        }
    }
}

fn round(state: &mut [u64; 4]) {
    let [mut a, mut b, mut c, mut d] = *state;
    a = a.wrapping_add(b);
    c = c.wrapping_add(d);
    b = b.rotate_left(13) ^ a;
    d = d.rotate_left(16) ^ c;
    a = a.rotate_left(32);
    c = c.wrapping_add(b);
    a = a.wrapping_add(d);
    b = b.rotate_left(17) ^ c;
    d = d.rotate_left(21) ^ a;
    c = c.rotate_left(32);
    *state = [a, b, c, d];
}

fn siphash13(keys: [u64; 2], bytes: &[u8]) -> u64 {
    let mut state = [
        keys[0] ^ 0x736f_6d65_7073_6575,
        keys[1] ^ 0x646f_7261_6e64_6f6d,
        keys[0] ^ 0x6c79_6765_6e65_7261,
        keys[1] ^ 0x7465_6462_7974_6573,
    ];
    let mut chunks = bytes.chunks_exact(8);
    for chunk in &mut chunks {
        let word = u64::from_le_bytes(chunk.try_into().unwrap());
        state[3] ^= word;
        round(&mut state);
        state[0] ^= word;
    }
    let mut last = (bytes.len() as u64) << 56;
    for (index, byte) in chunks.remainder().iter().enumerate() {
        last |= u64::from(*byte) << (8 * index);
    }
    state[3] ^= last;
    round(&mut state);
    state[0] ^= last;
    state[2] ^= 255;
    for _ in 0..3 {
        round(&mut state);
    }
    state.into_iter().fold(0, |value, word| value ^ word)
}

#[derive(Clone)]
struct Entry<T> {
    key: T,
    hash: i64,
}

struct Table<T> {
    slots: Vec<Option<Entry<T>>>,
    used: usize,
}

impl<T: Clone + Eq> Table<T> {
    fn new() -> Self {
        Self {
            slots: vec![None; 8],
            used: 0,
        }
    }

    fn locate(&self, key: &T, hash: i64) -> usize {
        let mask = self.slots.len() - 1;
        let mut perturb = hash as u64;
        let mut position = perturb as usize & mask;
        loop {
            let end = if position + 9 <= mask {
                position + 9
            } else {
                position
            };
            for at in position..=end {
                match &self.slots[at] {
                    None => return at,
                    Some(entry) if entry.hash == hash && &entry.key == key => return at,
                    _ => {}
                }
            }
            perturb >>= 5;
            position = position
                .wrapping_mul(5)
                .wrapping_add(1)
                .wrapping_add(perturb as usize)
                & mask;
        }
    }

    fn insert(&mut self, key: T, hash: i64) {
        let position = self.locate(&key, hash);
        if self.slots[position].is_some() {
            return;
        }
        self.slots[position] = Some(Entry { key, hash });
        self.used += 1;
        if self.used * 5 >= (self.slots.len() - 1) * 3 {
            let minimum = self.used * if self.used > 50_000 { 2 } else { 4 };
            let mut size = 8;
            while size <= minimum {
                size *= 2;
            }
            let old = std::mem::replace(&mut self.slots, vec![None; size]);
            // No deletion occurs in set(iterable); table-order reinsertion
            // and the first free probe slot reproduce set_insert_clean.
            for entry in old.into_iter().flatten() {
                let position = self.locate(&entry.key, entry.hash);
                self.slots[position] = Some(entry);
            }
        }
    }

    fn into_order(self) -> Vec<T> {
        self.slots
            .into_iter()
            .flatten()
            .map(|entry| entry.key)
            .collect()
    }

    fn merge(&mut self, other: &Self) {
        if other.used == 0 {
            return;
        }
        // CPython set_merge reserves for both sets before visiting the source
        // table. Repeated ordinary insertion alone has different resize/order
        // behavior, even when every final member is the same.
        if (self.used + other.used) * 5 >= (self.slots.len() - 1) * 3 {
            let minimum = (self.used + other.used) * 2;
            let mut size = 8;
            while size <= minimum {
                size *= 2;
            }
            if size != self.slots.len() {
                let old = std::mem::replace(&mut self.slots, vec![None; size]);
                for entry in old.into_iter().flatten() {
                    let position = self.locate(&entry.key, entry.hash);
                    self.slots[position] = Some(entry);
                }
            }
        }
        // There are no deleted entries in this restricted table. The empty
        // same-size fast path copies slot positions, not iteration insertion.
        if self.used == 0 && self.slots.len() == other.slots.len() {
            self.slots.clone_from(&other.slots);
            self.used = other.used;
        } else if self.used == 0 {
            for entry in other.slots.iter().flatten() {
                let position = self.locate(&entry.key, entry.hash);
                self.slots[position] = Some(entry.clone());
            }
            self.used = other.used;
        } else {
            for entry in other.slots.iter().flatten() {
                self.insert(entry.key.clone(), entry.hash);
            }
        }
    }
}

pub(crate) struct StringSet {
    table: Table<String>,
    hash: StringHash,
}

impl StringSet {
    pub(crate) fn new(seed: u32) -> Self {
        Self {
            table: Table::new(),
            hash: StringHash::seeded(seed),
        }
    }

    pub(crate) fn insert(&mut self, text: String) {
        let hash = self.hash.hash(&text);
        self.table.insert(text, hash);
    }

    pub(crate) fn into_order(self) -> Vec<String> {
        self.table.into_order()
    }

    #[allow(dead_code)] // Serialization observers need a nonmutating table-order view.
    pub(crate) fn iter(&self) -> impl Iterator<Item = &str> {
        self.table
            .slots
            .iter()
            .flatten()
            .map(|entry| entry.key.as_str())
    }

    /// `set.update(other_set)` for two explicit, identical Python hash profiles.
    /// A cross-process random hash secret is not inferred from observed values.
    #[allow(dead_code)] // Used by the separate native Alpha exon-history path.
    pub(crate) fn update_set(&mut self, other: &Self) -> Result<(), &'static str> {
        if self.hash.keys != other.hash.keys {
            return Err("set update requires identical explicit hash profiles");
        }
        self.table.merge(&other.table);
        Ok(())
    }

    /// Python `set.copy()` constructs an empty table and takes the set fast path;
    /// a Rust Clone of the slots would preserve too much allocation history.
    #[allow(dead_code)] // Used by the separate native Alpha exon-history path.
    pub(crate) fn copy(&self) -> Self {
        let mut copied = Self {
            table: Table::new(),
            hash: self.hash,
        };
        copied.table.merge(&self.table);
        copied
    }

    /// CPython's optimized constant set literal: BUILD_SET(0), then SET_UPDATE
    /// from the compiler's frozenset. Keep this separate from set(iterable).
    pub(crate) fn literal(seed: u32, values: &[&str]) -> Self {
        let mut frozen = Self::new(seed);
        for value in values {
            frozen.insert((*value).into());
        }
        let mut result = Self::new(seed);
        if frozen.table.used * 5 >= 7 * 3 {
            let minimum = frozen.table.used * 2;
            let mut size = 8;
            while size <= minimum {
                size *= 2;
            }
            result.table.slots = vec![None; size];
        }
        for entry in frozen.table.slots.into_iter().flatten() {
            let index = result.table.locate(&entry.key, entry.hash);
            result.table.slots[index] = Some(entry);
            result.table.used += 1;
        }
        result
    }

    /// Fresh set-difference iteration used by the installer's five-key constant
    /// minus one-key literal. No dummy-slot mutation or general set deletion.
    pub(crate) fn difference_small(self, excluded: &[&str]) -> Vec<String> {
        let mut result = Self {
            table: Table::new(),
            hash: self.hash,
        };
        for entry in self.table.slots.into_iter().flatten() {
            if !excluded.contains(&entry.key.as_str()) {
                result.table.insert(entry.key, entry.hash);
            }
        }
        result.into_order()
    }
}

/// A newly constructed Python set of exact two-text tuples. No deletion or
/// set algebra is modeled; normalization inserts the tuples in source order.
pub(crate) struct PairSet {
    table: Table<(String, String)>,
    hash: StringHash,
}

impl PairSet {
    pub(crate) fn new(seed: u32) -> Self {
        Self {
            table: Table::new(),
            hash: StringHash::seeded(seed),
        }
    }

    pub(crate) fn insert(&mut self, pair: (String, String)) {
        let hash = self.hash.pair_hash(&pair.0, &pair.1);
        self.table.insert(pair, hash);
    }

    pub(crate) fn into_order(self) -> Vec<(String, String)> {
        self.table.into_order()
    }
}
