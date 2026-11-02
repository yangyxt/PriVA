//! CPython 3.12's 64-bit little-endian, SipHash13 string-set iteration order.
//!
//! The gene-statistics producer constructs dictionaries by iterating newly built
//! sets of transcript/domain strings. Preserve that order for an explicit
//! `PYTHONHASHSEED`; a different process's unspecified random secret is unknown.
//! This small table supports construction from an iterable, not deletion or set
//! algebra. Unicode hashes use Python's one-, two- or four-byte character kind.
//!
//! Algorithm references (CPython v3.12.3): Objects/setobject.c,
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
struct Entry {
    text: String,
    hash: i64,
}

pub(crate) struct StringSet {
    slots: Vec<Option<Entry>>,
    used: usize,
    hash: StringHash,
}

impl StringSet {
    pub(crate) fn new(seed: u32) -> Self {
        Self {
            slots: vec![None; 8],
            used: 0,
            hash: StringHash::seeded(seed),
        }
    }

    fn locate(&self, text: &str, hash: i64) -> usize {
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
                    Some(entry) if entry.hash == hash && entry.text == text => return at,
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

    pub(crate) fn insert(&mut self, text: String) {
        let hash = self.hash.hash(&text);
        let position = self.locate(&text, hash);
        if self.slots[position].is_some() {
            return;
        }
        self.slots[position] = Some(Entry { text, hash });
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
                let position = self.locate(&entry.text, entry.hash);
                self.slots[position] = Some(entry);
            }
        }
    }

    pub(crate) fn into_order(self) -> Vec<String> {
        self.slots
            .into_iter()
            .flatten()
            .map(|entry| entry.text)
            .collect()
    }
}
