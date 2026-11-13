//! NumPy 1.26.4's DAS random stream: SeedSequence(pool_size=4, no spawn),
//! PCG64 and unweighted Generator.choice(..., shuffle=True).
//!
//! Adapted to Rust from numpy/random/{bit_generator.pyx,_generator.pyx,
//! src/pcg64/pcg64.{c,h},src/distributions/distributions.c}. Preserve seed-word
//! order, the cached upper uint32, bounded-draw rejection and both choice
//! algorithms. This is not a cryptographic generator or a general NumPy API.
//! See NUMPY_LICENSE and NUMPY_RANDOM_LICENSE for the upstream notices.
use std::collections::BTreeSet;

fn hashmix(mut value: u32, constant: &mut u32) -> u32 {
    value ^= *constant;
    *constant = constant.wrapping_mul(0x931e_8875);
    value = value.wrapping_mul(*constant);
    value ^ (value >> 16)
}

fn mix(left: u32, right: u32) -> u32 {
    let value = left
        .wrapping_mul(0xca01_f9dd)
        .wrapping_sub(right.wrapping_mul(0x4973_f715));
    value ^ (value >> 16)
}

pub(crate) fn seed_sequence(entropy: &[u32], words: usize) -> Vec<u32> {
    let mut constant = 0x43b0_d7e5;
    let mut pool = [0; 4];
    for (at, value) in pool.iter_mut().enumerate() {
        *value = hashmix(entropy.get(at).copied().unwrap_or(0), &mut constant);
    }
    for source in 0..4 {
        for target in 0..4 {
            if source != target {
                pool[target] = mix(pool[target], hashmix(pool[source], &mut constant));
            }
        }
    }
    for &word in entropy.iter().skip(4) {
        for value in &mut pool {
            *value = mix(*value, hashmix(word, &mut constant));
        }
    }
    constant = 0x8b51_f9dd;
    (0..words)
        .map(|at| {
            let mut value = pool[at % 4] ^ constant;
            constant = constant.wrapping_mul(0x58f3_8ded);
            value = value.wrapping_mul(constant);
            value ^ (value >> 16)
        })
        .collect()
}

#[derive(Clone, Debug)]
pub(crate) struct Pcg64 {
    state: u128,
    increment: u128,
    has_uint32: bool,
    buffered_uint32: u32,
}

impl Pcg64 {
    pub(crate) fn new(seed: u64) -> Self {
        let entropy = if seed >> 32 == 0 {
            vec![seed as u32]
        } else {
            vec![seed as u32, (seed >> 32) as u32]
        };
        Self::from_entropy(&entropy)
    }

    pub(crate) fn from_entropy(entropy: &[u32]) -> Self {
        let words = seed_sequence(entropy, 8);
        let pair = |at| words[at] as u64 | ((words[at + 1] as u64) << 32);
        let initial = ((pair(0) as u128) << 64) | pair(2) as u128;
        let sequence = ((pair(4) as u128) << 64) | pair(6) as u128;
        let mut value = Self {
            state: 0,
            increment: sequence.wrapping_shl(1) | 1,
            has_uint32: false,
            buffered_uint32: 0,
        };
        value.step();
        value.state = value.state.wrapping_add(initial);
        value.step();
        value
    }

    fn step(&mut self) {
        const MULTIPLIER: u128 =
            ((2_549_297_995_355_413_924_u128) << 64) | 4_865_540_595_714_422_341;
        self.state = self
            .state
            .wrapping_mul(MULTIPLIER)
            .wrapping_add(self.increment);
    }

    pub(crate) fn state(&self) -> (u128, u128, bool, u32) {
        (
            self.state,
            self.increment,
            self.has_uint32,
            self.buffered_uint32,
        )
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.step();
        ((self.state >> 64) as u64 ^ self.state as u64).rotate_right((self.state >> 122) as u32)
    }

    pub(crate) fn next_u32(&mut self) -> u32 {
        if self.has_uint32 {
            self.has_uint32 = false;
            return self.buffered_uint32;
        }
        let value = self.next_u64();
        self.has_uint32 = true;
        self.buffered_uint32 = (value >> 32) as u32;
        value as u32
    }

    /// An inclusive upper bound, as used by NumPy's choice and int64 fill.
    pub(crate) fn bounded(&mut self, maximum: u64) -> u64 {
        match maximum {
            0 => 0,
            0xffff_ffff => self.next_u32() as u64,
            u64::MAX => self.next_u64(),
            1..=0xffff_fffe => {
                let span = maximum as u32 + 1;
                let threshold = span.wrapping_neg() % span;
                loop {
                    let product = self.next_u32() as u64 * span as u64;
                    if product as u32 >= threshold {
                        break product >> 32;
                    }
                }
            }
            _ => {
                let span = maximum + 1;
                let threshold = span.wrapping_neg() % span;
                loop {
                    let product = self.next_u64() as u128 * span as u128;
                    if product as u64 >= threshold {
                        break (product >> 64) as u64;
                    }
                }
            }
        }
    }

    fn shuffle_tail(&mut self, values: &mut [usize], first: usize) {
        for index in (first..values.len()).rev() {
            let other = self.bounded(index as u64) as usize;
            values.swap(index, other);
        }
    }

    pub(crate) fn choice(
        &mut self,
        population: usize,
        size: usize,
        replace: bool,
    ) -> Result<Vec<usize>, &'static str> {
        if population == 0 && size != 0 {
            return Err("a must be a positive integer unless no samples are taken");
        }
        if !replace && size > population {
            return Err("Cannot take a larger sample than population when replace is False");
        }
        if size == 0 {
            return Ok(vec![]);
        }
        if replace {
            return Ok((0..size)
                .map(|_| self.bounded((population - 1) as u64) as usize)
                .collect());
        }
        if population > 10_000 && size > population / 50 {
            let mut indices: Vec<_> = (0..population).collect();
            self.shuffle_tail(&mut indices, (population - size).max(1));
            return Ok(indices.split_off(population - size));
        }
        // Floyd only queries membership; the set's own storage/iteration order
        // does not affect any draw, collision decision or emitted index.
        let mut seen = BTreeSet::new();
        let mut indices = Vec::with_capacity(size);
        for end in population - size..population {
            let draw = self.bounded(end as u64) as usize;
            if seen.insert(draw) {
                indices.push(draw);
            } else {
                seen.insert(end);
                indices.push(end);
            }
        }
        self.shuffle_tail(&mut indices, 1);
        Ok(indices)
    }
}
