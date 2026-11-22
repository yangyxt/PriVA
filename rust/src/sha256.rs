//! SHA-256 (FIPS 180-4) for checksums of cache files, computed in blocks so
//! files of any size can be hashed without reading them into memory.
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

const ROUND_CONSTANTS: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

const INITIAL_STATE: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

/// Block size of the file reads, as the Python script copies and hashes.
const READ_BLOCK: usize = 8 * 1024 * 1024;

pub(crate) struct Sha256 {
    state: [u32; 8],
    pending: [u8; 64],
    pending_len: usize,
    total_len: u64,
}

fn compress(state: &mut [u32; 8], block: &[u8]) {
    let mut words = [0u32; 64];
    for (index, chunk) in block.chunks_exact(4).enumerate() {
        words[index] = u32::from_be_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
    }
    for index in 16..64 {
        let w15 = words[index - 15];
        let w2 = words[index - 2];
        let small0 = w15.rotate_right(7) ^ w15.rotate_right(18) ^ (w15 >> 3);
        let small1 = w2.rotate_right(17) ^ w2.rotate_right(19) ^ (w2 >> 10);
        words[index] = words[index - 16]
            .wrapping_add(small0)
            .wrapping_add(words[index - 7])
            .wrapping_add(small1);
    }
    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
    for index in 0..64 {
        let big1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let choose = (e & f) ^ (!e & g);
        let temp1 = h
            .wrapping_add(big1)
            .wrapping_add(choose)
            .wrapping_add(ROUND_CONSTANTS[index])
            .wrapping_add(words[index]);
        let big0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let majority = (a & b) ^ (a & c) ^ (b & c);
        let temp2 = big0.wrapping_add(majority);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(temp1);
        d = c;
        c = b;
        b = a;
        a = temp1.wrapping_add(temp2);
    }
    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
    state[5] = state[5].wrapping_add(f);
    state[6] = state[6].wrapping_add(g);
    state[7] = state[7].wrapping_add(h);
}

impl Sha256 {
    pub(crate) fn new() -> Sha256 {
        Sha256 {
            state: INITIAL_STATE,
            pending: [0u8; 64],
            pending_len: 0,
            total_len: 0,
        }
    }

    pub(crate) fn update(&mut self, mut data: &[u8]) {
        self.total_len = self.total_len.wrapping_add(data.len() as u64);
        if self.pending_len > 0 {
            let room = 64 - self.pending_len;
            let take = room.min(data.len());
            self.pending[self.pending_len..self.pending_len + take].copy_from_slice(&data[..take]);
            self.pending_len += take;
            data = &data[take..];
            if self.pending_len == 64 {
                let block = self.pending;
                compress(&mut self.state, &block);
                self.pending_len = 0;
            }
        }
        while data.len() >= 64 {
            compress(&mut self.state, &data[..64]);
            data = &data[64..];
        }
        if !data.is_empty() {
            self.pending[..data.len()].copy_from_slice(data);
            self.pending_len = data.len();
        }
    }

    /// The digest as 64 lowercase hexadecimal characters.
    pub(crate) fn finish_hex(mut self) -> String {
        let bit_length = self.total_len.wrapping_mul(8);
        let mut padding = vec![0x80u8];
        let used = (self.pending_len + 1) % 64;
        let zeros = if used <= 56 { 56 - used } else { 120 - used };
        padding.extend(std::iter::repeat_n(0u8, zeros));
        padding.extend_from_slice(&bit_length.to_be_bytes());
        // `update` also counts the padding, which no longer matters.
        self.update(&padding);
        let mut hex = String::with_capacity(64);
        for word in self.state.iter() {
            hex.push_str(&format!("{word:08x}"));
        }
        hex
    }
}

pub(crate) fn sha256_hex(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finish_hex()
}

/// The SHA-256 of a whole file, read in blocks.
pub(crate) fn sha256_file(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut block = vec![0u8; READ_BLOCK];
    loop {
        let read = file.read(&mut block)?;
        if read == 0 {
            break;
        }
        hasher.update(&block[..read]);
    }
    Ok(hasher.finish_hex())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_digests() {
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    #[test]
    fn chunked_updates_match_one_update() {
        let data: Vec<u8> = (0..1000u32).map(|n| (n % 251) as u8).collect();
        let mut hasher = Sha256::new();
        for part in data.chunks(37) {
            hasher.update(part);
        }
        assert_eq!(hasher.finish_hex(), sha256_hex(&data));
        let million_a = vec![b'a'; 1_000_000];
        assert_eq!(
            sha256_hex(&million_a),
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
        );
    }
}
