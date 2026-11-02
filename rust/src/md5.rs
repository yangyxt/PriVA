//! Streaming MD5 for matching existing publisher checksums in installer builders.
//! Implements the byte-oriented algorithm in RFC 1321, section 3:
//! <https://www.rfc-editor.org/rfc/rfc1321#section-3>.
//! This does not replace the SHA-256 identities used for native-cache provenance.
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

const INITIAL: [u32; 4] = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476];
const ROTATE: [[u32; 4]; 4] = [
    [7, 12, 17, 22],
    [5, 9, 14, 20],
    [4, 11, 16, 23],
    [6, 10, 15, 21],
];
const CONSTANTS: [u32; 64] = [
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
    0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
    0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
    0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
    0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
    0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
    0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391,
];

pub(crate) struct Md5 {
    words: [u32; 4],
    pending: [u8; 64],
    pending_len: usize,
    // Byte count modulo 2^64; finish uses its low 64 bits after multiplying by 8.
    length: u64,
}

fn compress(state: &mut [u32; 4], block: &[u8]) {
    let mut message = [0u32; 16];
    for (slot, bytes) in message.iter_mut().zip(block.chunks_exact(4)) {
        *slot = u32::from_le_bytes(bytes.try_into().unwrap());
    }
    let [mut a, mut b, mut c, mut d] = *state;
    for (index, constant) in CONSTANTS.iter().enumerate() {
        let (value, word) = match index / 16 {
            0 => ((b & c) | (!b & d), index),
            1 => ((b & d) | (c & !d), (5 * index + 1) % 16),
            2 => (b ^ c ^ d, (3 * index + 5) % 16),
            _ => (c ^ (b | !d), (7 * index) % 16),
        };
        let next = b.wrapping_add(
            a.wrapping_add(value)
                .wrapping_add(*constant)
                .wrapping_add(message[word])
                .rotate_left(ROTATE[index / 16][index % 4]),
        );
        (a, b, c, d) = (d, next, b, c);
    }
    for (word, update) in state.iter_mut().zip([a, b, c, d]) {
        *word = word.wrapping_add(update);
    }
}

impl Md5 {
    pub(crate) fn new() -> Self {
        Self {
            words: INITIAL,
            pending: [0; 64],
            pending_len: 0,
            length: 0,
        }
    }

    pub(crate) fn update(&mut self, mut bytes: &[u8]) {
        self.length = self.length.wrapping_add(bytes.len() as u64);
        if self.pending_len != 0 {
            let count = (64 - self.pending_len).min(bytes.len());
            self.pending[self.pending_len..self.pending_len + count]
                .copy_from_slice(&bytes[..count]);
            self.pending_len += count;
            bytes = &bytes[count..];
            if self.pending_len == 64 {
                compress(&mut self.words, &self.pending);
                self.pending_len = 0;
            }
        }
        while bytes.len() >= 64 {
            compress(&mut self.words, &bytes[..64]);
            bytes = &bytes[64..];
        }
        if !bytes.is_empty() {
            self.pending[..bytes.len()].copy_from_slice(bytes);
            self.pending_len = bytes.len();
        }
    }

    pub(crate) fn finish_hex(mut self) -> String {
        let bits = self.length.wrapping_mul(8);
        // At most two blocks: delimiter, zeros, then the original bit length.
        let mut tail = [0u8; 128];
        tail[0] = 0x80;
        let count = if self.pending_len < 56 { 64 } else { 128 } - self.pending_len;
        tail[count - 8..count].copy_from_slice(&bits.to_le_bytes());
        self.update(&tail[..count]);
        let mut result = String::with_capacity(32);
        for byte in self.words.iter().flat_map(|word| word.to_le_bytes()) {
            use std::fmt::Write;
            write!(&mut result, "{byte:02x}").unwrap();
        }
        result
    }
}

pub(crate) fn md5_hex(bytes: &[u8]) -> String {
    let mut digest = Md5::new();
    digest.update(bytes);
    digest.finish_hex()
}

/// Hash a readable stream without retaining its contents. An interrupted read
/// is retried; any other I/O error is returned instead of a partial digest.
pub(crate) fn md5_reader(mut input: impl Read) -> io::Result<String> {
    let mut digest = Md5::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        match input.read(&mut buffer) {
            Ok(0) => return Ok(digest.finish_hex()),
            Ok(count) => digest.update(&buffer[..count]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        }
    }
}

pub(crate) fn md5_file(path: &Path) -> io::Result<String> {
    md5_reader(File::open(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rfc_vectors() {
        for (input, expected) in [
            ("", "d41d8cd98f00b204e9800998ecf8427e"),
            ("a", "0cc175b9c0f1b6a831c399e269772661"),
            ("abc", "900150983cd24fb0d6963f7d28e17f72"),
            ("message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
            (
                "abcdefghijklmnopqrstuvwxyz",
                "c3fcd3d76192e4007dfb496cca67e13b",
            ),
            (
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789",
                "d174ab98d277d9f5a5611c2c9f419d9f",
            ),
            (
                "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
                "57edf4a22be3c955ac49da2e2107b67a",
            ),
        ] {
            assert_eq!(md5_hex(input.as_bytes()), expected);
        }
    }

    #[test]
    fn partitioning_and_padding_boundaries() {
        let input: Vec<_> = (0..1024).map(|x| (x % 251) as u8).collect();
        for length in 0..=256 {
            let expected = md5_hex(&input[..length]);
            for chunk_size in [1, 2, 7, 55, 56, 57, 63, 64, 65, 127, 128] {
                let mut state = Md5::new();
                for part in input[..length].chunks(chunk_size) {
                    state.update(part);
                    state.update(&[]);
                }
                assert_eq!(
                    state.finish_hex(),
                    expected,
                    "length={length}, chunk={chunk_size}"
                );
            }
        }
        assert_eq!(
            md5_hex(&vec![b'a'; 1_000_000]),
            "7707d6ae4e027c70eea2a935c2296f21"
        );
    }

    #[test]
    fn read_errors_do_not_publish_partial_digest() {
        struct InterruptedOnce(bool);
        impl Read for InterruptedOnce {
            fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
                if !self.0 {
                    self.0 = true;
                    return Err(io::ErrorKind::Interrupted.into());
                }
                b"".as_slice().read(output)
            }
        }
        assert_eq!(md5_reader(InterruptedOnce(false)).unwrap(), md5_hex(b""));
        struct Broken;
        impl Read for Broken {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::ErrorKind::InvalidData.into())
            }
        }
        assert_eq!(
            md5_reader(Broken).unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
    }

    #[test]
    fn failure_after_a_successful_prefix_stays_an_error() {
        struct PrefixThenError {
            left: usize,
            calls: usize,
        }
        impl Read for PrefixThenError {
            fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
                self.calls += 1;
                if self.left == 0 {
                    return Err(io::ErrorKind::PermissionDenied.into());
                }
                let count = self.left.min(output.len()).min(17);
                output[..count].fill(0xa5);
                self.left -= count;
                Ok(count)
            }
        }
        for length in [1, 55, 56, 63, 64, 65, 127, 128, 1025] {
            let mut source = PrefixThenError {
                left: length,
                calls: 0,
            };
            assert_eq!(
                md5_reader(&mut source).unwrap_err().kind(),
                io::ErrorKind::PermissionDenied
            );
            assert_eq!(source.left, 0);
            assert!(source.calls > 1);
        }
    }

    #[test]
    fn byte_count_wraps_modulo_u64_while_preserving_buffer_offset() {
        let mut digest = Md5::new();
        digest.length = u64::MAX - 63;
        digest.update(&[0; 64]);
        assert_eq!(digest.length, 0);
        assert_eq!(digest.pending_len, 0);
        digest.update(&[1; 7]);
        assert_eq!(digest.length, 7);
        assert_eq!(digest.pending_len, 7);
    }
}
