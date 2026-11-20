//! `gzip.GzipFile(filename="", fileobj=..., mode="wb", mtime=0)` byte for
//! byte: the header GzipFile writes (no name, mtime 0, XFL from the level, OS
//! 255), the raw DEFLATE stream of the same zlib Python's `zlib` module uses
//! (linked statically from the build environment, as Expat and libyaml are)
//! and the CRC-32/ISIZE trailer.
//!
//! zlib's output does not depend on how the input is split between `deflate`
//! calls without a flush, so one pass over the whole text equals GzipFile's
//! buffered writes. A flush does change the stream: a `TextIOWrapper` over
//! the GzipFile flushes it (`Z_SYNC_FLUSH`) when it closes, before GzipFile's
//! own `Z_FINISH`, and `sync_flush` reproduces that.
use std::io;
use std::os::raw::{c_char, c_int, c_uint, c_ulong, c_void};

#[repr(C)]
struct ZStream {
    next_in: *const u8,
    avail_in: c_uint,
    total_in: c_ulong,
    next_out: *mut u8,
    avail_out: c_uint,
    total_out: c_ulong,
    msg: *const c_char,
    state: *mut c_void,
    zalloc: *const c_void,
    zfree: *const c_void,
    opaque: *mut c_void,
    data_type: c_int,
    adler: c_ulong,
    reserved: c_ulong,
}

unsafe extern "C" {
    fn zlibVersion() -> *const c_char;
    fn deflateInit2_(
        stream: *mut ZStream,
        level: c_int,
        method: c_int,
        window_bits: c_int,
        mem_level: c_int,
        strategy: c_int,
        version: *const c_char,
        stream_size: c_int,
    ) -> c_int;
    fn deflate(stream: *mut ZStream, flush: c_int) -> c_int;
    fn deflateEnd(stream: *mut ZStream) -> c_int;
    fn crc32(crc: c_ulong, buffer: *const u8, length: c_uint) -> c_ulong;
}

const Z_OK: c_int = 0;
const Z_STREAM_END: c_int = 1;
const Z_BUF_ERROR: c_int = -5;
const Z_NO_FLUSH: c_int = 0;
const Z_SYNC_FLUSH: c_int = 2;
const Z_FINISH: c_int = 4;
const Z_DEFLATED: c_int = 8;
const CHUNK: usize = 1 << 20;

fn failure(what: &str, code: c_int) -> io::Error {
    io::Error::other(format!("zlib {what} failed with status {code}"))
}

struct Deflater {
    stream: Box<ZStream>,
    buffer: Vec<u8>,
}

impl Deflater {
    fn new(level: u32) -> io::Result<Self> {
        // SAFETY: a zeroed z_stream (null allocators) is zlib's documented
        // initial state; the boxed stream keeps the address zlib records.
        let mut stream: Box<ZStream> = Box::new(unsafe { std::mem::zeroed() });
        let status = unsafe {
            deflateInit2_(
                &mut *stream,
                level as c_int,
                Z_DEFLATED,
                -15,
                8,
                0,
                zlibVersion(),
                std::mem::size_of::<ZStream>() as c_int,
            )
        };
        if status != Z_OK {
            return Err(failure("deflateInit2", status));
        }
        Ok(Self {
            stream,
            buffer: vec![0; CHUNK],
        })
    }

    /// `deflate` over `input` with `flush` until zlib has nothing more to
    /// write for it, appending to `out`.
    fn run(&mut self, input: &[u8], flush: c_int, out: &mut Vec<u8>) -> io::Result<()> {
        self.stream.next_in = input.as_ptr();
        self.stream.avail_in = input.len() as c_uint;
        loop {
            self.stream.next_out = self.buffer.as_mut_ptr();
            self.stream.avail_out = CHUNK as c_uint;
            // SAFETY: next_in/next_out point into live buffers of the stated sizes.
            let status = unsafe { deflate(&mut *self.stream, flush) };
            out.extend_from_slice(&self.buffer[..CHUNK - self.stream.avail_out as usize]);
            match status {
                Z_STREAM_END => return Ok(()),
                Z_OK | Z_BUF_ERROR => {}
                code => return Err(failure("deflate", code)),
            }
            if flush != Z_FINISH && self.stream.avail_in == 0 && self.stream.avail_out != 0 {
                return Ok(());
            }
        }
    }
}

impl Drop for Deflater {
    fn drop(&mut self) {
        // SAFETY: the stream was initialized by deflateInit2_.
        unsafe {
            deflateEnd(&mut *self.stream);
        }
    }
}

/// The complete gzip member GzipFile writes for `body` at `level` (0 to 9),
/// with one `Z_SYNC_FLUSH` before the end when `sync_flush` is set.
pub(crate) fn compress(body: &[u8], level: u32, sync_flush: bool) -> io::Result<Vec<u8>> {
    compress_named(body, level, sync_flush, b"", 0)
}

/// `compress` with the name and modification time `gzip.open(path, "wb")`
/// stores in the header; an empty `name` stores none.
pub(crate) fn compress_named(
    body: &[u8],
    level: u32,
    sync_flush: bool,
    name: &[u8],
    mtime: u32,
) -> io::Result<Vec<u8>> {
    if level > 9 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "gzip level must be in 0..=9",
        ));
    }
    let xfl = match level {
        9 => 2,
        1 => 4,
        _ => 0,
    };
    let mut out = vec![0x1f, 0x8b, 8, if name.is_empty() { 0 } else { 8 }];
    out.extend_from_slice(&mtime.to_le_bytes());
    out.extend_from_slice(&[xfl, 0xff]);
    if !name.is_empty() {
        out.extend_from_slice(name);
        out.push(0);
    }
    let mut deflater = Deflater::new(level)?;
    let mut crc: c_ulong = 0;
    for part in body.chunks(CHUNK) {
        deflater.run(part, Z_NO_FLUSH, &mut out)?;
        // SAFETY: part is a live slice of the stated length.
        crc = unsafe { crc32(crc, part.as_ptr(), part.len() as c_uint) };
    }
    if sync_flush {
        deflater.run(&[], Z_SYNC_FLUSH, &mut out)?;
    }
    deflater.run(&[], Z_FINISH, &mut out)?;
    out.extend_from_slice(&(crc as u32).to_le_bytes());
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn header_and_round_trip() {
        let body = b"{\"a\":1}\n".repeat(500_000);
        for sync_flush in [false, true] {
            let member = compress(&body, 9, sync_flush).unwrap();
            assert_eq!(&member[..10], &[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 2, 0xff]);
            let mut decoded = Vec::new();
            flate2::read::GzDecoder::new(&member[..])
                .read_to_end(&mut decoded)
                .unwrap();
            assert_eq!(decoded, body);
        }
        assert_eq!(compress(b"", 9, false).unwrap().len(), 20);
    }
}
