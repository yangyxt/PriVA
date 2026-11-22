//! Streaming gzip for JSON/TSV builder outputs, using the pinned flate2 backend.
//! Defaults are level 9, mtime 0 and no filename. Callers preserving Python or
//! pandas metadata supply their selected timestamp and filename header bytes.
//! Text encoding, newline handling and publication remain caller-owned. Gzip
//! compatibility means identical decompressed bytes, not identical DEFLATE bytes.
use flate2::write::GzEncoder;
use flate2::{Compression, GzBuilder};
use std::io::{self, Write};

#[derive(Clone, Debug)]
pub(crate) struct GzipOptions {
    pub(crate) level: u32,
    pub(crate) mtime: u32,
    /// Already-selected header bytes; no basename or text-encoding conversion.
    /// An empty name is omitted, as with Python GzipFile. NUL is invalid.
    pub(crate) filename: Option<Vec<u8>>,
}

impl Default for GzipOptions {
    fn default() -> Self {
        Self {
            level: 9,
            mtime: 0,
            filename: None,
        }
    }
}

/// Success requires calling `finish`; dropping a writer cannot report errors.
/// A failed write/flush is retained so a later finish cannot certify partial data.
#[must_use = "call finish and check its result before publishing the output"]
pub(crate) struct GzipWriter<W: Write> {
    encoder: GzEncoder<HeaderWriter<W>>,
    failure: Option<(io::ErrorKind, String)>,
}

struct HeaderWriter<W> {
    writer: W,
    prefix_written: usize,
    level_zero: bool,
}

impl<W: Write> Write for HeaderWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        // flate2 1.1.10 uses XFL=4 for level 0 as well as level 1. GzipFile
        // uses XFL=0 for level 0. Correct only that header byte, preserving the
        // DEFLATE level/payload and the underlying writer's short-write rules.
        if self.level_zero && self.prefix_written < 9 {
            let count = bytes.len().min(9 - self.prefix_written);
            let mut prefix = [0; 9];
            prefix[..count].copy_from_slice(&bytes[..count]);
            if self.prefix_written + count == 9 {
                prefix[8 - self.prefix_written] = 0;
            }
            let written = self.writer.write(&prefix[..count])?;
            self.prefix_written += written;
            Ok(written)
        } else {
            self.writer.write(bytes)
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

impl<W: Write> GzipWriter<W> {
    pub(crate) fn new(writer: W, options: GzipOptions) -> io::Result<Self> {
        if options.level > 9 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "gzip level must be in 0..=9",
            ));
        }
        if options
            .filename
            .as_ref()
            .is_some_and(|name| name.contains(&0))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "gzip filename contains NUL",
            ));
        }
        // GzipFile emits OS=255 independently of the host; keep this explicit.
        let mut builder = GzBuilder::new().mtime(options.mtime).operating_system(255);
        if let Some(name) = options.filename.filter(|name| !name.is_empty()) {
            builder = builder.filename(name);
        }
        Ok(Self {
            encoder: builder.write(
                HeaderWriter {
                    writer,
                    prefix_written: 0,
                    level_zero: options.level == 0,
                },
                Compression::new(options.level),
            ),
            failure: None,
        })
    }

    fn previous_error(&self) -> Option<io::Error> {
        self.failure
            .as_ref()
            .map(|(kind, message)| io::Error::new(*kind, message.clone()))
    }

    fn remember<T>(&mut self, result: &io::Result<T>) {
        if let Err(error) = result {
            // Write::write_all retries Interrupted. It does not mean data failed.
            if error.kind() != io::ErrorKind::Interrupted && self.failure.is_none() {
                self.failure = Some((error.kind(), error.to_string()));
            }
        }
    }

    /// Complete the member and trailer, then flush the underlying writer.
    /// No extra DEFLATE sync-flush is inserted here. Interrupted I/O is retried.
    /// This does not fsync, close a file descriptor, or publish/rename a file.
    pub(crate) fn finish(mut self) -> io::Result<W> {
        if let Some(error) = self.previous_error() {
            return Err(error);
        }
        loop {
            match self.encoder.try_finish() {
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                result => {
                    result?;
                    break;
                }
            }
        }
        // flate2 explicitly permits finish after successful try_finish.
        let mut writer = self.encoder.finish()?.writer;
        loop {
            match writer.flush() {
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Ok(()) => return Ok(writer),
                Err(error) => return Err(error),
            }
        }
    }
}

impl<W: Write> Write for GzipWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if let Some(error) = self.previous_error() {
            return Err(error);
        }
        let result = self.encoder.write(bytes);
        self.remember(&result);
        result
    }

    /// A requested sync-flush can change compressed bytes, never the payload.
    fn flush(&mut self) -> io::Result<()> {
        if let Some(error) = self.previous_error() {
            return Err(error);
        }
        let result = self.encoder.flush();
        self.remember(&result);
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::io::Read;
    use std::rc::Rc;

    fn encoded(data: &[u8], options: GzipOptions) -> Vec<u8> {
        let mut writer = GzipWriter::new(Vec::new(), options).unwrap();
        writer.write_all(data).unwrap();
        writer.finish().unwrap()
    }

    #[derive(Default, Debug)]
    struct Sink {
        bytes: Vec<u8>,
        fault: Option<(usize, io::ErrorKind)>,
        fired: bool,
        flush_fault: Option<io::ErrorKind>,
        zero: bool,
        max_write: usize,
    }

    #[derive(Clone, Default, Debug)]
    struct Shared(Rc<RefCell<Sink>>);

    impl Write for Shared {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            let mut state = self.0.borrow_mut();
            if state.zero {
                return Ok(0);
            }
            let mut count = bytes.len();
            if state.max_write > 0 {
                count = count.min(state.max_write);
            }
            if let Some((at, kind)) = state.fault.filter(|_| !state.fired) {
                if state.bytes.len() == at {
                    state.fired = true;
                    return Err(kind.into());
                }
                count = count.min(at - state.bytes.len());
            }
            state.bytes.extend_from_slice(&bytes[..count]);
            Ok(count)
        }
        fn flush(&mut self) -> io::Result<()> {
            self.0
                .borrow_mut()
                .flush_fault
                .take()
                .map_or(Ok(()), |kind| Err(kind.into()))
        }
    }

    #[test]
    fn default_header_and_payload_are_deterministic() {
        for data in [b"".as_slice(), b"a\0\xff\n", "gene\tβ\r\n".as_bytes()] {
            let bytes = encoded(data, GzipOptions::default());
            assert_eq!(&bytes[..10], &[31, 139, 8, 0, 0, 0, 0, 0, 2, 255]);
            assert_eq!(bytes, encoded(data, GzipOptions::default()));
            let mut restored = Vec::new();
            flate2::read::GzDecoder::new(bytes.as_slice())
                .read_to_end(&mut restored)
                .unwrap();
            assert_eq!(restored, data);
        }
    }

    #[test]
    fn explicit_header_options_and_empty_name() {
        let options = GzipOptions {
            level: 1,
            mtime: 0x01020304,
            filename: Some(b"caf\xe9.tsv".to_vec()),
        };
        let bytes = encoded(b"payload", options);
        assert_eq!(&bytes[..10], &[31, 139, 8, 8, 4, 3, 2, 1, 4, 255]);
        assert_eq!(&bytes[10..19], b"caf\xe9.tsv\0");
        let bytes = encoded(
            b"",
            GzipOptions {
                level: 0,
                filename: Some(Vec::new()),
                ..GzipOptions::default()
            },
        );
        assert_eq!(bytes[3], 0);
        assert_eq!(bytes[8], 0);
    }

    #[test]
    fn invalid_options_fail_without_writing() {
        for options in [
            GzipOptions {
                level: 10,
                ..GzipOptions::default()
            },
            GzipOptions {
                filename: Some(b"a\0b".to_vec()),
                ..GzipOptions::default()
            },
        ] {
            let sink = Shared::default();
            let error = match GzipWriter::new(sink.clone(), options) {
                Ok(_) => panic!("invalid options accepted"),
                Err(error) => error,
            };
            assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
            assert!(sink.0.borrow().bytes.is_empty());
        }
    }

    #[test]
    fn short_writes_and_interruptions_complete_exactly() {
        let payload: Vec<u8> = (0..4096).map(|n| (n % 251) as u8).collect();
        let expected = encoded(&payload, GzipOptions::default());
        for at in [0, 3, 10, expected.len() - 8, expected.len() - 1] {
            let sink = Shared::default();
            sink.0.borrow_mut().fault = Some((at, io::ErrorKind::Interrupted));
            sink.0.borrow_mut().max_write = 3;
            sink.0.borrow_mut().flush_fault = Some(io::ErrorKind::Interrupted);
            let mut writer = GzipWriter::new(sink.clone(), GzipOptions::default()).unwrap();
            writer.write_all(&payload).unwrap();
            writer.finish().unwrap();
            assert!(sink.0.borrow().fired);
            assert_eq!(sink.0.borrow().bytes, expected);
        }
    }

    #[test]
    fn header_body_and_trailer_failures_remain_errors() {
        let payload = vec![b'x'; 4096];
        let length = encoded(&payload, GzipOptions::default()).len();
        for at in [0, 3, 10, length - 8, length - 1] {
            let sink = Shared::default();
            sink.0.borrow_mut().fault = Some((at, io::ErrorKind::BrokenPipe));
            let mut writer = GzipWriter::new(sink.clone(), GzipOptions::default()).unwrap();
            let write = writer.write_all(&payload);
            let finish = writer.finish();
            assert!(sink.0.borrow().fired);
            assert!(finish.is_err());
            assert_eq!(finish.unwrap_err().kind(), io::ErrorKind::BrokenPipe);
            if let Err(error) = write {
                assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
            }
        }
    }

    #[test]
    fn finish_reports_underlying_flush_failure() {
        let sink = Shared::default();
        sink.0.borrow_mut().flush_fault = Some(io::ErrorKind::PermissionDenied);
        let mut writer = GzipWriter::new(sink, GzipOptions::default()).unwrap();
        writer.write_all(b"payload").unwrap();
        assert_eq!(
            writer.finish().unwrap_err().kind(),
            io::ErrorKind::PermissionDenied
        );
    }

    #[test]
    fn zero_progress_and_prior_flush_failures_cannot_finish_successfully() {
        let sink = Shared::default();
        sink.0.borrow_mut().zero = true;
        let mut writer = GzipWriter::new(sink, GzipOptions::default()).unwrap();
        assert_eq!(
            writer.write_all(b"payload").unwrap_err().kind(),
            io::ErrorKind::WriteZero
        );
        assert_eq!(
            writer.finish().unwrap_err().kind(),
            io::ErrorKind::WriteZero
        );
        let sink = Shared::default();
        sink.0.borrow_mut().flush_fault = Some(io::ErrorKind::BrokenPipe);
        let mut writer = GzipWriter::new(sink, GzipOptions::default()).unwrap();
        writer.write_all(b"payload").unwrap();
        assert_eq!(
            writer.flush().unwrap_err().kind(),
            io::ErrorKind::BrokenPipe
        );
        assert_eq!(
            writer.finish().unwrap_err().kind(),
            io::ErrorKind::BrokenPipe
        );
    }
}
