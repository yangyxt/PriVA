//! Failure boundary for the opt-in Python-equivalent builders.
//!
//! Argument admission keeps the CLI's exit 2. Once admitted, an original
//! procedure exception is an exit-1 runtime failure. Keep its exception
//! category explicit without presenting a fabricated Python traceback.
use crate::cli::Result;
use crate::pandas_polars::NonUniqueRightKeys;
use std::error::Error;
use std::fmt;
use std::io;

#[derive(Debug)]
pub(crate) struct RuntimeError {
    category: &'static str,
    cause: Box<dyn Error>,
}

impl fmt::Display for RuntimeError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(out, "{}: {}", self.category, self.cause)
    }
}

impl RuntimeError {
    /// The Python exception class name this failure stands for.
    pub(crate) fn category_name(&self) -> &'static str {
        self.category
    }
    /// `str()` of the exception, without the category prefix.
    pub(crate) fn message(&self) -> String {
        self.cause.to_string()
    }
}

impl Error for RuntimeError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.cause.as_ref())
    }
}

pub(crate) fn category(category: &'static str, cause: Box<dyn Error>) -> Box<dyn Error> {
    if cause.is::<RuntimeError>() {
        cause
    } else {
        Box::new(RuntimeError { category, cause })
    }
}

fn runtime(cause: Box<dyn Error>) -> Box<dyn Error> {
    let name = if cause.is::<NonUniqueRightKeys>() {
        "MergeError"
    } else if let Some(error) = cause.downcast_ref::<io::Error>() {
        match error.kind() {
            io::ErrorKind::NotFound => "FileNotFoundError",
            io::ErrorKind::PermissionDenied => "PermissionError",
            io::ErrorKind::IsADirectory => "IsADirectoryError",
            io::ErrorKind::NotADirectory => "NotADirectoryError",
            io::ErrorKind::AlreadyExists => "FileExistsError",
            io::ErrorKind::BrokenPipe => "BrokenPipeError",
            // `cli::invalid` represents the original explicit ValueError.
            // An actual OS EINVAL remains OSError, identified by its errno.
            io::ErrorKind::InvalidInput if error.raw_os_error().is_none() => "ValueError",
            _ => "OSError",
        }
    } else {
        // Unexpected implementation failures must not be mislabeled as an
        // anticipated Python input exception merely to pass a comparison.
        "RuntimeError"
    };
    category(name, cause)
}

pub(crate) fn finish<T>(result: Result<T>) -> Result<T> {
    result.map_err(runtime)
}

pub(crate) fn gzip(error: io::Error) -> Box<dyn Error> {
    // flate2 1.1.10 exposes both its gzip framing errors and zio's DEFLATE
    // failure as InvalidInput. The pinned zio diagnostic distinguishes the
    // latter; Python exposes zlib.error, whose class __name__ is "error".
    let name = match error.kind() {
        io::ErrorKind::UnexpectedEof => "EOFError",
        io::ErrorKind::InvalidInput if error.to_string() == "corrupt deflate stream" => "error",
        io::ErrorKind::InvalidInput => "BadGzipFile",
        _ => return runtime(error.into()),
    };
    category(name, error.into())
}
