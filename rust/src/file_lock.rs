//! Exclusive locks used by installer helpers on Unix.
//!
//! Callers open the file or directory themselves: Python's append-file,
//! read-only directory and O_NOFOLLOW publication routes are distinct and must
//! retain their own opening checks. The owned descriptor keeps the lock alive;
//! dropping it releases the lock without unlinking or replacing the lock inode.
//! Rust 1.90 uses flock on Unix, matching the pinned Python fcntl.flock calls.
use std::fs::File;
use std::io;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Wait {
    Blocking,
    Nonblocking,
}

/// An exclusively locked, owned descriptor. It is deliberately not Clone and
/// does not expose its File: reacquiring a duplicated descriptor has
/// platform-dependent behavior, and the lock must outlive the protected work.
pub(crate) struct ExclusiveLock {
    file: File,
}

impl ExclusiveLock {
    pub(crate) fn acquire(file: File, wait: Wait) -> io::Result<Self> {
        match wait {
            Wait::Blocking => file.lock()?,
            Wait::Nonblocking => file.try_lock()?,
        }
        Ok(Self { file })
    }

    /// Explicitly report an unlock failure, then close the descriptor. Ordinary
    /// scope exit also releases the lock through File's descriptor lifetime.
    pub(crate) fn release(self) -> io::Result<()> {
        self.file.unlock()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, OpenOptions};
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    #[test]
    fn lock_lifetime_retains_inode_and_contents() {
        let path = std::env::temp_dir().join(format!(
            "priva-file-lock-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::write(&path, b"preserved").unwrap();
        let open = || OpenOptions::new().append(true).open(&path).unwrap();
        let first = ExclusiveLock::acquire(open(), Wait::Blocking).unwrap();
        let blocked = ExclusiveLock::acquire(open(), Wait::Nonblocking);
        assert!(matches!(blocked, Err(error) if error.kind() == io::ErrorKind::WouldBlock));
        first.release().unwrap();
        drop(ExclusiveLock::acquire(open(), Wait::Nonblocking).unwrap());
        let after_drop = ExclusiveLock::acquire(open(), Wait::Nonblocking).unwrap();
        after_drop.release().unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"preserved");
        fs::remove_file(path).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn directory_lock_uses_existing_directory_descriptor() {
        let path = std::env::temp_dir().join(format!(
            "priva-directory-lock-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        let lock = ExclusiveLock::acquire(File::open(&path).unwrap(), Wait::Nonblocking).unwrap();
        let blocked = ExclusiveLock::acquire(File::open(&path).unwrap(), Wait::Nonblocking);
        assert!(matches!(blocked, Err(error) if error.kind() == io::ErrorKind::WouldBlock));
        drop(lock);
        assert_eq!(fs::read_dir(&path).unwrap().count(), 0);
        drop(ExclusiveLock::acquire(File::open(&path).unwrap(), Wait::Blocking).unwrap());
        fs::remove_dir(path).unwrap();
    }
}
