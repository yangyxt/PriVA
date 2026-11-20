//! The content key of a step's results.
//!
//! A step's lineage says how its results were made; its content says what
//! they hold. A step downstream keys on the content of its inputs, so a
//! change upstream that leaves the results byte for byte the same (a new
//! comment in a function, a faster command) stops there instead of rebuilding
//! everything it reaches.
//!
//! The content is the SHA-256 of a listing of the step's data folder and its
//! output keys, one line per entry, in byte order of the path:
//!
//! ```text
//! d REL                          a folder
//! f REL SHA256                   a file, by its bytes
//! l REL input ID REL2 CONTENT    a link into an input step's folder
//! l REL own REL2                 a link inside this step's folder
//! l REL tree REL2 SHA256         a link to a file shipped in the checkout
//! l REL other IDENTITY           a link elsewhere (path size mtime)
//! k KEY <the same forms>         an output key's value
//! ```
//!
//! No path in it names the step's own folder, which holds the lineage, so the
//! same results under another lineage have the same content.
use super::{Result, fail, identity, realpath_m};
use crate::sha256::{sha256_file, sha256_hex};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// An input step's folder and its content (or lineage, when its receipt has
/// no content).
pub(crate) struct Input {
    pub(crate) id: String,
    pub(crate) folder: PathBuf,
    pub(crate) content: String,
}

pub(crate) struct Places<'a> {
    pub(crate) own: PathBuf,
    pub(crate) tree: PathBuf,
    pub(crate) inputs: &'a [Input],
}

impl Places<'_> {
    /// How a path outside the step's data walk is named: by input, by its
    /// place in this folder, by a shipped file's bytes, or by identity.
    fn name(&self, target: &Path) -> String {
        let target = realpath_m(target);
        for input in self.inputs {
            if let Ok(rel) = target.strip_prefix(&input.folder) {
                return format!("input {} {} {}", input.id, rel.display(), input.content);
            }
        }
        if let Ok(rel) = target.strip_prefix(&self.own) {
            return format!("own {}", rel.display());
        }
        if let Ok(rel) = target.strip_prefix(&self.tree) {
            if target.is_file() {
                let sum = sha256_file(&target).unwrap_or_else(|_| "unreadable".into());
                return format!("tree {} {sum}", rel.display());
            }
            return format!("tree {}", rel.display());
        }
        format!("other {}", identity(&target.to_string_lossy()))
    }
}

/// One entry of the walk; files get their sum later, on several threads.
enum Entry {
    Folder(String),
    File(String, PathBuf),
    Link(String, String),
}

fn walk(base: &Path, folder: &Path, places: &Places, out: &mut Vec<Entry>) -> Result<()> {
    let mut children: Vec<_> = std::fs::read_dir(folder)
        .map_err(|e| fail(format!("cannot list {}: {e}", folder.display())))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    children.sort_by(|a, b| {
        a.as_os_str()
            .as_encoded_bytes()
            .cmp(b.as_os_str().as_encoded_bytes())
    });
    for path in children {
        let rel = path
            .strip_prefix(base)
            .unwrap_or(&path)
            .to_string_lossy()
            .into_owned();
        let meta = std::fs::symlink_metadata(&path)
            .map_err(|e| fail(format!("cannot read {}: {e}", path.display())))?;
        if meta.file_type().is_symlink() {
            let target = std::fs::read_link(&path).unwrap_or_default();
            let target = if target.is_absolute() {
                target
            } else {
                path.parent().unwrap_or(base).join(target)
            };
            out.push(Entry::Link(rel, places.name(&target)));
        } else if meta.is_dir() {
            out.push(Entry::Folder(rel));
            walk(base, &path, places, out)?;
        } else {
            out.push(Entry::File(rel, path));
        }
    }
    Ok(())
}

/// The SHA-256 of every file, `threads` at a time, in the order given.
fn sums(files: &[&Path], threads: usize) -> Result<Vec<String>> {
    let next = AtomicUsize::new(0);
    let results = Mutex::new(vec![String::new(); files.len()]);
    let failure = Mutex::new(None);
    std::thread::scope(|scope| {
        for _ in 0..threads.clamp(1, files.len().max(1)) {
            scope.spawn(|| {
                loop {
                    let at = next.fetch_add(1, Ordering::Relaxed);
                    if at >= files.len() {
                        break;
                    }
                    match sha256_file(files[at]) {
                        Ok(sum) => results.lock().unwrap()[at] = sum,
                        Err(e) => {
                            *failure.lock().unwrap() =
                                Some(format!("cannot read {}: {e}", files[at].display()));
                        }
                    }
                }
            });
        }
    });
    if let Some(message) = failure.into_inner().unwrap() {
        return Err(fail(message));
    }
    Ok(results.into_inner().unwrap())
}

/// The listing described above, for the data folder `data` and the output
/// `keys` (key, value) in the order given.
pub(crate) fn listing(
    data: &Path,
    keys: &[(String, String)],
    places: &Places,
    threads: usize,
) -> Result<String> {
    let mut entries = Vec::new();
    if data.is_dir() {
        walk(data, data, places, &mut entries)?;
    }
    let files: Vec<&Path> = entries
        .iter()
        .filter_map(|e| match e {
            Entry::File(_, path) => Some(path.as_path()),
            _ => None,
        })
        .collect();
    let mut sums = sums(&files, threads)?.into_iter();
    let mut text = String::new();
    for entry in &entries {
        match entry {
            Entry::Folder(rel) => text.push_str(&format!("d {rel}\n")),
            Entry::File(rel, _) => {
                text.push_str(&format!("f {rel} {}\n", sums.next().unwrap_or_default()))
            }
            Entry::Link(rel, name) => text.push_str(&format!("l {rel} {name}\n")),
        }
    }
    for (key, value) in keys {
        text.push_str(&format!("k {key} {}\n", places.name(Path::new(value))));
    }
    Ok(text)
}

pub(crate) fn content(text: &str) -> String {
    sha256_hex(text.as_bytes())
}
