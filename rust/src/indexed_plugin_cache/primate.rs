//! The publisher's physical-line-12/header/sort procedure, streamed in order.
use super::*;
use crate::file_lock::{ExclusiveLock, Wait};
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::os::fd::AsRawFd;
use std::process::{Child, ChildStdin};

fn binary_line(reader: &mut impl BufRead, limit: usize) -> Result<Vec<u8>> {
    let mut result = Vec::new();
    reader
        .take(limit as u64)
        .read_until(b'\n', &mut result)
        .map_err(|e| gzip_error(e.into()))?;
    Ok(result)
}
fn raw_format(path: &Path) -> Result<bool> {
    let mut reader = BufReader::new(gzip(path)?);
    let mut line = Vec::new();
    for _ in 0..12 {
        line = binary_line(&mut reader, PROBE)?;
    }
    let head = line.split(|b| *b == b'\t').next().unwrap_or(&[]);
    let lower: Vec<u8> = head.iter().map(u8::to_ascii_lowercase).collect();
    let trimmed = lower.trim_ascii();
    if [
        b"chr".as_slice(),
        b"chrom".as_slice(),
        b"chromosome".as_slice(),
    ]
    .contains(&trimmed)
    {
        return Ok(true);
    }
    first_record(path, default_settings("primateai"))?;
    Ok(false)
}
fn sha(path: &Path) -> Result<String> {
    crate::sha256::sha256_file(path).map_err(|e| os_error(e, path, None))
}
fn read_json(path: &Path) -> Result<Json> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    let text = crate::python_io::decode_utf8(&bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    json_text::parse(text, false).map_err(|e| error("JSONDecodeError", e))
}
fn dict_get<'a>(doc: &'a Json, key: &str) -> Result<&'a Json> {
    if doc.as_object().is_none() {
        let kind = match doc {
            Json::Array(_) => "list",
            Json::Str(_) => "str",
            Json::Null => "NoneType",
            Json::Int(_) => "int",
            Json::Float(_) => "float",
            Json::Bool(_) => "bool",
            Json::Object(_) => unreachable!(),
        };
        return Err(error(
            "AttributeError",
            format!("'{kind}' object has no attribute 'get'"),
        ));
    }
    Ok(doc.get(key).unwrap_or(&Json::Null))
}
fn mkdir(path: &Path) -> Result<()> {
    fs::create_dir_all(path).map_err(|e| os_error(e, path, None))
}
fn write_json(path: &Path, value: &Json) -> Result<()> {
    let payload = json_text::dumps_with_options(
        value,
        json_text::DumpsOptions {
            sort_keys: true,
            indent: Some(2),
            ..Default::default()
        },
    )
    .map_err(invalid)?
        + "\n";
    fs::write(path, payload).map_err(|e| os_error(e, path, None))
}
/// Python's finally blocks kill and reap a still-running child even after a
/// stream failure. The parent closes stdin before a normal wait; child stderr
/// remains inherited rather than hidden inside a new logging policy.
struct Reaped(Child, Option<BufWriter<ChildStdin>>, usize);
impl Reaped {
    fn new(mut child: Child) -> Self {
        let pipe = child.stdin.take().unwrap();
        // CPython 3.12 io.open uses the descriptor's positive st_blksize,
        // falling back to DEFAULT_BUFFER_SIZE. A pipe can be 4096 bytes even
        // though the default constant is 8192. Preserve that failure boundary.
        let mut info = std::mem::MaybeUninit::<libc::stat>::uninit();
        let capacity = if unsafe { libc::fstat(pipe.as_raw_fd(), info.as_mut_ptr()) } == 0 {
            let block = unsafe { info.assume_init() }.st_blksize;
            if block > 1 { block as usize } else { 8192 }
        } else {
            8192
        };
        // Python buffers an exactly full write; Rust directly writes a slice
        // at capacity. One extra byte plus an explicit flush keeps Python's
        // boundary while using the standard buffered writer.
        let input = BufWriter::with_capacity(capacity + 1, pipe);
        Self(child, Some(input), capacity)
    }
    fn write_input(&mut self, bytes: &[u8]) -> Result<()> {
        let input = self.1.as_mut().unwrap();
        if bytes.len() > self.2 - input.buffer().len() {
            input.flush().map_err(pipe_error)?;
        }
        input.write_all(bytes).map_err(pipe_error)
    }
    fn close_input(&mut self) -> Result<()> {
        let mut input = self.1.take().unwrap();
        let outcome = input.flush().map_err(pipe_error);
        // BufferedWriter.close closes the raw descriptor even after a failed
        // flush. Do not let Rust's Drop retry that failed buffered write.
        let (pipe, _) = input.into_parts();
        drop(pipe);
        outcome
    }
}
fn pipe_error(cause: io::Error) -> Box<dyn std::error::Error> {
    if cause.kind() == io::ErrorKind::BrokenPipe {
        error("BrokenPipeError", "[Errno 32] Broken pipe")
    } else {
        cause.into()
    }
}
impl Drop for Reaped {
    fn drop(&mut self) {
        if matches!(self.0.try_wait(), Ok(None)) {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
}
fn stream(source: &Path, output: &Path, threads: &str, tmp: &Path) -> Result<()> {
    let body = output.parent().unwrap().join("sorted.body");
    let sorted = File::create(&body).map_err(|e| os_error(e, &body, None))?;
    let mut sorter = Reaped::new(
        Command::new("sort")
            .args([
                "-k1,1",
                "-k",
                "2,2n",
                "-V",
                &format!("--parallel={threads}"),
                "--buffer-size=256M",
                "-T",
            ])
            .arg(tmp)
            .stdin(Stdio::piped())
            .stdout(sorted)
            .spawn()
            .map_err(|e| os_error(e, Path::new("sort"), None))?,
    );
    let mut header = Vec::new();
    {
        let mut reader = BufReader::new(gzip(source)?);
        let mut number = 0_usize;
        let mut count = 0_usize;
        loop {
            let mut line = binary_line(&mut reader, usize::MAX)?;
            if line.is_empty() {
                break;
            }
            number += 1;
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            if number == 12 {
                line.insert(0, b'#');
            }
            if line.is_empty() {
                continue;
            }
            count += 1;
            line.push(b'\n');
            if count < 12 {
                header.extend_from_slice(&line);
            } else {
                sorter.write_input(&line)?;
            }
        }
    }
    sorter.close_input()?;
    require(sorter.0.wait()?.success(), "PrimateAI sort failed")?;
    drop(sorter);
    let compressed = File::create(output).map_err(|e| os_error(e, output, None))?;
    let mut compressor = Reaped::new(
        Command::new("bgzip")
            .args(["-@", threads, "-c"])
            .stdin(Stdio::piped())
            .stdout(compressed)
            .spawn()
            .map_err(|e| os_error(e, Path::new("bgzip"), None))?,
    );
    compressor.write_input(&header)?;
    let mut source = file(&body)?;
    let mut block = vec![0_u8; PROBE];
    loop {
        let count = source.read(&mut block)?;
        if count == 0 {
            break;
        }
        compressor.write_input(&block[..count])?;
    }
    compressor.close_input()?;
    require(
        compressor.0.wait()?.success(),
        "PrimateAI BGZF compression failed",
    )?;
    drop(compressor);
    fs::remove_file(&body).map_err(|e| os_error(e, &body, None))?;
    Ok(())
}

pub(super) fn format(
    input: &Path,
    cache: &Path,
    assembly: &str,
    threads: &str,
    tmp: Option<&Path>,
) -> Result<Json> {
    let source = PathBuf::from(absolute(input)?);
    require(
        is_file(&source)? && metadata(&source)?.len() != 0,
        format!(
            "Missing or empty selected PrimateAI source: {}",
            source.display()
        ),
    )?;
    if !raw_format(&source)? {
        return ensure_index(&source, "primateai");
    }
    let before = identity(&source)?;
    let source_identity = object([
        ("path", Json::Str(source.display().to_string())),
        ("sha256", Json::Str(sha(&source)?)),
        ("size", Json::Int(before.2 as i128)),
        ("assembly", Json::Str(assembly.into())),
    ]);
    let mut hash = crate::sha256::Sha256::new();
    hash.update(
        json_text::dumps(&source_identity, true)
            .map_err(invalid)?
            .as_bytes(),
    );
    let key = hash.finish_hex()[..20].to_owned();
    let parent = PathBuf::from(absolute(cache)?).join(assembly);
    mkdir(&parent)?;
    let destination = parent.join(&key);
    let filename = source.file_name().unwrap().to_string_lossy();
    let name = filename.strip_suffix(".gz").unwrap_or(&filename);
    let name = format!(
        "{}.sorted.tsv.bgz",
        name.strip_suffix(".tsv").unwrap_or(name)
    );
    let lock_path = parent.join(format!(".{key}.lock"));
    let lock = OpenOptions::new()
        .append(true)
        .create(true)
        .open(&lock_path)
        .map_err(|e| os_error(e, &lock_path, None))?;
    let _lock = ExclusiveLock::acquire(lock, Wait::Blocking)?;
    if exists(&destination)? {
        let receipt = read_json(&destination.join("source.ready.json"))?;
        let output = destination.join(&name);
        // Preserve Python's short-circuit order: a conflicting source/algorithm
        // never opens the output merely to compute its hash.
        require(
            json_text::py_equal(dict_get(&receipt, "source")?, &source_identity)
                && dict_get(&receipt, "algorithm")?.as_str() == Some("upstream-primate-format-v1")
                && dict_get(&receipt, "output_sha256")?.as_str() == Some(sha(&output)?.as_str()),
            format!(
                "Conflicting managed PrimateAI cache: {}",
                destination.display()
            ),
        )?;
        let result = ensure_index(&output, "primateai")?;
        require(
            dict_get(&receipt, "index_sha256")?.as_str()
                == Some(sha(&suffix(&output, ".tbi"))?.as_str()),
            format!(
                "Conflicting managed PrimateAI index: {}.tbi",
                output.display()
            ),
        )?;
        require(
            identity(&source)? == before,
            "Selected PrimateAI source changed during cache reuse",
        )?;
        return Ok(result);
    }
    let temporary = TemporaryDirectory::new(&parent, &format!(".{key}.build-"))?;
    let outcome = (|| {
        let output = temporary.0.join(&name);
        let temp_sort = tmp.unwrap_or(&temporary.0);
        mkdir(temp_sort)?;
        stream(&source, &output, threads, temp_sort)?;
        ensure_index(&output, "primateai")?;
        require(
            identity(&source)? == before,
            "Selected PrimateAI source changed during formatting",
        )?;
        write_json(
            &temporary.0.join("source.ready.json"),
            &object([
                ("source", source_identity),
                ("algorithm", Json::Str("upstream-primate-format-v1".into())),
                ("output_sha256", Json::Str(sha(&output)?)),
                ("index_sha256", Json::Str(sha(&suffix(&output, ".tbi"))?)),
                ("scientific_qualified", Json::Bool(false)),
            ]),
        )?;
        fs::rename(&temporary.0, &destination)
            .map_err(|e| os_error(e, &temporary.0, Some(&destination)))?;
        Ok(object([
            (
                "path",
                Json::Str(destination.join(&name).display().to_string()),
            ),
            (
                "index",
                Json::Str(
                    destination
                        .join(format!("{name}.tbi"))
                        .display()
                        .to_string(),
                ),
            ),
            ("formatted", Json::Bool(true)),
        ]))
    })();
    temporary.finish(outcome)
}
