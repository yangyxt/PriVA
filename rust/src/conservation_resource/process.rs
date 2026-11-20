use super::*;
use std::ffi::CString;
use std::io::Read;
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

// Byte-for-byte copy of the original embedded program; retain the mature
// reader, its chromosome iteration and its bounded signal/index probes.
const PERL_READER: &str = include_str!("reader.pl");

fn which(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH").unwrap_or_else(|| "/bin:/usr/bin".into());
    std::env::split_paths(&path).find_map(|folder| {
        let path = folder.join(name);
        let c = CString::new(path.as_os_str().as_encoded_bytes()).ok()?;
        // shutil.which requires an executable non-directory; it does not
        // resolve the returned name or choose a different environment.
        (path.exists() && !path.is_dir() && unsafe { libc::access(c.as_ptr(), libc::X_OK) } == 0)
            .then_some(path)
    })
}

fn text(bytes: &[u8]) -> Result<String> {
    crate::python_io::decode_utf8(bytes)
        .map(|s| s.replace("\r\n", "\n").replace('\r', "\n"))
        .map_err(|message| error("UnicodeDecodeError", message))
}
fn command_repr(program: &Path, args: &[String]) -> String {
    format!(
        "[{}]",
        std::iter::once(repr(&program.display().to_string()))
            .chain(args.iter().map(|v| repr(v)))
            .collect::<Vec<_>>()
            .join(", ")
    )
}
fn captured(
    program: &Path,
    args: &[String],
    timeout: Option<Duration>,
    c_locale: bool,
) -> Result<(i32, String, String)> {
    let mut command = Command::new(program);
    command
        .args(args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if c_locale {
        command.env("LC_ALL", "C");
    }
    let mut child = command.spawn().map_err(|e| os_error(e, program, None))?;
    let mut stdout = child.stdout.take().unwrap();
    let mut stderr = child.stderr.take().unwrap();
    let out = thread::spawn(move || {
        let mut bytes = Vec::new();
        stdout.read_to_end(&mut bytes).map(|_| bytes)
    });
    let err = thread::spawn(move || {
        let mut bytes = Vec::new();
        stderr.read_to_end(&mut bytes).map(|_| bytes)
    });
    let started = Instant::now();
    let (status, timed_out) = loop {
        if let Some(status) = child.try_wait()? {
            break (status, false);
        }
        if timeout.is_some_and(|limit| started.elapsed() >= limit) {
            let _ = child.kill();
            break (child.wait()?, true);
        }
        thread::sleep(Duration::from_millis(5));
    };
    let out = out
        .join()
        .map_err(|_| invalid("stdout reader panicked"))??;
    let err = err
        .join()
        .map_err(|_| invalid("stderr reader panicked"))??;
    if timed_out {
        return Err(error(
            "TimeoutExpired",
            format!(
                "Command '{}' timed out after {} seconds",
                command_repr(program, args),
                timeout.unwrap().as_secs()
            ),
        ));
    }
    use std::os::unix::process::ExitStatusExt;
    let code = status
        .code()
        .unwrap_or_else(|| -status.signal().unwrap_or(1));
    Ok((code, text(&out)?, text(&err)?))
}

pub(super) fn bigwig(path: &Path, prefix: &Path) -> Result<Json> {
    let selected = prefix.join("perl");
    let executable = if is_file(&selected)? {
        Some(selected)
    } else {
        which("perl")
    };
    let executable = executable.ok_or_else(|| {
        error(
            "ValueError",
            "Install the declared perl-bio-bigfile runtime before Conservation validation",
        )
    })?;
    let (code, stdout, stderr) = captured(
        &executable,
        &["-e".into(), PERL_READER.into(), path.display().to_string()],
        Some(Duration::from_secs(120)),
        false,
    )?;
    require(
        code == 0,
        format!(
            "Cannot validate Conservation BigWig with declared Bio::DB::BigFile: {}: {}",
            path.display(),
            crate::python_text::python_strip(&stderr)
        ),
    )?;
    let value =
        json_text::parse(&stdout, false).map_err(|message| error("JSONDecodeError", message))?;
    if value.as_object().is_none() {
        let kind = match value {
            Json::Null => "NoneType",
            Json::Array(_) => "list",
            Json::Str(_) => "str",
            Json::Bool(_) => "bool",
            Json::Int(_) => "int",
            Json::Float(_) => "float",
            Json::Object(_) => unreachable!(),
        };
        return Err(error(
            "AttributeError",
            format!("'{kind}' object has no attribute 'get'"),
        ));
    }
    require(
        value
            .get("chroms")
            .and_then(Json::as_object)
            .is_some_and(|v| !v.is_empty()),
        "Invalid BigWig chromosome dictionary",
    )?;
    Ok(value)
}
pub(super) fn bsd_sum(path: &Path) -> Result<[i128; 2]> {
    let args = vec!["-r".into(), "--".into(), path.display().to_string()];
    let (code, stdout, _) = captured(Path::new("sum"), &args, None, true)?;
    if code != 0 {
        let message = if code < 0 {
            let signal = crate::python_subprocess::signal_repr(-code);
            format!(
                "Command '{}' died with {signal}.",
                command_repr(Path::new("sum"), &args)
            )
        } else {
            format!(
                "Command '{}' returned non-zero exit status {code}.",
                command_repr(Path::new("sum"), &args)
            )
        };
        return Err(error("CalledProcessError", message));
    }
    let columns: Vec<_> = stdout
        .split(crate::python_text::is_python_space)
        .filter(|s| !s.is_empty())
        .collect();
    require(
        columns.len() >= 2 && columns[..2].iter().all(|s| reference::isdigit(s)),
        "Invalid BSD checksum output",
    )?;
    let mut out = [0; 2];
    for (index, value) in columns[..2].iter().enumerate() {
        out[index] = reference::decimal(value)?.parse().map_err(|_| {
            error(
                "ValueError",
                "BSD checksum integer exceeds native receipt range",
            )
        })?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reader_program_remains_original() {
        let source = include_str!("../../../scripts/conservation_resource_check.py");
        let original = source
            .split("PERL_READER = r'''")
            .nth(1)
            .unwrap()
            .split("'''")
            .next()
            .unwrap();
        assert_eq!(PERL_READER, original);
    }
}
