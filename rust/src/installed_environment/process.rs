//! Captured subprocesses retain binary Git output and text-manager decode order.
use super::*;
use std::io::Read;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};
pub(super) fn run(
    mut command: Command,
    timeout: Option<Duration>,
) -> Result<(i32, Vec<u8>, Vec<u8>)> {
    let program = PathBuf::from(command.get_program());
    command.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = command.spawn().map_err(|e| os_error(e, &program, None))?;
    let mut stdout = child.stdout.take().unwrap();
    let mut stderr = child.stderr.take().unwrap();
    let out = std::thread::spawn(move || {
        let mut b = vec![];
        stdout.read_to_end(&mut b).map(|_| b)
    });
    let err = std::thread::spawn(move || {
        let mut b = vec![];
        stderr.read_to_end(&mut b).map(|_| b)
    });
    let start = Instant::now();
    let (status, timed_out) = loop {
        if let Some(status) = child.try_wait()? {
            break (status, false);
        }
        if timeout.is_some_and(|limit| start.elapsed() >= limit) {
            let _ = child.kill();
            break (child.wait()?, true);
        }
        std::thread::sleep(Duration::from_millis(5));
    };
    let stdout = out
        .join()
        .map_err(|_| invalid("stdout reader panicked"))??;
    let stderr = err
        .join()
        .map_err(|_| invalid("stderr reader panicked"))??;
    if timed_out {
        return Err(error(
            "TimeoutExpired",
            "Installed-package Git provenance command timed out",
        ));
    }
    use std::os::unix::process::ExitStatusExt;
    Ok((
        status
            .code()
            .unwrap_or_else(|| -status.signal().unwrap_or(1)),
        stdout,
        stderr,
    ))
}
pub(super) fn git(directory: &Path, args: &[&std::ffi::OsStr]) -> Result<Vec<u8>> {
    let mut command = Command::new("git");
    command
        .arg("-C")
        .arg(directory)
        .arg("--literal-pathspecs")
        .args(args);
    command
        .env("GIT_OPTIONAL_LOCKS", "0")
        .env("GIT_TERMINAL_PROMPT", "0");
    for name in [
        "GIT_DIR",
        "GIT_COMMON_DIR",
        "GIT_WORK_TREE",
        "GIT_INDEX_FILE",
    ] {
        command.env_remove(name);
    }
    let (code, out, _) = run(command, Some(Duration::from_secs(30)))?;
    require(
        code == 0,
        "Installed-package link has no committed git-annex provenance",
    )?;
    Ok(out)
}
pub(super) fn manager(manager: &str, selection: &[String]) -> Result<String> {
    let mut command = Command::new(manager);
    command.args(["env", "export"]).args(selection);
    let (code, stdout, stderr) = run(command, None)?;
    // subprocess.run(text=True) decodes both pipes before return-code handling.
    let stdout = text(&stdout)?;
    let _ = text(&stderr)?;
    require(
        code == 0,
        "The package manager could not export the installed environment",
    )?;
    Ok(stdout)
}
