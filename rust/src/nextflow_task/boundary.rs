//! Run a Nextflow task with only its task directory writable on the host.
//!
//!   priva task-boundary check --task-dir /absolute/task
//!   priva task-boundary exec --task-dir . -- /bin/bash -e .command.sh
//!
//! This opt-in Linux Bubblewrap view does not qualify annotation, scoring or
//! family results. It never retries a denied write outside the task
//! directory. Nextflow's driver, scheduler wrapper and publication stay
//! outside this view.
use std::fs;
use std::os::unix::fs::{FileTypeExt, MetadataExt};
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};

struct Refusal(String);

fn refuse<T>(message: impl Into<String>) -> Result<T, Refusal> {
    Err(Refusal(message.into()))
}

fn os(error: std::io::Error) -> Refusal {
    Refusal(error.to_string())
}

/// Check the real Nextflow controls and reject pre-existing writable
/// aliases. A read-only cache bind does not protect a hard link reached
/// through a different writable bind. Walk only this task's physical
/// directories; never follow staged symlinks.
fn check_task(task_dir: &Path) -> Result<PathBuf, Refusal> {
    if fs::symlink_metadata(task_dir).is_ok_and(|m| m.file_type().is_symlink()) {
        return refuse("Task directory must not be a symlink");
    }
    let task = fs::canonicalize(task_dir).map_err(os)?;
    let current = std::env::current_dir()
        .and_then(fs::canonicalize)
        .map_err(os)?;
    if !task.is_dir() || task != current {
        return refuse("Task directory must be the current physical directory");
    }
    for name in [".command.sh", ".command.run"] {
        let value = match fs::symlink_metadata(task.join(name)) {
            Ok(value) => value,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return refuse(format!("Missing Nextflow task control: {name}"));
            }
            Err(error) => return Err(os(error)),
        };
        if !value.file_type().is_file() || value.nlink() != 1 {
            return refuse(format!(
                "Nextflow task control must be an independent regular file: {name}"
            ));
        }
    }
    let mut pending = vec![task.clone()];
    while let Some(directory) = pending.pop() {
        // os.walk ignores a directory it cannot list.
        let Ok(entries) = fs::read_dir(&directory) else {
            continue;
        };
        for entry in entries {
            let path = entry.map_err(os)?.path();
            let value = fs::symlink_metadata(&path).map_err(os)?;
            let kind = value.file_type();
            if kind.is_symlink() {
                continue;
            }
            if kind.is_dir() {
                pending.push(path);
                continue;
            }
            if !kind.is_file() || kind.is_fifo() || kind.is_socket() {
                return refuse(format!(
                    "Unsupported writable task entry: {}",
                    path.display()
                ));
            }
            if value.nlink() != 1 {
                return refuse(format!(
                    "Writable task contains a hard-linked file: {}",
                    path.display()
                ));
            }
        }
    }
    Ok(task)
}

/// Private scratch, without following an input's directory symlink.
fn private_directory(task: &Path, name: &str) -> Result<PathBuf, Refusal> {
    let path = task.join(name);
    let link = fs::symlink_metadata(&path).is_ok_and(|m| m.file_type().is_symlink());
    if link || (path.exists() && !path.is_dir()) {
        return refuse(format!(
            "Private task scratch is not a physical directory: {name}"
        ));
    }
    if !path.is_dir() {
        fs::create_dir(&path).map_err(os)?;
    }
    Ok(path)
}

fn command_argv(task: &Path, command: &[String]) -> Result<Vec<String>, Refusal> {
    if command.is_empty() {
        return refuse("Pass a command after --");
    }
    let temp = private_directory(task, ".priva_task_tmp")?;
    let cache = private_directory(task, ".priva_task_cache")?;
    let task = task.display().to_string();
    // Keep scientific inputs, HOME/configuration and package/model paths
    // intact. Only disposable scratch/cache defaults are redirected. No model
    // is fetched at runtime: network access is disabled inside this view.
    let mut argv: Vec<String> = [
        "--die-with-parent",
        "--new-session",
        "--unshare-pid",
        "--unshare-net",
        "--cap-drop",
        "ALL",
        "--ro-bind",
        "/",
        "/",
        "--bind",
        &task,
        &task,
        "--tmpfs",
        "/tmp",
        "--proc",
        "/proc",
        "--dev",
        "/dev",
        "--chdir",
        &task,
        "--setenv",
        "TMPDIR",
        &temp.display().to_string(),
        "--setenv",
        "XDG_CACHE_HOME",
        &cache.display().to_string(),
        "--setenv",
        "PYTHONDONTWRITEBYTECODE",
        "1",
        "--setenv",
        "PRIVA_TASK_BOUNDARY",
        "bwrap-readonly-v1",
        "--",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    argv.extend(command.iter().cloned());
    Ok(argv)
}

fn which(tool: &str) -> Option<PathBuf> {
    let search = std::env::var("PATH").unwrap_or_default();
    search.split(':').find_map(|directory| {
        let directory = if directory.is_empty() { "." } else { directory };
        let candidate = Path::new(directory).join(tool);
        let usable = fs::metadata(&candidate).is_ok_and(|m| !m.is_dir())
            && super::contract::access_executable(&candidate);
        usable.then_some(candidate)
    })
}

fn usage(message: &str) -> i32 {
    eprintln!(
        "usage: priva task-boundary {{check,exec}} --task-dir DIR [-- COMMAND ...]\npriva task-boundary: error: {message}"
    );
    2
}

/// The exit code: 0 for a passing check, 125 for a refusal, 2 for usage.
pub(crate) fn main(args: &[String]) -> i32 {
    let Some(procedure) = args.first().map(String::as_str) else {
        return usage("the following arguments are required: procedure");
    };
    if !matches!(procedure, "check" | "exec") {
        return usage(&format!(
            "argument procedure: invalid choice: '{procedure}' (choose from 'check', 'exec')"
        ));
    }
    let mut task_dir: Option<PathBuf> = None;
    let mut command: Vec<String> = Vec::new();
    let mut index = 1;
    while index < args.len() {
        let argument = &args[index];
        if argument == "--task-dir" {
            let Some(value) = args.get(index + 1) else {
                return usage("argument --task-dir: expected one argument");
            };
            task_dir = Some(PathBuf::from(value));
            index += 2;
        } else if let Some(value) = argument.strip_prefix("--task-dir=") {
            task_dir = Some(PathBuf::from(value));
            index += 1;
        } else if procedure == "exec" {
            command = args[index..].to_vec();
            break;
        } else {
            return usage(&format!("unrecognized arguments: {argument}"));
        }
    }
    let Some(task_dir) = task_dir else {
        return usage("the following arguments are required: --task-dir");
    };
    if command.first().map(String::as_str) == Some("--") {
        command.remove(0);
    }
    let outcome = (|| -> Result<(), Refusal> {
        let task = check_task(&task_dir)?;
        if procedure == "check" {
            return Ok(());
        }
        let Some(executable) = which("bwrap") else {
            return refuse("Bubblewrap (bwrap) is required for the task_readonly profile");
        };
        let argv = command_argv(&task, &command)?;
        let error = std::process::Command::new(&executable)
            .arg0(&executable)
            .args(&argv)
            .exec();
        Err(os(error))
    })();
    match outcome {
        Ok(()) => 0,
        Err(Refusal(message)) => {
            eprintln!("PriVA task boundary refused: {message}");
            125
        }
    }
}
