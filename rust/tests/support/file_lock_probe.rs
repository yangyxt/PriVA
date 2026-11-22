//! Process-boundary probe for the shared installer lock, not a production CLI.
#[path = "../../src/file_lock.rs"]
mod file_lock;
use file_lock::{ExclusiveLock, Wait};
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};

fn run() -> io::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() != 6 || args[0] != "--path" || args[2] != "--open" || args[4] != "--wait" {
        return Err(io::Error::other(
            "usage: file_lock_probe --path PATH --open append|directory --wait blocking|nonblocking",
        ));
    }
    let file = match args[3].as_str() {
        "append" => OpenOptions::new()
            .append(true)
            .create(true)
            .open(&args[1])?,
        "directory" => File::open(&args[1])?,
        _ => return Err(io::Error::other("invalid --open")),
    };
    let wait = match args[5].as_str() {
        "blocking" => Wait::Blocking,
        "nonblocking" => Wait::Nonblocking,
        _ => return Err(io::Error::other("invalid --wait")),
    };
    println!("opened");
    io::stdout().flush()?;
    let lock = ExclusiveLock::acquire(file, wait)?;
    println!("acquired");
    io::stdout().flush()?;
    let mut line = String::new();
    io::stdin().lock().read_line(&mut line)?;
    match line.trim_end() {
        "release" => lock.release()?,
        "drop" | "" => drop(lock),
        "exit" => std::process::exit(0),
        _ => return Err(io::Error::other("invalid lifetime command")),
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("errno={:?}; kind={:?}", error.raw_os_error(), error.kind());
        std::process::exit(1);
    }
}
