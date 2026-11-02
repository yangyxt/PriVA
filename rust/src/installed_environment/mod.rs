//! Optional package-only export with the original source's publication order.
//! The manager is read-only; no solve, installation or activation occurs here.
mod annex;
mod ffi;
mod files;
mod process;
mod yaml;
use crate::cli::{Result, invalid};
use crate::python_path::{is_file, os_error, path_value, repr};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use yaml::Value;
const HELP: &str = "usage: priva installed-environment export --config CONFIG --manager MANAGER [--name NAME] [--prefix PREFIX]\n       priva installed-environment export-installation --config CONFIG --manager MANAGER\n\nExport installed packages only. No solve, package installation, or environment activation.";
fn error(kind: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(kind, invalid(message))
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn text(bytes: &[u8]) -> Result<String> {
    crate::python_io::decode_utf8(bytes)
        .map(|v| v.replace("\r\n", "\n").replace('\r', "\n"))
        .map_err(|m| error("UnicodeDecodeError", m))
}
fn package(mut document: yaml::Yaml, name: &str, prefix: Option<&Path>) -> Result<yaml::Yaml> {
    let root = document.root;
    require(
        document.mapping(root).is_some(),
        "The package manager exported a different environment",
    )?;
    let reported_name = document.get(root, "name");
    if let Some(prefix) = prefix {
        let reported = document
            .get(root, "prefix")
            .and_then(|n| document.string(n));
        let same = if let Some(reported) = reported {
            files::realpath(&path_value(reported), true)? == files::realpath(prefix, true)?
        } else {
            false
        };
        require(
            same,
            "The package manager exported a different environment prefix",
        )?;
    }
    require(
        reported_name.is_none_or(|n| {
            document.string(n) == Some(name)
                || (matches!(document.value(n), Value::Null) && prefix.is_some())
        }),
        "The package manager exported a different environment",
    )?;
    let deps = document.get(root, "dependencies");
    require(
        deps.and_then(|n| document.sequence(n))
            .is_some_and(|v| !v.is_empty()),
        "The installed environment export has no dependencies",
    )?;
    let deps = deps.unwrap();
    let string_ok = |n| {
        document
            .string(n)
            .is_some_and(|s| !crate::python_text::python_strip(s).is_empty())
    };
    for entry in document.sequence(deps).unwrap() {
        let pip_ok = document.mapping(*entry).is_some_and(|map| {
            map.len() == 1
                && document.string(map[0].0) == Some("pip")
                && document
                    .sequence(map[0].1)
                    .is_some_and(|v| v.iter().all(|n| string_ok(*n)))
        });
        require(
            string_ok(*entry) || pip_ok,
            "The installed environment export has an unsupported dependency",
        )?;
    }
    let channels = document.get(root, "channels");
    require(
        channels.is_none_or(|n| {
            document
                .sequence(n)
                .is_some_and(|v| v.iter().all(|n| document.string(*n).is_some()))
        }),
        "The installed environment export has invalid channels",
    )?;
    if let Some(channels) = channels {
        for node in document.sequence(channels).unwrap() {
            let channel = document.string(*node).unwrap();
            let url = crate::python_url::split(channel)?;
            let private = url.netloc.rsplit_once('@').is_some_and(|(auth, _)| {
                auth.split_once(':')
                    .map_or(!auth.is_empty(), |(user, password)| {
                        !user.is_empty() || !password.is_empty()
                    })
            });
            require(
                !private
                    && url.query.is_empty()
                    && url.fragment.is_empty()
                    && url.scheme != "file"
                    && !channel.starts_with('/'),
                "A channel cannot be exported as a portable public package source",
            )?;
        }
    }
    document.package_root(name, channels, deps);
    Ok(document)
}
fn previous_matches(document: &yaml::Yaml, name: &str) -> bool {
    let root = document.root;
    document.mapping(root).is_some_and(|map| {
        map.len() == 3
            && map.iter().all(|(k, _)| {
                document
                    .string(*k)
                    .is_some_and(|s| ["name", "channels", "dependencies"].contains(&s))
            })
    }) && document.get(root, "name").and_then(|n| document.string(n)) == Some(name)
}
fn export(config_file: &Path, manager: &str, name: &str, prefix: Option<&Path>) -> Result<PathBuf> {
    require(
        name.bytes()
            .next()
            .is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_')
            && name
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-')),
        "The environment name must be a simple name",
    )?;
    let config = yaml::load(&files::read_text(config_file)?)?;
    let selected = config
        .get(config.root, "cache_root")
        .and_then(|n| config.string(n));
    require(
        selected.is_some_and(|s| Path::new(s).is_absolute()),
        "Set cache_root to the existing absolute cache directory",
    )?;
    let root = path_value(selected.unwrap());
    require(
        files::is_dir(&root)? && files::realpath(&root, true)? != Path::new("/"),
        "cache_root must be an existing dedicated directory",
    )?;
    let output = root.join(format!("{name}.installed.yml"));
    let before = if files::is_link(&output)? {
        Some(annex::snapshot(&output)?)
    } else {
        None
    };
    let selection = if let Some(prefix) = prefix {
        require(
            prefix.is_absolute() && is_file(&prefix.join("conda-meta/history"))?,
            "Select a completed absolute Mamba environment prefix",
        )?;
        vec!["--prefix".into(), prefix.display().to_string()]
    } else {
        vec!["--name".into(), name.into()]
    };
    let reported = process::manager(manager, &selection)?;
    let document = package(yaml::load(&reported)?, name, prefix)?;
    let payload = document.dump()?;
    if let Some(before) = before {
        let after = annex::snapshot(&output)?;
        let previous = yaml::load_bytes(&after.payload)?;
        require(
            before == after
                && previous_matches(&previous, name)
                && after.payload == payload.as_bytes(),
            "Linked installed-package export changed or differs from the installed packages; retained",
        )?;
        return Ok(output);
    }
    require(
        !files::is_link(&output)?,
        "Refusing a linked installed-package export destination",
    )?;
    if files::exists(&output)? {
        let previous = yaml::load(&files::read_text(&output)?)?;
        require(
            previous_matches(&previous, name),
            "Existing file is not a matching installed-package export; retained",
        )?;
    }
    if files::exists(&output)? && files::read_text(&output)? == payload {
        return Ok(output);
    }
    files::publish(&output, &payload)?;
    Ok(output)
}
fn installation(config_file: &Path, manager: &str) -> Result<Vec<PathBuf>> {
    let config = yaml::load(&files::read_text(config_file)?)?;
    let mut exports = vec![export(config_file, manager, "priva_acmg", None)?];
    let base = config
        .get(config.root, "cadd_base_dir")
        .and_then(|n| config.string(n));
    if let Some(base) = base.filter(|s| !s.is_empty() && *s != "null") {
        let recipes = path_value(base).join("envs");
        let root = files::realpath(&recipes.join("conda"), false)?;
        let entries = match fs::read_dir(&recipes) {
            Ok(entries) => Some(entries),
            Err(e)
                if matches!(
                    e.raw_os_error(),
                    Some(libc::ENOENT | libc::ENOTDIR | libc::EACCES)
                ) =>
            {
                None
            }
            Err(e) => return Err(os_error(e, &recipes, None)),
        };
        let mut paths = vec![];
        if let Some(entries) = entries {
            for entry in entries {
                let entry = entry?;
                if entry.file_name().to_string_lossy().ends_with(".yml") {
                    paths.push(entry.path());
                }
            }
        }
        paths.sort();
        for recipe in paths {
            use std::os::unix::ffi::OsStrExt;
            let mut bytes = root.as_os_str().as_bytes().to_vec();
            bytes.extend(fs::read(&recipe).map_err(|e| os_error(e, &recipe, None))?);
            let prefix = recipes
                .join("conda")
                .join(format!("{}_", crate::md5::md5_hex(&bytes)));
            if is_file(&prefix.join("conda-meta/history"))?
                && files::exists(&PathBuf::from(format!(
                    "{}.env_setup_done",
                    prefix.display()
                )))?
            {
                let name = format!("cadd-{}", recipe.file_stem().unwrap().to_string_lossy());
                exports.push(export(config_file, manager, &name, Some(&prefix))?);
            }
        }
    }
    Ok(exports)
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.first().is_some_and(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let command = args
        .first()
        .ok_or_else(|| invalid("the following arguments are required: procedure"))?;
    let names: &[&str] = match command.as_str() {
        "export" => &["config", "manager", "name", "prefix"],
        "export-installation" => &["config", "manager"],
        _ => {
            return Err(invalid(format!(
                "argument procedure: invalid choice: {} (choose from 'export', 'export-installation')",
                repr(command)
            )));
        }
    };
    let Some(options) = crate::python_args::parse_long_option_occurrences(&args[1..], names, &[])
        .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let options: HashMap<_, _> = options.into_iter().collect();
    crate::python_args::require_options(&options, &["config", "manager"]).map_err(invalid)?;
    let config = path_value(&options["config"]);
    let manager = &options["manager"];
    crate::builder_error::finish((|| {
        let paths = if command == "export-installation" {
            installation(&config, manager)?
        } else {
            let prefix = options.get("prefix").map(|s| path_value(s));
            vec![export(
                &config,
                manager,
                options.get("name").map_or("priva_acmg", String::as_str),
                prefix.as_deref(),
            )?]
        };
        for path in paths {
            println!("{}", path.display());
        }
        Ok(())
    })())
}
