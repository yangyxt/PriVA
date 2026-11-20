use super::*;
use files::*;
use plan::Adoption;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::os::unix::fs::{OpenOptionsExt, PermissionsExt};

fn require_managed_parent(destination: &Path, root: &Path) -> Result<()> {
    let mut ancestor = destination.parent().unwrap_or(Path::new("."));
    while !exists(ancestor)? && !is_link(ancestor)? {
        ancestor = ancestor.parent().unwrap_or(Path::new("/"));
    }
    require(
        resolve(ancestor)?.starts_with(resolve(root)?),
        format!(
            "Refusing to copy through a managed-directory link into an external source: {}",
            destination.display()
        ),
    )
}
fn copy_file(
    source: &Path,
    destination: &Path,
    root: &Path,
    records: &mut Vec<Json>,
) -> Result<()> {
    let before = cache_identity(source)?;
    require_managed_parent(destination, root)?;
    mkdir(destination.parent().unwrap_or(Path::new(".")))?;
    let (file, temporary) =
        temporary(destination.parent().unwrap_or(Path::new(".")), ".adopting-")?;
    drop(file);
    let result = (|| {
        copy2(source, &temporary)?;
        let sha = sha256_file(source)?;
        require(
            json_text::py_equal(&cache_identity(source)?, &before)
                && sha256_file(&temporary)? == sha,
            format!(
                "Source changed or copy differs during cache adoption: {}",
                source.display()
            ),
        )?;
        if exists(destination)?
            && !is_link(destination)?
            && is_file(destination)?
            && sha256_file(destination)? == sha
        {
            unlink_missing(&temporary)?;
        } else {
            if exists(destination)? || is_link(destination)? {
                let relative = destination
                    .strip_prefix(root)
                    .map_err(|_| error("ValueError", "destination is not in cache_root"))?;
                let retained = root
                    .join(".installer-retained")
                    .join(stamp()?)
                    .join(relative);
                mkdir(retained.parent().unwrap())?;
                rename(destination, &retained)?;
            }
            rename(&temporary, destination)?;
        }
        let bytes = if let Json::Array(v) = &before {
            v[1].clone()
        } else {
            unreachable!()
        };
        records.push(Json::Object(vec![
            ("source".into(), Json::Str(source.display().to_string())),
            (
                "destination".into(),
                Json::Str(destination.display().to_string()),
            ),
            ("bytes".into(), bytes),
            ("sha256".into(), Json::Str(sha)),
            ("source_identity".into(), before),
        ]));
        Ok(())
    })();
    unlink_missing(&temporary)?;
    result
}
fn destination_links(path: &Path) -> Result<()> {
    // os.walk does not traverse directory links. Its scandir order is preserved.
    let children = match fs::read_dir(native(path)?) {
        Ok(v) => v,
        Err(_) => return Ok(()),
    };
    let mut dirs = vec![];
    for child in children.flatten() {
        if child.path().is_dir() {
            dirs.push(crate::python_unicode::from_native(&child.path()));
        }
    }
    for child in &dirs {
        require(
            !is_link(child)?,
            format!(
                "Refusing mutable cache adoption through an existing directory link: {}",
                child.display()
            ),
        )?;
    }
    for child in dirs {
        destination_links(&child)?;
    }
    Ok(())
}
fn caught_os(cause: &dyn std::error::Error) -> Option<String> {
    let text = cause.to_string();
    [
        "OSError",
        "FileNotFoundError",
        "FileExistsError",
        "PermissionError",
        "IsADirectoryError",
        "NotADirectoryError",
        "SpecialFileError",
    ]
    .into_iter()
    .find_map(|k| text.strip_prefix(&format!("{k}: ")).map(str::to_owned))
}
fn copytree(
    source: &Path,
    dest: &Path,
    root: &Path,
    prepared: &[PathBuf],
    records: &mut Vec<Json>,
    depth: usize,
) -> Result<Vec<(String, String, String)>> {
    if depth > 990 {
        return Err(error("RecursionError", "maximum recursion depth exceeded"));
    }
    let entries = fs::read_dir(native(source)?)
        .map_err(|e| os_error(e, source, None))?
        .collect::<std::io::Result<Vec<_>>>()
        .map_err(|e| os_error(e, source, None))?;
    let mut ignored = vec![];
    for entry in &entries {
        if prepared.contains(&resolve(&crate::python_unicode::from_native(
            &entry.path(),
        ))?) {
            ignored.push(entry.file_name());
        }
    }
    mkdir(dest)?;
    let mut errors = vec![];
    for entry in entries {
        if ignored.contains(&entry.file_name()) {
            continue;
        }
        let src = crate::python_unicode::from_native(&entry.path());
        let dst = dest.join(crate::python_unicode::from_native(Path::new(
            &entry.file_name(),
        )));
        let result = if is_dir(&src)? {
            copytree(&src, &dst, root, prepared, records, depth + 1)
                .map(|nested| errors.extend(nested))
        } else {
            copy_file(&src, &dst, root, records)
        };
        if let Err(e) = result {
            if let Some(message) = caught_os(e.as_ref()) {
                errors.push((
                    src.display().to_string(),
                    dst.display().to_string(),
                    message,
                ));
            } else {
                return Err(e);
            }
        }
    }
    if let Err(e) = copystat(source, dest) {
        if let Some(message) = caught_os(e.as_ref()) {
            errors.push((
                source.display().to_string(),
                dest.display().to_string(),
                message,
            ));
        } else {
            return Err(e);
        }
    }
    Ok(errors)
}
pub(crate) fn adopt_sources(ctx: &Context, plan: &[Adoption], root: &Path) -> Result<()> {
    if plan.is_empty() {
        return Ok(());
    }
    mkdir(root)?;
    let mut records = vec![];
    let prepared = plan
        .iter()
        .filter(|i| i.kind == "prepared_projection")
        .map(|i| resolve(&i.source))
        .collect::<Result<Vec<_>>>()?;
    let lock_path = root.join(".installer-adoption.lock");
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .mode(0o666)
        .open(native(&lock_path)?)
        .map_err(|e| os_error(e, &lock_path, None))?;
    let _lock = crate::file_lock::ExclusiveLock::acquire(file, crate::file_lock::Wait::Blocking)?;
    for item in plan {
        let (source, destination) = (&item.source, &item.destination);
        require(
            destination.starts_with(root),
            format!(
                "Adoption destination escapes cache_root: {}",
                destination.display()
            ),
        )?;
        if item.kind == "prepared_projection" {
            require_managed_parent(destination, root)?;
            plan::adopt_prepared(ctx, source, destination, root)?;
            records.push(Json::Object(vec![
                ("source".into(), Json::Str(source.display().to_string())),
                (
                    "destination".into(),
                    Json::Str(destination.display().to_string()),
                ),
                ("kind".into(), Json::Str("prepared_projection".into())),
                ("source_identity".into(), cache_identity(source)?),
            ]));
        } else if item.kind == "directory" {
            require(
                !resolve(destination)?.starts_with(resolve(source)?),
                format!(
                    "Cannot adopt a cache directory into itself: {}",
                    source.display()
                ),
            )?;
            require_managed_parent(destination, root)?;
            require(
                !is_link(destination)?,
                format!(
                    "Refusing mutable directory adoption through an existing link: {}",
                    destination.display()
                ),
            )?;
            if is_dir(destination)? {
                destination_links(destination)?;
            }
            let errors = copytree(source, destination, root, &prepared, &mut records, 0)?;
            if !errors.is_empty() {
                return Err(error(
                    "Error",
                    format!(
                        "[{}]",
                        errors
                            .iter()
                            .map(|(a, b, c)| format!("({}, {}, {})", repr(a), repr(b), repr(c)))
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                ));
            }
        } else {
            copy_file(source, destination, root, &mut records)?;
        }
    }
    for item in &records {
        let source = Path::new(item.get("source").unwrap().as_str().unwrap());
        require(
            json_text::py_equal(
                &cache_identity(source)?,
                item.get("source_identity").unwrap(),
            ),
            "An explicit cache source changed during adoption; original config is retained",
        )?;
    }
    let receipt = root.join("installer-source-adoptions.json");
    let mut previous = if is_file(&receipt)? {
        receipt_json::RecordJson::parse(&read_text(&receipt)?)?
    } else {
        receipt_json::RecordJson::from_json(&Json::Object(vec![("schema_version".into(),Json::Int(1)),("copies".into(),Json::Array(vec![])),("scientific_acceptance".into(),Json::Bool(false)),("scope".into(),Json::Str("Byte-preserving source adoption; resource-specific scientific reuse checks remain required.".into()))]))
    };
    previous.get_mut("copies")?.extend_unique(&records)?;
    let (mut stream, temporary) = temporary(root, ".adoptions-")?;
    // The source leaves a receipt temporary on dump/write/replace failure.
    let data = previous.dumps(true)? + "\n";
    stream
        .write_all(data.as_bytes())
        .map_err(|e| os_error(e, &temporary, None))?;
    drop(stream);
    rename(&temporary, &receipt)
}
pub(crate) fn render_config(original: &str, updates: &Paths) -> Result<String> {
    if original
        .trim_start_matches(crate::python_text::is_python_space)
        .starts_with('{')
    {
        let mut value = receipt_json::RecordJson::parse(original)?;
        value.update(updates)?;
        return Ok(value.dumps(false)? + "\n");
    }
    crate::config_text::render_yaml_scalars(
        original,
        &updates
            .iter()
            .map(|(k, v)| {
                Ok((
                    k.clone(),
                    receipt_json::RecordJson::Str(crate::python_unicode::points(v)).dumps(true)?,
                ))
            })
            .collect::<Result<Vec<_>>>()?,
        crate::config_text::Writer::Paths,
    )
    .map_err(|e| error(e.kind, e.message))
}
enum ConfigText {
    Unicode(String),
    EncodingFailure(String),
}
impl ConfigText {
    fn render(original: &str, updates: &Paths) -> Result<Self> {
        match render_config(original, updates) {
            Ok(value) => Ok(Self::Unicode(value)),
            Err(cause) => match cause.to_string().strip_prefix("UnicodeEncodeError: ") {
                Some(message) => Ok(Self::EncodingFailure(message.into())),
                None => Err(cause),
            },
        }
    }
    fn equals(&self, text: &str) -> bool {
        matches!(self, Self::Unicode(value) if value == text)
    }
    fn encode(&self) -> Result<&str> {
        match self {
            Self::Unicode(value) => Ok(value),
            Self::EncodingFailure(message) => Err(error("UnicodeEncodeError", message)),
        }
    }
}
fn publish_config(
    config: &Path,
    destination: &Path,
    original: &str,
    updated: &ConfigText,
    message: &str,
) -> Result<()> {
    let prefix = format!(".{}.", name(destination));
    let (mut file, temporary) = temporary(
        destination
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new(".")),
        &prefix,
    )?;
    let result = (|| {
        // Python renders Unicode first, creates the temporary, and only then
        // encodes it in TextIOWrapper.write. Preserve intervening I/O failures
        // and final cleanup when an ignored JSON scalar contains a surrogate.
        file.write_all(updated.encode()?.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        drop(file);
        use std::os::unix::fs::MetadataExt;
        fs::set_permissions(
            native(&temporary)?,
            fs::Permissions::from_mode(metadata(config)?.mode() & 0o777),
        )
        .map_err(|e| os_error(e, &temporary, None))?;
        require(read_text(config)? == original, message)?;
        rename(&temporary, destination)
    })();
    if exists(&temporary)? {
        unlink_missing(&temporary)?;
    }
    result
}
pub(crate) fn resolve_config(
    ctx: &Context,
    path: &Path,
    source_root: &Path,
    check: bool,
    managed: Option<&str>,
    output: Option<&Path>,
) -> Result<Option<Paths>> {
    let original = read_text(path)?;
    let config = Config::parse(&original)?;
    let managed = if let Some(s) = managed.filter(|s| !s.is_empty()) {
        let value = json_text::parse(s, false).map_err(|e| error("JSONDecodeError", e))?;
        let pairs = value.as_object().ok_or_else(|| {
            error(
                "ValueError",
                "Managed path updates must be a JSON mapping of path strings",
            )
        })?;
        pairs
            .iter()
            .map(|(k, v)| {
                v.as_str()
                    .map(|v| {
                        (
                            crate::python_unicode::encode(k),
                            crate::python_unicode::encode(v),
                        )
                    })
                    .ok_or_else(|| {
                        error(
                            "ValueError",
                            "Managed path updates must be a JSON mapping of path strings",
                        )
                    })
            })
            .collect::<Result<Paths>>()?
    } else {
        vec![]
    };
    let source_root = PathBuf::from(absolute(source_root)?);
    let mut updates = defaults::portable_updates(ctx, &config, &source_root)?;
    require(
        !managed
            .iter()
            .any(|(k, v)| get(&updates, k).is_ok_and(|u| u != v)),
        "Managed updates differ from the original configuration/cache_root plan",
    )?;
    let plan = plan::adoption_plan(ctx, &config, &updates, &source_root)?;
    if check {
        return Ok(Some(updates));
    }
    if output.is_some() {
        let dirs = FOLDERS
            .iter()
            .map(|(_, _, k)| Ok((k.to_string(), get(&updates, k)?.to_owned())))
            .collect::<Result<Paths>>()?;
        for (key, _) in plan::resource_groups(ctx, &config, &dirs)? {
            let value = config.get(&key);
            let value = if key == "hg38_hg19_chain" && !value.empty() {
                Some(
                    plan::reverse_chain_source(ctx, value, Some(&source_root))?
                        .display()
                        .to_string(),
                )
            } else {
                value.string().map(str::to_owned)
            };
            if let Some(value) =
                value.filter(|v| !v.is_empty() && !FOLDERS.iter().any(|(_, _, k)| *k == key))
            {
                let source = tempfile_parent(Path::new(&value))?;
                if is_file(&source)?
                    || is_dir(&source)?
                    || (!GENERATED.contains(&key.as_str())
                        && plan::bundled_resource(ctx, &key, &value, Some(&source_root))?.is_none())
                {
                    set(&mut updates, &key, value);
                }
            }
        }
    } else {
        adopt_sources(
            ctx,
            &plan,
            &plan::checked_storage(config.index("cache_root")?, "cache_root")?,
        )?;
    }
    let updated = ConfigText::render(&original, &updates)?;
    let destination = output.unwrap_or(path);
    if destination == path && updated.equals(&original) {
        return Ok(None);
    }
    publish_config(
        path,
        destination,
        &original,
        &updated,
        "Configuration changed during path resolution; retry after its writer closes",
    )?;
    Ok(None)
}
pub(crate) fn configure_cache_paths(
    ctx: &Context,
    path: &Path,
    config: &Config,
    mode: &str,
    source_root: Option<&Path>,
) -> Result<Option<Paths>> {
    require(
        matches!(mode, "write" | "--check"),
        "Use priva_configure_cache_paths CONFIG [--check]",
    )?;
    let root = match source_root {
        Some(p) => p.to_path_buf(),
        None => ctx.source_root()?,
    };
    let (directories, updates) = plan::resolve_cache_paths(ctx, config, Some(&root))?;
    if mode == "--check" {
        return Ok(Some(updates));
    }
    if updates.is_empty() {
        return Ok(None);
    }
    let original = read_text(path)?;
    let updated = ConfigText::render(&original, &updates)?;
    adopt_sources(
        ctx,
        &plan::adoption_plan(ctx, config, &updates, &root)?,
        &plan::checked_storage(config.index("cache_root")?, "cache_root")?,
    )?;
    for (_, dir) in &directories {
        mkdir(Path::new(dir))?;
    }
    if !updated.equals(&original) {
        publish_config(
            path,
            path,
            &original,
            &updated,
            "Configuration changed during source adoption; original selections retained",
        )?;
    }
    for (k, d) in directories {
        println!(
            "{}: {}",
            crate::python_unicode::display(&k),
            crate::python_unicode::display(&d)
        );
    }
    Ok(None)
}
