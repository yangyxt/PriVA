//! Observe the actual selected producer, VEP/Perl and scientific Python runtime.
//! This module never invokes annotation or numerical model generation.
use super::io::*;
use super::{error, object, require};
use crate::cli::Result;
use crate::installed_environment::yaml::{self, Value, Yaml};
use crate::json_text::Json;
use crate::python_path::{absolute, metadata, os_error, path_value};
use std::collections::BTreeSet;
use std::fs;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

pub(super) struct Settings {
    pub config: Json,
    pub source: PathBuf,
    pub paths: Json,
}
fn yaml_value(doc: &Yaml, node: usize) -> Result<Json> {
    Ok(match doc.value(node) {
        Value::Null => Json::Null,
        Value::Bool(v) => Json::Bool(*v),
        Value::Str(v) => string(v),
        Value::Number(v) if doc.tag(node) == "int" => Json::Int(v.parse().map_err(|_| {
            error(
                "UnsupportedConfiguration",
                "selected YAML integer exceeds the native metadata range",
            )
        })?),
        Value::Number(v) => Json::Float(v.parse().map_err(|_| {
            error(
                "UnsupportedConfiguration",
                "selected YAML float is not representable",
            )
        })?),
        // These fields are consumed only as scalars. Keep the container type
        // and emptiness for the original truth/path/int checks, without
        // recursively unfolding valid cyclic YAML aliases.
        Value::Sequence(v) => Json::Array(if v.is_empty() {
            vec![]
        } else {
            vec![Json::Null]
        }),
        Value::Mapping(v) => Json::Object(if v.is_empty() {
            vec![]
        } else {
            vec![("".into(), Json::Null)]
        }),
        Value::Other(_) | Value::Codepoints { .. } => {
            return Err(error(
                "UnsupportedConfiguration",
                format!(
                    "selected YAML {} is outside the scalar metadata profile",
                    doc.tag(node)
                ),
            ));
        }
    })
}
fn config(path: &Path) -> Result<Json> {
    let doc = yaml::load(&read_text(path)?)?;
    require(
        doc.mapping(doc.root).is_some(),
        "Installer configuration must be a mapping",
    )?;
    let mut values = Json::Object(Vec::new());
    // SafeLoader still constructs and validates the entire document. Unused
    // configuration values, including aliases/typed collections, are ignored
    // exactly as in the original helper; they are never emitted or rewritten.
    for key in [
        "cache_root",
        "alphamissense_vep_vcf",
        "alphamissense_vep_base_vcf",
        "alphamissense_vcf",
        "assembly",
        "vep_installed_version",
        "threads",
        "conservation_file",
        "vep_plugins_dir",
        "vep_cache_dir",
        "ref_genome",
        "alphamissense_prescore",
    ] {
        if let Some(node) = doc.get(doc.root, key) {
            values.set(key, yaml_value(&doc, node)?);
        }
    }
    Ok(values)
}
pub(super) fn settings(config_file: &Path, source: &Path, base: Option<&Path>) -> Result<Settings> {
    let config = config(config_file)?;
    require(
        config.as_object().is_some(),
        "Installer configuration must be a mapping",
    )?;
    let root = optional(&config, "cache_root")?;
    let root = if truth(root) {
        path(root)?
    } else {
        PathBuf::from(".")
    };
    require(
        root.is_absolute(),
        "Configure the managed cache_root before AlphaMissense annotation",
    )?;
    let canonical = path(get(&config, "alphamissense_vep_vcf")?)?;
    let selected = optional(&config, "alphamissense_vep_base_vcf")?;
    let base = if let Some(base) = base {
        base.to_owned()
    } else if truth(selected) {
        path(selected)?
    } else {
        let name = canonical.file_name().unwrap_or_default().to_string_lossy();
        require(
            !name.is_empty(),
            format!("{} has an empty name", canonical.display()),
        )?;
        canonical.with_file_name(format!(
            "{}.base.vcf.gz",
            name.strip_suffix(".vcf.gz").unwrap_or(&name)
        ))
    };
    let paths = object([
        (
            "converted",
            pstr(&path(get(&config, "alphamissense_vcf")?)?),
        ),
        ("base", pstr(&base)),
        ("canonical", pstr(&canonical)),
    ]);
    let mut selected = Vec::new();
    let mut present = Vec::new();
    for (_, value) in members(&paths)? {
        let p = path(value)?;
        require(
            p.is_absolute()
                && p.file_name()
                    .is_some_and(|n| n.to_string_lossy().ends_with(".vcf.gz"))
                && !p.components().any(|c| c == std::path::Component::ParentDir),
            "Use absolute managed AlphaMissense .vcf.gz destinations",
        )?;
        let resolved = resolve(&p)?;
        require(
            resolved.starts_with(resolve(&root)?) && !symlink(&p)?,
            format!(
                "AlphaMissense destinations must be independent files inside cache_root: {}",
                p.display()
            ),
        )?;
        selected.push(resolved);
    }
    require(
        selected.iter().collect::<BTreeSet<_>>().len() == selected.len(),
        "Converted, complete base and runtime VCF paths must differ",
    )?;
    for (_, v) in members(&paths)? {
        let p = path(v)?;
        if exists(&p)? {
            present.push(p);
        }
    }
    for i in 0..present.len() {
        for j in i + 1..present.len() {
            let a = metadata(&present[i])?;
            let b = metadata(&present[j])?;
            require(
                (a.dev(), a.ino()) != (b.dev(), b.ino()),
                "AlphaMissense base and runtime files must not be hardlinks",
            )?;
        }
    }
    Ok(Settings {
        config,
        source: source.to_owned(),
        paths,
    })
}
pub(super) fn assembly(value: &Json) -> Result<&'static str> {
    match value.as_str() {
        Some("hg19" | "GRCh37") => Ok("hg19"),
        Some("hg38" | "GRCh38") => Ok("hg38"),
        _ if matches!(value, Json::Array(_) | Json::Object(_)) => Err(error(
            "TypeError",
            format!("unhashable type: '{}'", type_name(value)),
        )),
        _ => Err(error(
            "ValueError",
            "AlphaMissense requires hg19/GRCh37 or hg38/GRCh38",
        )),
    }
}
pub(super) fn major(value: &str) -> bool {
    value
        .as_bytes()
        .first()
        .is_some_and(|b| matches!(b, b'1'..=b'9'))
        && value.bytes().all(|b| b.is_ascii_digit())
}
pub(super) fn procedure(path: &Path, name: &str) -> Result<String> {
    let value = read_text(path)?;
    let prefix = format!("function {name}() {{\n");
    let start = value
        .match_indices(&prefix)
        .find(|(i, _)| *i == 0 || value.as_bytes()[i - 1] == b'\n')
        .map(|v| v.0);
    if let Some(start) = start {
        if let Some(end) = value[start + prefix.len()..].find("\n}") {
            return Ok(value[start..start + prefix.len() + end + 2].to_owned());
        }
    }
    Err(error(
        "ValueError",
        format!("Cannot identify scientific shell procedure: {name}"),
    ))
}
fn annotation_recipe(path: &Path) -> Result<String> {
    let body = procedure(path, "AlphaMissense_anno")?;
    let lines: Vec<&str> = body.lines().collect();
    let converter = lines.iter().find(|line| {
        line.trim_start()
            .strip_prefix("python \"$convert_py\"")
            .is_some_and(|s| !s.is_empty())
    });
    let start = lines
        .iter()
        .position(|line| line.trim_start_matches([' ', '\t']) == "basic_vep_annotation \\");
    require(
        converter.is_some() && start.is_some(),
        "Cannot identify actual AlphaMissense converter/VEP invocation",
    )?;
    let mut selected = Vec::new();
    for line in &lines[start.unwrap()..] {
        selected.push(*line);
        let trimmed = line.trim_end();
        let before = trimmed.strip_suffix('\\').unwrap_or(trimmed).trim_end();
        if before.ends_with("&&")
            && before[..before.len() - 2]
                .chars()
                .last()
                .is_some_and(char::is_whitespace)
            || !trimmed.ends_with('\\')
        {
            break;
        }
    }
    Ok(crate::sha256::sha256_hex(
        format!(
            "{}\n{}",
            converter.unwrap().trim(),
            selected.join("\n").trim()
        )
        .as_bytes(),
    ))
}
fn execution(path: &Path) -> Result<String> {
    let body = procedure(path, "basic_vep_annotation")?;
    let mut start = None;
    let mut offset = 0;
    let mut whitespace_start = 0;
    for line in body.split_inclusive('\n') {
        let trim = line.trim_start();
        if trim.strip_prefix("vep").is_some_and(|v| {
            v.chars().next().is_some_and(char::is_whitespace) && !v.trim().is_empty()
        }) {
            start = Some(whitespace_start);
            break;
        }
        offset += line.len();
        if !line.trim().is_empty() {
            whitespace_start = offset;
        }
    }
    require(
        start.is_some(),
        "Cannot identify the actual VEP execution block",
    )?;
    let start = body
        .find("    if [[ -f ${spliceai_snv_prescore}")
        .unwrap_or(start.unwrap());
    let execution = body[start..]
        .lines()
        .filter(|l| !l.trim_start().starts_with("log "))
        .collect::<Vec<_>>()
        .join("\n");
    Ok(crate::sha256::sha256_hex(execution.as_bytes()))
}
pub(super) fn tree(directory: &Path, code: bool) -> Result<Json> {
    let root = resolve(directory)?;
    require(
        root.is_dir(),
        format!("Missing selected VEP dependency tree: {}", root.display()),
    )?;
    let mut rows = Vec::new();
    let mut total = 0i128;
    let mut latest_m = 0i128;
    let mut latest_c = 0i128;
    fn visit(
        root: &Path,
        dir: &Path,
        code: bool,
        rows: &mut Vec<Json>,
        total: &mut i128,
        latest_m: &mut i128,
        latest_c: &mut i128,
    ) -> Result<()> {
        // os.walk's default onerror ignores an unreadable directory.
        let Ok(entries) = fs::read_dir(dir) else {
            return Ok(());
        };
        let mut dirs = Vec::new();
        let mut files = Vec::new();
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_dir() {
                dirs.push(p)
            } else {
                files.push(p)
            }
        }
        dirs.sort();
        files.sort();
        for p in files {
            let m = metadata(&p)?;
            require(
                m.is_file(),
                format!("Nonregular selected VEP dependency: {}", p.display()),
            )?;
            let name = p.strip_prefix(root).unwrap().display().to_string();
            rows.push(if code {
                object([("name", string(name)), ("sha256", string(digest(&p)?))])
            } else {
                object([
                    ("name", string(name)),
                    ("size", Json::Int(m.len() as i128)),
                    ("inode", Json::Int(m.ino() as i128)),
                    ("mtime_ns", time(m.mtime(), m.mtime_nsec())),
                    ("ctime_ns", time(m.ctime(), m.ctime_nsec())),
                ])
            });
            *total += m.len() as i128;
            *latest_m = (*latest_m).max(m.mtime() as i128 * 1_000_000_000 + m.mtime_nsec() as i128);
            *latest_c = (*latest_c).max(m.ctime() as i128 * 1_000_000_000 + m.ctime_nsec() as i128);
        }
        for p in dirs {
            if !matches!(
                p.file_name().and_then(|n| n.to_str()),
                Some(".git" | "__pycache__")
            ) && !symlink(&p)?
            {
                visit(root, &p, code, rows, total, latest_m, latest_c)?;
            }
        }
        Ok(())
    }
    visit(
        &root,
        &root,
        code,
        &mut rows,
        &mut total,
        &mut latest_m,
        &mut latest_c,
    )?;
    require(
        !rows.is_empty(),
        format!("Empty selected VEP dependency tree: {}", root.display()),
    )?;
    Ok(object([
        ("path", pstr(&root)),
        ("files", Json::Int(rows.len() as i128)),
        ("bytes", Json::Int(total)),
        ("inventory_sha256", string(json_digest(&Json::Array(rows))?)),
        ("latest_mtime_ns", Json::Int(latest_m)),
        ("latest_ctime_ns", Json::Int(latest_c)),
    ]))
}
fn version(help: &str) -> Option<String> {
    for (start, _) in help.match_indices("ensembl-vep") {
        let raw = &help[start + 11..];
        let tail = if let Some(rest) = raw.trim_start().strip_prefix(':') {
            if !rest.chars().next().is_some_and(char::is_whitespace) {
                continue;
            }
            rest
        } else {
            if !raw.chars().next().is_some_and(char::is_whitespace) {
                continue;
            }
            raw
        };
        let tail = tail.trim_start();
        let mut parts = tail.split('.');
        let first = parts
            .next()?
            .chars()
            .take_while(char::is_ascii_digit)
            .collect::<String>();
        if first.is_empty() {
            continue;
        }
        let mut out = first;
        let mut tail = &tail[out.len()..];
        while let Some(next) = tail.strip_prefix('.') {
            let part = next
                .chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            if part.is_empty() {
                break;
            }
            out.push('.');
            out.push_str(&part);
            tail = &next[part.len()..];
        }
        return Some(out);
    }
    None
}
const PERL_QUERY: &str = "if (@ARGV) { require $ARGV[0]; } my %loaded = %INC; my @shared = @DynaLoader::dl_shared_objects; require JSON::PP; print JSON::PP::encode_json({version => $], shared => [@shared], loaded => {map {$_ => $loaded{$_}} grep {defined $loaded{$_} && !ref $loaded{$_}} keys %loaded}})";
pub(super) fn science(settings: &Settings) -> Result<Json> {
    let cfg = &settings.config;
    let source = &settings.source;
    let assembly = assembly(get(cfg, "assembly")?)?;
    let genome = if assembly == "hg19" {
        "GRCh37"
    } else {
        "GRCh38"
    };
    let major = text(get(cfg, "vep_installed_version")?);
    require(self::major(&major), "VEP major must be numeric")?;
    let executable = which("vep").ok_or_else(|| {
        error(
            "ValueError",
            "Activate the configured PriVA environment before recording annotation inputs",
        )
    })?;
    let executable = resolve(&executable)?;
    let help = run(&[executable.display().to_string(), "--help".into()], true)?;
    let version = version(&help);
    require(
        version
            .as_ref()
            .is_some_and(|v| v.split('.').next() == Some(major.as_str())),
        "Runnable VEP does not match the configured major",
    )?;
    let version = version.unwrap();
    let api = executable.parent().unwrap();
    let perl = which("perl")
        .ok_or_else(|| error("ValueError", "The active VEP Perl runtime is missing"))?;
    let conservation = optional(cfg, "conservation_file")?;
    let plugin = path(get(cfg, "vep_plugins_dir")?)?.join("Conservation.pm");
    let no_conservation = matches!(conservation, Json::Null)
        || matches!(conservation,Json::Str(s) if s.is_empty()||s=="null");
    let mut argv = vec![
        perl.display().to_string(),
        format!("-I{}", api.display()),
        format!("-I{}", api.join("modules").display()),
        "-MGetopt::Long".into(),
        "-MFindBin".into(),
        "-MBio::EnsEMBL::VEP::Runner".into(),
        "-MBio::EnsEMBL::VEP::Utils".into(),
        "-e".into(),
        PERL_QUERY.into(),
    ];
    if !no_conservation {
        argv.push(plugin.display().to_string());
    }
    let loaded = parse(&run(&argv, false)?)?;
    let files = get(&loaded, "loaded")?;
    require(
        files.get("Bio/EnsEMBL/VEP/Runner.pm").is_some(),
        "Cannot observe the API actually loaded by VEP",
    )?;
    let mut roots = BTreeSet::new();
    for (name, location) in members(files)? {
        if name.starts_with("Bio/EnsEMBL/") {
            let mut p = resolve(&path(location)?)?;
            for _ in path_value(name).components() {
                p.pop();
            }
            roots.insert(p.join("Bio"));
        }
    }
    let mut trees = Vec::new();
    for p in roots {
        trees.push(tree(&p, true)?);
    }
    let mut file_members = members(files)?.clone();
    file_members.sort_by(|a, b| a.0.cmp(&b.0));
    let mut loaded_files = Json::Object(Vec::new());
    for (name, location) in file_members {
        let p = path(&location)?;
        if p.is_file() {
            loaded_files.set(
                &name,
                object([
                    ("path", pstr(&resolve(&p)?)),
                    ("sha256", string(digest(&p)?)),
                ]),
            );
        }
    }
    if let Some(shared) = loaded.get("shared") {
        for location in list(shared)? {
            let p = path(location)?;
            let resolved = resolve(&p)?;
            loaded_files.set(
                &format!("XS:{}", resolved.display()),
                object([("path", pstr(&resolved)), ("sha256", string(digest(&p)?))]),
            );
        }
    }
    let cache = path(get(cfg, "vep_cache_dir")?)?
        .join("homo_sapiens_merged")
        .join(format!("{major}_{genome}"));
    let reference = path(get(cfg, "ref_genome")?)?;
    let threads = get(cfg, "threads")?;
    let threads = match threads {
        Json::Int(n) => *n,
        Json::Bool(v) => i128::from(*v),
        Json::Float(f) if f.is_nan() => {
            return Err(error("ValueError", "cannot convert float NaN to integer"));
        }
        Json::Float(f) if f.is_infinite() => {
            return Err(error(
                "OverflowError",
                "cannot convert float infinity to integer",
            ));
        }
        Json::Float(f) if *f >= i128::MIN as f64 && *f < -(i128::MIN as f64) => f.trunc() as i128,
        Json::Float(_) => {
            return Err(error(
                "UnsupportedConfiguration",
                "selected thread integer exceeds the native metadata range",
            ));
        }
        Json::Str(v) => {
            if !crate::python_text::python_int_is_valid(v) {
                return Err(error(
                    "ValueError",
                    format!(
                        "invalid literal for int() with base 10: {}",
                        crate::python_path::repr(v)
                    ),
                ));
            }
            let normalized: String = v
                .chars()
                .filter(|c| *c != '_')
                .map(|c| {
                    if let Some(digit) = crate::python_text::decimal_digit_value(c) {
                        char::from(b'0' + digit as u8)
                    } else if !c.is_ascii() && crate::python_text::is_python_space(c) {
                        ' '
                    } else {
                        c
                    }
                })
                .collect();
            normalized
                .trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r'))
                .parse()
                .map_err(|_| {
                    error(
                        "UnsupportedConfiguration",
                        "selected thread integer exceeds the native metadata range",
                    )
                })?
        }
        _ => {
            return Err(error(
                "TypeError",
                format!(
                    "int() argument must be a string, a bytes-like object or a real number, not '{}'",
                    type_name(threads)
                ),
            ));
        }
    };
    let converter = digest(&source.join("scripts/alphmis_tsv2vcf.py"))?;
    let execution = execution(&source.join("scripts/install_utils.sh"))?;
    let recipe = annotation_recipe(&source.join("scripts/install_utils.sh"))?;
    let prescore = path(get(cfg, "alphamissense_prescore")?)?;
    let prescore = Json::Array(vec![
        identity(&prescore, false)?,
        identity(&suffix(&prescore, ".tbi"), false)?,
    ]);
    let mut reference_files = vec![
        identity(&reference, false)?,
        identity(&suffix(&reference, ".fai"), false)?,
    ];
    if exists(&suffix(&reference, ".gzi"))? {
        reference_files.push(identity(&suffix(&reference, ".gzi"), false)?);
    }
    let conservation = if no_conservation {
        Json::Null
    } else {
        object([
            ("track", identity(&path(conservation)?, false)?),
            ("plugin_sha256", string(digest(&plugin)?)),
        ])
    };
    Ok(object([
        ("assembly", string(assembly)),
        ("vep_version", string(version)),
        ("threads", Json::Int(threads)),
        ("converter_sha256", string(converter)),
        ("basic_vep_execution_sha256", string(execution)),
        ("annotation_recipe_sha256", string(recipe)),
        ("prescore", prescore),
        ("reference", Json::Array(reference_files)),
        ("conservation", conservation),
        (
            "vep",
            object([
                ("executable", pstr(&executable)),
                ("sha256", string(digest(&executable)?)),
                ("api", Json::Array(trees)),
                ("perl", identity(&perl, false)?),
                ("perl_version", get(&loaded, "version")?.clone()),
                (
                    "loaded",
                    object([
                        ("files", Json::Int(members(&loaded_files)?.len() as i128)),
                        ("inventory_sha256", string(json_digest(&loaded_files)?)),
                        (
                            "runner_path",
                            pstr(&resolve(&path(get(files, "Bio/EnsEMBL/VEP/Runner.pm")?)?)?),
                        ),
                    ]),
                ),
            ]),
        ),
        ("cache", tree(&cache, false)?),
    ]))
}

// Python is observed, not used as a fallback for the native helper. Its import
// resolution and sys.version are properties of that actual scientific runtime;
// guessing them from a conda directory would change the original contract.
const PYTHON_OBSERVER: &str = r#"import importlib,importlib.metadata,json,sys
packages=[]
for name,dist in [('numpy','numpy'),('pandas','pandas'),('scipy','scipy'),('sklearn','scikit-learn'),('pysam','pysam'),('joblib','joblib')]:
 m=importlib.import_module(name); d=importlib.metadata.distribution(dist)
 packages.append(dict(name=name,distribution=dist,version=str(m.__version__),distribution_version=d.version,module_file=m.__file__,owned=[dict(name=str(f),path=str(d.locate_file(f))) for f in d.files or () if '__pycache__' not in f.parts and not str(f).endswith(('.pyc','.pyo'))]))
print(json.dumps(dict(python=sys.executable,version=sys.version,packages=packages)))
"#;
pub(super) fn motif_runtime(python: &str) -> Result<Json> {
    let observed = parse(&run(
        &[python.into(), "-c".into(), PYTHON_OBSERVER.into()],
        false,
    )?)?;
    let mut packages = Json::Object(Vec::new());
    for package in list(get(&observed, "packages")?)? {
        let mut owned = Vec::new();
        let mut latest = i128::MIN;
        for item in list(get(package, "owned")?)? {
            let p = PathBuf::from(absolute(&path(get(item, "path")?)?)?);
            let m = match fs::metadata(&p) {
                Ok(m) => m,
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => continue,
                Err(e) => return Err(os_error(e, &p, None)),
            };
            if m.is_file() {
                latest = latest.max(m.ctime() as i128 * 1_000_000_000 + m.ctime_nsec() as i128);
                owned.push(object([
                    ("name", get(item, "name")?.clone()),
                    ("path", pstr(&p)),
                    ("inode", Json::Int(m.ino() as i128)),
                    ("size", Json::Int(m.len() as i128)),
                    ("mtime_ns", time(m.mtime(), m.mtime_nsec())),
                    ("ctime_ns", time(m.ctime(), m.ctime_nsec())),
                ]));
            }
        }
        require(
            !owned.is_empty(),
            format!(
                "Cannot identify scientific Python package files: {}",
                text(get(package, "distribution")?)
            ),
        )?;
        owned.sort_by(|a, b| {
            a.get("name")
                .and_then(Json::as_str)
                .cmp(&b.get("name").and_then(Json::as_str))
        });
        packages.set(
            &text(get(package, "name")?),
            object([
                ("version", get(package, "version")?.clone()),
                (
                    "distribution_version",
                    get(package, "distribution_version")?.clone(),
                ),
                (
                    "module_path",
                    pstr(&resolve(&path(get(package, "module_file")?)?)?),
                ),
                ("files", Json::Int(owned.len() as i128)),
                (
                    "inventory_sha256",
                    string(json_digest(&Json::Array(owned))?),
                ),
                ("latest_ctime_ns", Json::Int(latest)),
            ]),
        );
    }
    Ok(object([
        (
            "python",
            identity(&path(get(&observed, "python")?)?, false)?,
        ),
        ("version", get(&observed, "version")?.clone()),
        ("packages", packages),
    ]))
}
