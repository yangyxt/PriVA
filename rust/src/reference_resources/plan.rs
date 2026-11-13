//! `checked_path` and `resource_plan`: where each reference resource lives,
//! and whether it is a managed default (pinned, prepared when missing) or an
//! explicitly supplied original that is only checked.
use super::{Result, require, value_error};
use crate::installer_prerequisites::checks::assembly_alias;
use crate::installer_prerequisites::paths::{
    abspath, exists, is_relative_to, is_symlink, join, lexists, parent, parts,
};
use crate::installer_prerequisites::{Config, Value};
use crate::python_path::{metadata, os_error};
use crate::python_text::is_python_space;
use std::collections::VecDeque;
use std::fs;
use std::os::unix::fs::MetadataExt;
use std::path::Path;

/// `CHR1_LENGTH` keys, in the script's order.
pub(super) const ASSEMBLIES: [&str; 2] = ["hg19", "hg38"];
const CHAIN_KEY: &str = "gofcards_hg19_to_hg38_chain";
const CHAIN_FILE: &str = "liftOver/hg19ToHg38.over.chain";
const PATH_TEXT: &str =
    "Reference resource paths must be nonempty path strings without tabs/newlines";
const DISTINCT: &str = "GoFCards hg19 and hg38 FASTAs must be distinct assembly-correct inputs";
/// Index, source, staging and lock companions that keep a legacy default in place.
const LEGACY_STATE: [&str; 11] = [
    ".fai",
    ".fai.partial",
    ".gzi",
    ".gz",
    ".gz.partial",
    ".partial",
    ".partial.fai",
    ".partial.fai.partial",
    ".source.json",
    ".source.json.partial",
    ".install.lock",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    Fasta,
    Chain,
}

impl Kind {
    pub(super) fn name(self) -> &'static str {
        match self {
            Kind::Fasta => "fasta",
            Kind::Chain => "chain",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Scope {
    Selected,
    Both,
    Gofcards,
}

/// One planned resource: its path, its assembly (or `hg19_to_hg38` for the
/// chain), the configuration key that selected it, and whether it is managed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Resource {
    pub(super) path: String,
    pub(super) assembly: &'static str,
    pub(super) key: String,
    pub(super) managed: bool,
    pub(super) kind: Kind,
}

/// The absolute roots named by PRIVA_INSTALL_INSPECT_ROOTS.
fn inspection_roots() -> Result<Vec<String>> {
    let setting = std::env::var_os("PRIVA_INSTALL_INSPECT_ROOTS").unwrap_or_default();
    setting
        .to_string_lossy()
        .split(':')
        .filter(|root| !root.is_empty())
        .map(abspath)
        .collect()
}

/// `checked_path(value)`: the absolute path of a path text. When
/// PRIVA_INSTALL_INSPECT_ROOTS names roots, the path must be inside one, and
/// each symbolic link on the way is checked before it is followed, so that no
/// target outside the roots is inspected.
pub(super) fn checked_path(value: &str) -> Result<String> {
    require(
        !value.is_empty()
            && value.trim_matches(is_python_space) == value
            && !value.contains(['\t', '\r', '\n']),
        PATH_TEXT,
    )?;
    require(
        !value.starts_with("/paedyl01/"),
        format!(
            "Unresolved old-host reference resource path: {value}; supply its actual local path"
        ),
    )?;
    let path = abspath(value)?;
    let roots = inspection_roots()?;
    if roots.is_empty() {
        return Ok(path);
    }
    let inside = |candidate: &str| roots.iter().any(|root| is_relative_to(candidate, root));
    require(
        inside(&path),
        format!("Reference resource path outside permitted inspection roots: {value}"),
    )?;
    let all = parts(&path);
    let mut cursor = all[0].to_owned();
    let mut pending = all[1..]
        .iter()
        .map(|part| (*part).to_owned())
        .collect::<VecDeque<_>>();
    let mut links = 0;
    while let Some(part) = pending.pop_front() {
        cursor = join(&cursor, &part);
        if is_symlink(&cursor)? {
            links += 1;
            require(
                links <= 40,
                format!("Too many symlinks in reference resource path: {value}"),
            )?;
            let link = fs::read_link(&cursor)
                .map_err(|failure| os_error(failure, Path::new(&cursor), None))?;
            let target = abspath(&join(&parent(&cursor), &link.to_string_lossy()))?;
            require(
                inside(&target),
                format!("Reference resource symlink leaves permitted inspection roots: {value}"),
            )?;
            let target_parts = parts(&target);
            for part in target_parts[1..].iter().rev() {
                pending.push_front((*part).to_owned());
            }
            cursor = target_parts[0].to_owned();
        }
    }
    Ok(path)
}

/// `checked_path(value if present(value) else default)`
fn checked_value(value: &Value, default: &str) -> Result<String> {
    if value.is_unset() {
        return checked_path(default);
    }
    match value {
        Value::Str(text) => checked_path(text),
        _ => Err(value_error(PATH_TEXT)),
    }
}

/// `destination(value, old_default, new_default)`: an absent exact legacy
/// default moves to the new default unless any companion state remains at
/// the old path. Existing old defaults stay managed: changing cache_root must
/// never relabel them as custom inputs.
fn destination(value: &Value, old: &str, new: &str) -> Result<(String, bool)> {
    let mut path = checked_value(value, new)?;
    if path == old && path != new && !lexists(&path) {
        for suffix in LEGACY_STATE {
            let artifact = checked_path(&format!("{path}{suffix}"))?;
            require(
                !lexists(&artifact),
                format!(
                    "Legacy managed reference has retained index/source/staging state: \
                     {artifact}; resolve the selected original destination explicitly before \
                     redirecting it"
                ),
            )?;
        }
        path = new.to_owned();
    }
    let managed = path == old || path == new;
    Ok((path, managed))
}

/// `os.path.samefile(first, second)`
fn same_file(first: &str, second: &str) -> Result<bool> {
    let one = metadata(Path::new(first))?;
    let other = metadata(Path::new(second))?;
    Ok(one.dev() == other.dev() && one.ino() == other.ino())
}

#[derive(Default)]
struct Planner {
    resources: Vec<Resource>,
    updates: Vec<(String, String)>,
}

impl Planner {
    fn add(
        &mut self,
        key: &str,
        value: &Value,
        assembly: &'static str,
        old: &str,
        new: &str,
    ) -> Result<String> {
        let (path, managed) = destination(value, old, new)?;
        self.updates.push((key.to_owned(), path.clone()));
        match self.resources.iter().find(|known| known.path == path) {
            Some(existing) => {
                require(existing.assembly == assembly, format!("{DISTINCT}: {path}"))?
            }
            None => self.resources.push(Resource {
                path: path.clone(),
                assembly,
                key: key.to_owned(),
                managed,
                kind: Kind::Fasta,
            }),
        }
        Ok(path)
    }
}

/// `resource_plan(config, data_dir, scope, use_cache_root=...)`: the planned
/// resources and the configuration updates naming them, in the script's order.
pub(crate) fn resource_plan(
    config: &Config,
    data_dir: &str,
    scope: Scope,
    use_cache_root: bool,
) -> Result<(Vec<Resource>, Vec<(String, String)>)> {
    let selected = assembly_alias(config.get("assembly"))?.ok_or_else(|| {
        value_error(
            "Set assembly to hg19/GRCh37 or hg38/GRCh38; reference filenames do not select an \
             assembly",
        )
    })?;
    let data_dir = checked_path(data_dir)?;
    let mut managed_dir = data_dir.clone();
    let root = config.get("cache_root");
    if use_cache_root && !root.is_unset() {
        let Value::Str(text) = root else {
            return Err(value_error(CACHE_ROOT));
        };
        require(cache_root_text(text), CACHE_ROOT)?;
        managed_dir = checked_path(text)?;
        require(
            parts(&managed_dir).len() > 1,
            "cache_root must name a dedicated cache directory",
        )?;
    }
    let fasta =
        |folder: &str, assembly: &str| join(folder, &format!("reference/ucsc.{assembly}.fasta"));
    let mut planner = Planner::default();
    let chosen = planner.add(
        "ref_genome",
        config.get("ref_genome"),
        selected,
        &fasta(&data_dir, selected),
        &fasta(&managed_dir, selected),
    )?;
    if scope != Scope::Selected {
        for assembly in ASSEMBLIES {
            let key = format!("gofcards_{assembly}_fasta");
            let mut value = config.get(&key).clone();
            if value.is_unset() && assembly == selected {
                value = Value::Str(chosen.clone());
            }
            planner.add(
                &key,
                &value,
                assembly,
                &fasta(&data_dir, assembly),
                &fasta(&managed_dir, assembly),
            )?;
        }
        // Catch distinct symbolic or hard links selecting the same wrong assembly.
        let resources = &planner.resources;
        for (number, first) in resources.iter().enumerate() {
            for second in &resources[number + 1..] {
                if exists(&first.path)?
                    && exists(&second.path)?
                    && same_file(&first.path, &second.path)?
                {
                    require(first.assembly == second.assembly, DISTINCT)?;
                }
            }
        }
    }
    if scope == Scope::Gofcards {
        let (path, managed) = destination(
            config.get(CHAIN_KEY),
            &join(&data_dir, CHAIN_FILE),
            &join(&managed_dir, CHAIN_FILE),
        )?;
        require(
            planner.resources.iter().all(|known| known.path != path),
            "The forward-chain destination must be distinct from every FASTA",
        )?;
        planner.resources.push(Resource {
            path: path.clone(),
            assembly: "hg19_to_hg38",
            key: CHAIN_KEY.to_owned(),
            managed,
            kind: Kind::Chain,
        });
        planner.updates.push((CHAIN_KEY.to_owned(), path));
    }
    Ok((planner.resources, planner.updates))
}

const CACHE_ROOT: &str = "cache_root must be an absolute path without spaces or shell characters";

/// `re.fullmatch(r"/[A-Za-z0-9_./-]+", root)`
fn cache_root_text(root: &str) -> bool {
    root.len() > 1
        && root.starts_with('/')
        && root
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || b"_./-".contains(&byte))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder_error::RuntimeError;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn message(failure: Box<dyn std::error::Error>) -> String {
        failure.downcast_ref::<RuntimeError>().unwrap().message()
    }

    fn folder() -> String {
        let path = std::env::temp_dir().join(format!(
            "priva-reference-plan-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        abspath(&path.display().to_string()).unwrap()
    }

    fn plan(text: &str, data: &str, scope: Scope, cache: bool) -> Result<Vec<Resource>> {
        let config = Config::load(text, false).unwrap();
        Ok(resource_plan(&config, data, scope, cache)?.0)
    }

    fn refused(text: &str, data: &str, scope: Scope, cache: bool) -> String {
        message(plan(text, data, scope, cache).unwrap_err())
    }

    #[test]
    fn path_texts_are_checked_before_any_file_test() {
        for bad in ["", " x", "x\t", "a\nb", "\u{1c}a"] {
            let failure = checked_path(bad).unwrap_err();
            assert_eq!(message(failure), PATH_TEXT, "{bad:?}");
        }
        let failure = checked_path("/paedyl01/x").unwrap_err();
        assert!(message(failure).starts_with("Unresolved old-host"));
        assert!(cache_root_text("/a/b_c.d-e"));
        for bad in ["/", "a", "/with space", "/\u{e9}"] {
            assert!(!cache_root_text(bad), "{bad}");
        }
    }

    #[test]
    fn blank_defaults_resolve_both_assemblies_and_the_chain() {
        let root = folder();
        let data = format!("{root}/data");
        let config = Config::load("assembly: GRCh38\n", false).unwrap();
        let (resources, updates) = resource_plan(&config, &data, Scope::Gofcards, false).unwrap();
        let paths: Vec<&str> = resources.iter().map(|r| r.path.as_str()).collect();
        let expected = [
            format!("{data}/reference/ucsc.hg38.fasta"),
            format!("{data}/reference/ucsc.hg19.fasta"),
            format!("{data}/liftOver/hg19ToHg38.over.chain"),
        ];
        assert_eq!(paths, expected);
        assert!(resources.iter().all(|r| r.managed));
        assert_eq!(resources[0].key, "ref_genome");
        let keys: Vec<&str> = updates.iter().map(|(key, _)| key.as_str()).collect();
        assert_eq!(keys[0], "ref_genome");
        assert_eq!(keys[1..3], ["gofcards_hg19_fasta", "gofcards_hg38_fasta"]);
        assert_eq!(keys[3], CHAIN_KEY);
        assert_eq!(updates[0].1, updates[2].1);
        assert!(!Path::new(&data).exists());
    }

    #[test]
    fn cache_root_routes_defaults_only_in_installer_mode() {
        let root = folder();
        let data = format!("{root}/data");
        let text = format!("assembly: hg19\ncache_root: {root}/cache\n");
        let routed = plan(&text, &data, Scope::Selected, true).unwrap();
        let expected = format!("{root}/cache/reference/ucsc.hg19.fasta");
        assert_eq!(routed[0].path, expected);
        let direct = plan(&text, &data, Scope::Selected, false).unwrap();
        assert_eq!(direct[0].path, format!("{data}/reference/ucsc.hg19.fasta"));
        for bad in ["relative", "/", "true"] {
            let text = format!("assembly: hg19\ncache_root: {bad}\n");
            let failure = refused(&text, &data, Scope::Selected, true);
            assert!(failure.starts_with("cache_root must"), "{bad}");
        }
    }

    #[test]
    fn legacy_defaults_move_only_without_retained_state() {
        let root = folder();
        let data = format!("{root}/data");
        let legacy = format!("{data}/reference/ucsc.hg38.fasta");
        let text = format!("assembly: hg38\ncache_root: {root}/cache\nref_genome: {legacy}\n");
        let moved = plan(&text, &data, Scope::Selected, true).unwrap();
        assert_eq!(
            moved[0].path,
            format!("{root}/cache/reference/ucsc.hg38.fasta")
        );
        assert!(moved[0].managed);
        fs::create_dir_all(format!("{data}/reference")).unwrap();
        fs::write(format!("{legacy}.install.lock"), "").unwrap();
        let failure = refused(&text, &data, Scope::Selected, true);
        assert!(failure.starts_with("Legacy managed reference has retained"));
    }

    #[test]
    fn one_path_cannot_serve_two_assemblies_or_the_chain() {
        let root = folder();
        let shared = format!("{root}/a.fa");
        let text = format!("assembly: hg38\nref_genome: {shared}\ngofcards_hg19_fasta: {shared}\n");
        let failure = refused(&text, &root, Scope::Both, false);
        assert_eq!(failure, format!("{DISTINCT}: {root}/a.fa"));
        let chain = format!("{root}/reference/ucsc.hg38.fasta");
        let text = format!("assembly: hg38\ngofcards_hg19_to_hg38_chain: {chain}\n");
        let failure = refused(&text, &root, Scope::Gofcards, false);
        assert!(failure.starts_with("The forward-chain destination"));
        let failure = refused("ref_genome: x.fa\n", &root, Scope::Selected, false);
        assert!(failure.contains("filenames do not select"));
        let failure = refused(
            "assembly: hg38\nref_genome: 5\n",
            &root,
            Scope::Selected,
            false,
        );
        assert_eq!(failure, PATH_TEXT);
        let explicit = format!("assembly: hg38\nref_genome: {root}/original.fa\n");
        let resources = plan(&explicit, &root, Scope::Selected, false).unwrap();
        assert!(!resources[0].managed);
    }
}
