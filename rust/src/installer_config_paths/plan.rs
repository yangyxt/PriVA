use super::*;
use files::{exists, is_dir, is_link};

pub(crate) fn checked_storage(value: Value<'_>, key: &str) -> Result<PathBuf> {
    let Some(s) = value.string().filter(|s| {
        s.starts_with('/')
            && s.len() > 1
            && s.bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"_./-".contains(&c))
    }) else {
        return Err(error(
            "ValueError",
            format!("Set {key} to an absolute storage path without spaces or shell characters"),
        ));
    };
    let path = tempfile_parent(Path::new(s))?;
    require(
        path != Path::new("/"),
        format!("{key} must name a dedicated storage directory"),
    )?;
    Ok(path)
}
pub(crate) fn validate_repeat_filename(value: Value<'_>) -> Result<()> {
    require(
        value.empty()
            || value.string().is_some_and(|v| {
                !v.is_empty() && v != "." && v != ".." && !v.contains(['/', '\\', '\r', '\n'])
            }),
        "repeat_region_file_name must be one filename beneath data/repeats, not an absolute or relative path",
    )
}
pub(crate) fn reverse_chain_source(
    ctx: &Context,
    value: Value<'_>,
    source_root: Option<&Path>,
) -> Result<PathBuf> {
    let root = match source_root {
        Some(p) => p.to_path_buf(),
        None => ctx.source_root()?,
    };
    if value.empty() {
        return Ok(join(&root, "data/liftover/hg38ToHg19.over.chain.gz"));
    }
    let s = value.string().ok_or_else(|| {
        error(
            "ValueError",
            "hg38_hg19_chain must be a bundled filename or an absolute supplied path",
        )
    })?;
    let p = path_value(s);
    if p.is_absolute() {
        return Ok(p);
    }
    require(
        folder_name(s),
        "hg38_hg19_chain must be a bundled filename or an absolute supplied path",
    )?;
    Ok(join(&root.join("data/liftover"), s))
}
pub(crate) fn resource_groups(ctx: &Context, config: &Config, dirs: &Paths) -> Result<Paths> {
    let root = checked_storage(config.index("cache_root")?, "cache_root")?;
    let plugin = path_value(get(dirs, "vep_plugins_cachedir")?);
    let mut groups: Vec<(PathBuf, Vec<&str>)> = vec![];
    let mut add = |p: PathBuf, keys: &[&'static str]| {
        if let Some((_, v)) = groups.iter_mut().find(|(k, _)| *k == p) {
            *v = keys.to_vec();
        } else {
            groups.push((p, keys.to_vec()));
        }
    };
    add(
        join(&plugin, "AlphaMissense"),
        &[
            "alphamissense_prescore",
            "alphamissense_vcf",
            "alphamissense_vep_vcf",
            "alphamissense_vep_base_vcf",
            "alphamissense_pd_stat",
            "alphamissense_tranx_domain_map",
            "alphamissense_intolerant_domains",
            "alphamissense_intolerant_motifs",
            "all_intolerant_domains",
            "pm1_regions_pkl",
        ],
    );
    add(
        path_value(get(dirs, "clinvar_vcf_dir")?),
        &[
            "clinvar_vcf",
            "clinvar_aa_stat",
            "clinvar_splice_stat",
            "clinvar_patho_af_stat",
            "clinvar_patho_exon_af_stat",
            "clinvar_pm5_index",
            "clinvar_gene_stat",
        ],
    );
    add(
        path_value(get(dirs, "gnomad_vcf_dir")?),
        &["gnomad_vcf_chrX"],
    );
    add(join(&plugin, "UTRannotator"), &["utr_annotator_file"]);
    add(
        join(&plugin, "LoFtee"),
        &[
            "human_ancestor_fasta",
            "loftee_conservation_file",
            "gerp_bigwig",
        ],
    );
    add(
        join(Path::new(get(dirs, "vep_plugins_dir")?), "LoFtee_repos"),
        &["loftee_parent_dir", "loftee_repo"],
    );
    for (folder, keys) in [
        ("Conservation", &["conservation_file"][..]),
        ("SpliceVault", &["splicevault_prescore"]),
        ("PrimateAI", &["primateai_prescore"]),
        (
            "SpliceAI",
            &["spliceai_snv_prescore", "spliceai_indel_prescore"],
        ),
        ("Ensembl_CDS", &["cds_fasta_file"]),
        ("MaveDB", &["mavedb_file"]),
        ("pext", &["pext_mean_bw", "pext_tissue_bw_dir"]),
    ] {
        add(join(&plugin, folder), keys);
    }
    add(
        path_value(get(dirs, "cadd_parent_dir")?),
        &["cadd_base_dir"],
    );
    for (folder, keys) in [
        ("InterPro", &["interpro_mapping_pickle"][..]),
        ("ClinGen", &["clingen_map", "gene_dosage_sensitivity"]),
        ("LOEUF", &["loeuf_prescore"]),
        (
            "reference",
            &["ref_genome", "gofcards_hg19_fasta", "gofcards_hg38_fasta"],
        ),
        (
            "liftOver",
            &["gofcards_hg19_to_hg38_chain", "hg38_hg19_chain"],
        ),
        (
            "patho_mechanism",
            &[
                "gene_mechanism_cache_dir",
                "gene_mechanism_raw_dir",
                "gene_nonlof_mechanism_json",
                "hpo_condition_mechanism_json",
                "ddg2p_mechanism_evidence",
                "clinvar_vcv_match_audit",
            ],
        ),
        (
            "gofcards",
            &["gofcards_workdir", "gofcards_exact_gof_cache"],
        ),
        ("patho_mechanism/raw/clinvar_vcv", &["clinvar_vcv_xml"]),
        ("hpo", &["hpo_raw_dir", "hpo_assertions"]),
        ("mondo", &["mondo_obo", "mondo_disease_scope_registry"]),
    ] {
        add(join(&root, folder), keys);
    }
    let mut result = vec![];
    for (folder, keys) in groups {
        for key in keys {
            set(&mut result, key, folder.display().to_string());
        }
    }
    for key in crate::python_set_order::StringSet::literal(ctx.seed, REPLACEMENT)
        .difference_small(&["hgnc_table"])
    {
        set(
            &mut result,
            &key,
            join(&root, Path::new(bundled_path(&key)).parent().unwrap())
                .display()
                .to_string(),
        );
    }
    Ok(result)
}
pub(crate) fn bundled_resource(
    ctx: &Context,
    key: &str,
    value: &str,
    source_root: Option<&Path>,
) -> Result<Option<PathBuf>> {
    let Some(root) = source_root else {
        return Ok(None);
    };
    let mut options = vec![];
    if BUNDLED_CACHES.iter().any(|(k, _)| *k == key) || REPLACEMENT.contains(&key) {
        options.push(bundled_path(key).to_owned());
    }
    match key {
        "interpro_mapping_pickle" => {
            options.push("data/InterPro/Interpro_entry_mapping.pkl.gz".into())
        }
        "clinvar_gene_stat" => options.push("data/ClinVar/clinvar_2star_stats.pkl.gz".into()),
        "hg38_hg19_chain" => options.push("data/liftover/hg38ToHg19.over.chain.gz".into()),
        "hgnc_table" => options.push("data/hgnc/non_alt_loci_set.tsv".into()),
        "loeuf_prescore" => {
            options = vec![
                "data/loeuf/loeuf_dataset.tsv.gz".into(),
                "data/loeuf/loeuf_dataset_hg38.tsv.gz".into(),
            ]
        }
        "clingen_map" | "gene_dosage_sensitivity" => {
            options = ["hg19", "hg38"]
                .iter()
                .map(|a| {
                    format!(
                        "data/clingen/{key}.{a}.{}",
                        if key == "clingen_map" {
                            "pkl.gz"
                        } else {
                            "tsv"
                        }
                    )
                })
                .collect()
        }
        _ => {}
    }
    let selected = if key == "hg38_hg19_chain" {
        reverse_chain_source(ctx, Value::Text(value), Some(root))?
    } else {
        path_value(value)
    };
    for relative in options {
        let bundled = join(root, &relative);
        let mut expected = join(&resolve(root)?, &relative);
        if key == "hgnc_table"
            && relative == "data/hgnc/non_alt_loci_set.tsv"
            && is_link(&bundled)?
            && crate::python_unicode::from_native(
                &std::fs::read_link(files::native(&bundled)?)
                    .map_err(|e| os_error(e, &bundled, None))?,
            ) == Path::new(bundled_path(key)).file_name().unwrap()
        {
            expected = join(&resolve(root)?, bundled_path(key));
        }
        let legacy = matches!(key, "clinvar_gene_stat" | "interpro_mapping_pickle")
            && selected.components().count() == 1
            && name(&selected) == name(&bundled)
            && !is_file(&selected)?;
        let alias = selected == Path::new(&relative) || legacy;
        let mut same = resolve(Path::new(&absolute(&selected)?))? == resolve(&bundled)?;
        if key == "hgnc_table" {
            let logical = if selected.is_absolute() {
                selected.clone()
            } else {
                join(root, &selected)
            };
            same = tempfile_parent(&logical)? == tempfile_parent(&bundled)?;
        }
        if (alias || same)
            && resolve(&bundled)? == expected
            && is_file(&bundled)?
            && metadata(&bundled)?.len() > 0
        {
            return Ok(Some(bundled));
        }
    }
    Ok(None)
}
pub(crate) fn clinvar_default_marker(
    value: &str,
    source_root: Option<&Path>,
) -> Result<Option<PathBuf>> {
    let Some(root) = source_root else {
        return Ok(None);
    };
    let relative = Path::new("data/ClinVar/clinvar_2star_stats.pkl.gz");
    let baseline = join(root, relative);
    let selected = path_value(value);
    let legacy = selected == Path::new(&name(&baseline)) && !is_file(&selected)?;
    Ok(
        ((selected == relative || absolute(&selected)? == absolute(&baseline)? || legacy)
            && resolve(&baseline)? == join(&resolve(root)?, relative))
        .then_some(baseline),
    )
}
pub(crate) fn source_is_needed(key: &str, config: &Config) -> bool {
    let old = matches!(config.get("assembly").string(), Some("hg19" | "GRCh37"));
    match key {
        "hg38_hg19_chain" => old,
        "pext_mean_bw" | "pext_tissue_bw_dir" => {
            config.get("pext_enabled").text().to_lowercase() == "true"
        }
        "mavedb_file" if old => false,
        _ => true,
    }
}
pub(crate) fn defer_missing_build_source(
    key: &str,
    source: &Path,
    root: Option<&Path>,
) -> Result<bool> {
    let Some(root) = root.filter(|_| key == "gofcards_mechanism_review_tsv") else {
        return Ok(false);
    };
    let relative = Path::new(bundled_path(key));
    if source == relative || tempfile_parent(source)? == tempfile_parent(&join(root, relative))? {
        return Ok(false);
    };
    Ok(!exists(source)?
        && !is_link(source)?
        && !["PRIVA_FORCE_GOFCARDS_CACHE", "PRIVA_FORCE_ALL_CACHES"]
            .iter()
            .any(|v| std::env::var(v).ok().as_deref() == Some("1")))
}
fn splice_failure(f: crate::spliceai_cache_input::Failure) -> Box<dyn std::error::Error> {
    use crate::spliceai_cache_input::Failure;
    match f {
        Failure::Rejected(s) => {
            let kind = if s.starts_with("[Errno 2]") {
                "FileNotFoundError"
            } else if s.starts_with("[Errno 13]") {
                "PermissionError"
            } else if s.starts_with("[Errno 21]") {
                "IsADirectoryError"
            } else if s.starts_with("[Errno ") {
                "OSError"
            } else {
                "ValueError"
            };
            error(kind, crate::python_unicode::encode(&s))
        }
        Failure::Crash(s) => error("RuntimeError", crate::python_unicode::encode(&s)),
        Failure::Helper(n) => error("CalledProcessError", format!("helper exited {n}")),
    }
}
pub(crate) fn source_format(ctx: &Context, key: &str, value: &str) -> Result<Option<&'static str>> {
    if !matches!(key, "spliceai_snv_prescore" | "spliceai_indel_prescore")
        || !is_file(Path::new(value))?
    {
        return Ok(None);
    }
    ctx.splice()?;
    Ok(Some(
        crate::spliceai_cache_input::inspect_header(&files::native(Path::new(value))?)
            .map_err(splice_failure)?
            .format,
    ))
}
pub(super) fn prepared(ctx: &Context, key: &str, value: &str) -> Result<bool> {
    if source_format(ctx, key, value)? != Some("vcf_spliceai2") {
        return Ok(false);
    }
    require(
        is_file(
            &resolve(Path::new(value))?
                .parent()
                .unwrap()
                .join("adapter.ready.json"),
        )?,
        "Prepared SpliceAI2 source requires its completed adjacent adapter.ready.json; ordinary VCF copying is not a substitute",
    )?;
    Ok(true)
}
pub(crate) fn resolve_cache_paths(
    ctx: &Context,
    config: &Config,
    source_root: Option<&Path>,
) -> Result<(Paths, Paths)> {
    require(config.mapping(), "The configuration must be a YAML mapping")?;
    validate_repeat_filename(config.get("repeat_region_file_name"))?;
    if config.get("cache_root").empty() {
        return Ok((vec![], vec![]));
    }
    let root = checked_storage(config.get("cache_root"), "cache_root")?;
    let folders = config.get("cache_folders");
    require(
        folders.null()
            || folders.mapping().is_some_and(|v| {
                v.iter()
                    .all(|(k, _)| k.is_some_and(|k| FOLDERS.iter().any(|(n, _, _)| *n == k)))
            }),
        "cache_folders must be a mapping with only the six supported resource names",
    )?;
    require(
        !config
            .keys()
            .iter()
            .any(|k| FOLDERS.iter().any(|(n, _, _)| *n == k)),
        "Indent resource folder names beneath cache_folders",
    )?;
    let (mut directories, mut previous) = (vec![], vec![]);
    for (resource, default, key) in FOLDERS {
        let selected = folders.get(resource);
        let selected = if matches!(selected, Value::Missing) {
            Value::Text(default)
        } else {
            selected
        };
        require(
            selected.string().is_some_and(folder_name),
            format!("cache_folders.{resource} must be one folder name"),
        )?;
        let dest = join(&root, selected.string().unwrap());
        set(&mut directories, key, dest.display().to_string());
        if let Some(old) = config
            .get(key)
            .string()
            .filter(|v| v.starts_with('/') && path_value(v) != Path::new("/"))
        {
            previous.push((tempfile_parent(Path::new(old))?, dest));
        }
    }
    previous.sort_by_key(|(old, _)| std::cmp::Reverse(old.components().count()));
    let (snv, indel) = (
        config.get("spliceai_snv_prescore"),
        config.get("spliceai_indel_prescore"),
    );
    if !config.get("assembly").null() {
        if let (Some(snv), Some(indel)) = (snv.string(), indel.string()) {
            if is_file(Path::new(snv))? && is_file(Path::new(indel))? {
                ctx.splice()?;
                let a = config.get("assembly").string().unwrap_or("");
                crate::spliceai_cache_input::validate_pair_input(
                    a,
                    &files::native(Path::new(snv))?,
                    &files::native(Path::new(indel))?,
                )
                .map_err(splice_failure)?;
            }
        }
    }
    let mut updates = directories.clone();
    for key in ["nm_to_enst_tsv", "hgnc_table"] {
        let mapping = config.get(key);
        if !mapping.empty() {
            let s = mapping
                .string()
                .ok_or_else(|| error("ValueError", format!("{key} must be a path string")))?;
            if let Some(r) = source_root.filter(|_| !Path::new(s).is_absolute()) {
                set(
                    &mut updates,
                    key,
                    tempfile_parent(&join(r, s))?.display().to_string(),
                );
            }
        }
    }
    for (key, folder) in resource_groups(ctx, config, &directories)? {
        let value = config.get(&key);
        if value.empty() {
            continue;
        }
        let value = value
            .string()
            .ok_or_else(|| error("ValueError", format!("{key} must be a path string")))?;
        let value = if key == "hg38_hg19_chain" {
            reverse_chain_source(ctx, Value::Text(value), source_root)?
                .display()
                .to_string()
        } else {
            value.into()
        };
        if !source_is_needed(&key, config)
            || defer_missing_build_source(&key, Path::new(&value), source_root)?
        {
            continue;
        }
        if source_format(ctx, &key, &value)? == Some("spliceai2_tsv") {
            continue;
        }
        if prepared(ctx, &key, &value)? {
            if resolve(Path::new(&value))?.starts_with(resolve(&root)?) {
                continue;
            }
            set(
                &mut updates,
                &key,
                join(
                    &path_value(get(&directories, "vep_plugins_cachedir")?)
                        .join("SpliceAI2/GRCh38"),
                    name(&resolve(Path::new(&value))?),
                )
                .display()
                .to_string(),
            );
            continue;
        }
        if key == "clinvar_gene_stat" {
            if let Some(b) = clinvar_default_marker(&value, source_root)? {
                set(&mut updates, &key, b.display().to_string());
                continue;
            }
        }
        if let Some(b) = bundled_resource(ctx, &key, &value, source_root)? {
            set(&mut updates, &key, b.display().to_string());
            continue;
        }
        let path = tempfile_parent(Path::new(&value))?;
        if matches!(key.as_str(), "alphamissense_prescore" | "gnomad_vcf_chrX")
            && !path.starts_with(&root)
        {
            set(
                &mut updates,
                &key,
                join(
                    &path_value(&folder).join("supplied").join(&key),
                    name(&path),
                )
                .display()
                .to_string(),
            );
            continue;
        }
        if let Some((old, new)) = previous.iter().find(|(old, _)| path.starts_with(old)) {
            set(
                &mut updates,
                &key,
                join(new, path.strip_prefix(old).unwrap())
                    .display()
                    .to_string(),
            );
        } else if !path.starts_with(&root) {
            set(
                &mut updates,
                &key,
                join(
                    &path_value(&folder).join("supplied").join(&key),
                    name(&path),
                )
                .display()
                .to_string(),
            );
        }
    }
    Ok((directories, updates))
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Adoption {
    pub key: String,
    pub source: PathBuf,
    pub destination: PathBuf,
    pub kind: String,
}
impl Adoption {
    pub fn json(&self) -> Json {
        Json::Object(vec![
            ("key".into(), Json::Str(self.key.clone())),
            (
                "source".into(),
                Json::Str(self.source.display().to_string()),
            ),
            (
                "destination".into(),
                Json::Str(self.destination.display().to_string()),
            ),
            ("kind".into(), Json::Str(self.kind.clone())),
        ])
    }
}
pub(crate) fn adoption_plan(
    ctx: &Context,
    config: &Config,
    updates: &Paths,
    root: &Path,
) -> Result<Vec<Adoption>> {
    for key in ["nm_to_enst_tsv", "hgnc_table"] {
        let m = config.get(key);
        if !m.empty() {
            let value = m
                .string()
                .ok_or_else(|| error("ValueError", format!("{key} must be a path string")))?;
            let selected = if Path::new(value).is_absolute() {
                path_value(value)
            } else {
                join(root, value)
            };
            require(
                is_file(&selected)? && metadata(&selected)?.len() > 0,
                format!(
                    "Explicit {key} source is missing/empty; refusing to substitute another dataset/release: {}",
                    selected.display()
                ),
            )?;
        }
    }
    let dirs = FOLDERS
        .iter()
        .map(|(_, _, k)| {
            Ok((
                k.to_string(),
                get(updates, k)
                    .map(str::to_owned)
                    .or_else(|_| config.get(k).path().map(|p| p.display().to_string()))?,
            ))
        })
        .collect::<Result<Paths>>()?;
    let keys = resource_groups(ctx, config, &dirs)?;
    let automatic = defaults::automatic_source_paths(ctx, config, updates, root)?;
    let mut candidates = updates.clone();
    for key in config.keys() {
        if keys.iter().any(|(k, _)| *k == key) {
            if let Ok(v) = get(updates, &key) {
                set(&mut candidates, &key, v);
            } else {
                set(&mut candidates, &key, config.get(&key).text());
            }
        }
    }
    let mut plan = vec![];
    for (key, destination) in candidates {
        if !keys.iter().any(|(k, _)| *k == key) && !directory_key(&key) {
            continue;
        }
        let v = config.get(&key);
        let value = if key == "hg38_hg19_chain" && !v.empty() {
            reverse_chain_source(ctx, v, Some(root))?
                .display()
                .to_string()
        } else {
            match v.string() {
                Some(s) if !v.empty() => s.into(),
                _ => continue,
            }
        };
        if source_format(ctx, &key, &value)? == Some("spliceai2_tsv")
            || !source_is_needed(&key, config)
        {
            continue;
        }
        if key == "clinvar_gene_stat" && clinvar_default_marker(&value, Some(root))?.is_some() {
            continue;
        }
        if bundled_resource(ctx, &key, &value, Some(root))?.is_some() {
            continue;
        }
        let source = tempfile_parent(Path::new(&value))?;
        let target = path_value(&destination);
        if defer_missing_build_source(&key, &source, Some(root))? {
            continue;
        }
        let present = if directory_key(&key) {
            is_dir(&source)?
        } else {
            is_file(&source)?
        };
        if !present {
            require(
                !exists(&source)? && !is_link(&source)?,
                format!(
                    "Explicit {key} source has the wrong file type or is a broken link: {}",
                    source.display()
                ),
            )?;
            if generated_directory(&key)
                || (GENERATED.contains(&key.as_str()) && key != "clinvar_gene_stat")
                || automatic
                    .iter()
                    .any(|(k, paths)| k == &key && paths.contains(&source))
            {
                continue;
            }
            return Err(error(
                "ValueError",
                format!(
                    "Explicit {key} source is missing; refusing to substitute another dataset/release: {}",
                    source.display()
                ),
            ));
        }
        if key == "gnomad_vcf_chrX" {
            ctx.import(
                "gnomad_cache.py",
                include_bytes!("../../../scripts/gnomad_cache.py"),
            )?;
            let assembly = config.index("assembly")?;
            let canonical = gnomad::default_path(
                Path::new(get(updates, "gnomad_vcf_dir")?),
                assembly.string().unwrap_or(""),
            )?;
            if source != canonical {
                gnomad::validate_family(
                    &source,
                    assembly.string().unwrap_or(""),
                    config.get("ref_genome"),
                )?;
            }
            if source != target && resolve(&source)? != resolve(&target)? {
                let a = gnomad::family_paths(&source)?;
                let b = gnomad::family_paths(&target)?;
                for ((_, member), (_, dest)) in a.iter().zip(&b) {
                    append_file(&mut plan, &key, member, dest)?;
                }
            }
            continue;
        }
        if source == target || resolve(&source)? == resolve(&target)? {
            continue;
        }
        if prepared(ctx, &key, &value)? {
            let item = Adoption {
                key: "spliceai2_prepared_pair".into(),
                source: resolve(&source)?.parent().unwrap().to_path_buf(),
                destination: target.parent().unwrap_or(Path::new(".")).to_path_buf(),
                kind: "prepared_projection".into(),
            };
            if !plan.contains(&item) {
                plan.push(item);
            }
            continue;
        }
        if directory_key(&key) {
            if !is_dir(&source)? || resolve(&source)?.starts_with(resolve(root)?.join("data")) {
                continue;
            }
            plan.push(Adoption {
                key,
                source,
                destination: target,
                kind: "directory".into(),
            });
        } else if is_file(&source)? {
            append_file(&mut plan, &key, &source, &target)?;
        }
    }
    let mut directories = plan
        .iter()
        .filter(|i| i.kind == "directory")
        .collect::<Vec<_>>();
    directories.sort_by_key(|i| i.source.components().count());
    let mut result = vec![];
    for item in &plan {
        if directories.iter().any(|parent| {
            item.kind != "prepared_projection"
                && item != *parent
                && item
                    .source
                    .strip_prefix(&parent.source)
                    .is_ok_and(|r| item.destination == join(&parent.destination, r))
        }) {
            continue;
        }
        if !result.contains(item) {
            result.push(item.clone());
        }
    }
    result.sort_by_key(|item| item.kind != "prepared_projection");
    Ok(result)
}
fn append_file(plan: &mut Vec<Adoption>, key: &str, source: &Path, target: &Path) -> Result<()> {
    plan.push(Adoption {
        key: key.into(),
        source: source.to_path_buf(),
        destination: target.to_path_buf(),
        kind: "file".into(),
    });
    for suffix in COMPANIONS {
        let p = files::suffix(source, suffix);
        if is_file(&p)? {
            plan.push(Adoption {
                key: key.into(),
                source: p,
                destination: files::suffix(target, suffix),
                kind: "file".into(),
            });
        }
    }
    Ok(())
}
pub(super) fn adopt_prepared(
    ctx: &Context,
    source: &Path,
    target: &Path,
    root: &Path,
) -> Result<()> {
    ctx.splice()?;
    crate::spliceai_cache_input::adopt_completed_projection(
        &files::native(source)?,
        &files::native(target)?,
        &files::native(root)?,
    )
    .map_err(splice_failure)?;
    Ok(())
}
