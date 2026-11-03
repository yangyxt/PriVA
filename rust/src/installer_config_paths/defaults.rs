use super::*;
use plan::{checked_storage, resource_groups};

fn fill(config: &mut Config, updates: &mut Paths, key: &str, value: impl Into<String>) {
    if config.get(key).empty() {
        let value = value.into();
        set(updates, key, value.clone());
        config.set(key, value);
    }
}
fn pathfill(config: &mut Config, updates: &mut Paths, key: &str, value: PathBuf) {
    fill(config, updates, key, value.display().to_string());
}
pub(crate) fn bundled_cache_default(
    key: &str,
    fallback: &Path,
    root: &Path,
    relative: Option<&str>,
) -> Result<PathBuf> {
    let relative = relative.unwrap_or_else(|| bundled_path(key));
    let bundled = join(root, relative);
    if !forced(key)
        && resolve(&bundled)? == join(&resolve(root)?, relative)
        && is_file(&bundled)?
        && metadata(&bundled)?.len() > 0
    {
        Ok(bundled)
    } else {
        Ok(fallback.to_path_buf())
    }
}
pub(crate) fn portable_defaults(config: &Config, source_root: &Path) -> Result<Paths> {
    require(config.mapping(), "The configuration must be a YAML mapping")?;
    let root = checked_storage(config.get("cache_root"), "cache_root")?;
    checked_storage(config.get("tmp_dir"), "tmp_dir")?;
    let a = config.get("assembly");
    if a.unhashable() {
        return Err(error(
            "TypeError",
            format!("unhashable type: '{}'", a.kind()),
        ));
    }
    let assembly = match a.string() {
        Some("hg19" | "GRCh37") => "hg19",
        Some("hg38" | "GRCh38") => "hg38",
        _ => {
            return Err(error(
                "ValueError",
                "Set assembly to hg19/GRCh37 or hg38/GRCh38",
            ));
        }
    };
    let version = config.get("vep_installed_version").text();
    require(
        version
            .bytes()
            .next()
            .is_some_and(|c| (b'1'..=b'9').contains(&c))
            && version.bytes().all(|c| c.is_ascii_digit()),
        "Set vep_installed_version to the requested numeric major release",
    )?;
    let f = config.get("cache_folders");
    let folders = if f.truth() { f } else { Value::Missing };
    require(
        matches!(folders, Value::Missing)
            || folders.mapping().is_some_and(|v| {
                v.iter()
                    .all(|(k, _)| k.is_some_and(|k| FOLDERS.iter().any(|(n, _, _)| *n == k)))
            }),
        "cache_folders must contain only the six supported resource names",
    )?;
    let (mut updates, mut effective) = (vec![], config.clone());
    for (resource, name, key) in FOLDERS {
        let v = folders.get(resource);
        let v = if matches!(v, Value::Missing) {
            Value::Text(name)
        } else {
            v
        };
        require(
            v.string().is_some_and(folder_name),
            format!("cache_folders.{resource} must be one folder name"),
        )?;
        pathfill(
            &mut effective,
            &mut updates,
            key,
            join(&root, v.string().unwrap()),
        );
    }
    for (key, relative) in BUNDLED_PATHS {
        let selected = effective.get(key);
        if selected.empty() {
            pathfill(
                &mut effective,
                &mut updates,
                key,
                join(source_root, relative),
            );
        } else if let Some(selected) = selected.string() {
            if !Path::new(selected).is_absolute() {
                let value = join(source_root, selected).display().to_string();
                set(&mut updates, key, value.clone());
                effective.set(key, value);
            }
        } else {
            return Err(error("ValueError", format!("{key} must be a path string")));
        }
    }
    let plugin = effective.index("vep_plugins_cachedir")?.path()?;
    let clinvar = effective.index("clinvar_vcf_dir")?.path()?;
    let genome = if assembly == "hg19" {
        "GRCh37"
    } else {
        "GRCh38"
    };
    let raw_assembly = config.index("assembly")?.text();
    let alpha = join(&plugin, "AlphaMissense");
    let gnomad = join(
        &effective.index("gnomad_vcf_dir")?.path()?,
        if assembly == "hg19" {
            "gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz"
        } else {
            "gnomad.joint.v4.1.sites.chrX.vcf.bgz"
        },
    );
    pathfill(&mut effective, &mut updates, "gnomad_vcf_chrX", gnomad);
    pathfill(
        &mut effective,
        &mut updates,
        "clinvar_vcf",
        join(&clinvar, format!("clinvar.{assembly}.vep.vcf.gz")),
    );
    for (key, suffix) in [
        ("clinvar_aa_stat", "aa_change.stats.pkl"),
        ("clinvar_splice_stat", "splice_change.stats.pkl"),
        ("clinvar_patho_af_stat", "patho_af_stat.pkl"),
        ("clinvar_patho_exon_af_stat", "patho_exon_af_stat.pkl"),
        ("clinvar_pm5_index", "pm5_index.tsv.gz"),
    ] {
        let v = effective.index("clinvar_vcf")?.text();
        fill(
            &mut effective,
            &mut updates,
            key,
            format!("{}.{suffix}", v.strip_suffix(".vcf.gz").unwrap_or(&v)),
        );
    }
    pathfill(
        &mut effective,
        &mut updates,
        "clinvar_gene_stat",
        join(source_root, "data/ClinVar/clinvar_2star_stats.pkl.gz"),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "utr_annotator_file",
        join(
            &plugin,
            format!("UTRannotator/uORF_5UTR_{genome}_PUBLIC.txt"),
        ),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "loeuf_prescore",
        join(
            &source_root.join("data/loeuf"),
            if assembly == "hg19" {
                "loeuf_dataset.tsv.gz"
            } else {
                "loeuf_dataset_hg38.tsv.gz"
            },
        ),
    );
    let parent = join(&effective.index("vep_plugins_dir")?.path()?, "LoFtee_repos");
    pathfill(&mut effective, &mut updates, "loftee_parent_dir", parent);
    let repo = join(
        &effective.index("loftee_parent_dir")?.path()?,
        format!("loftee-{raw_assembly}"),
    );
    pathfill(&mut effective, &mut updates, "loftee_repo", repo);
    let loftee = join(&plugin.join("LoFtee"), &raw_assembly);
    pathfill(
        &mut effective,
        &mut updates,
        "human_ancestor_fasta",
        join(&loftee, "human_ancestor.fa.gz"),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "loftee_conservation_file",
        join(&loftee, "phylocsf_gerp.sql"),
    );
    if assembly == "hg38" {
        pathfill(
            &mut effective,
            &mut updates,
            "gerp_bigwig",
            join(&loftee, "gerp_conservation_scores.homo_sapiens.bw"),
        );
    }
    for (key, name) in [
        (
            "alphamissense_prescore",
            format!("AlphaMissense_{assembly}.tsv.gz"),
        ),
        (
            "alphamissense_vcf",
            format!("AlphaMissense_{assembly}.vcf.gz"),
        ),
        (
            "alphamissense_vep_vcf",
            format!("AlphaMissense_{assembly}.vep.vcf.gz"),
        ),
        (
            "alphamissense_vep_base_vcf",
            format!("AlphaMissense_{assembly}.vep.base.vcf.gz"),
        ),
        (
            "alphamissense_pd_stat",
            format!("AlphaMissense_{assembly}.vep.prot.domain.stats.pkl"),
        ),
        (
            "alphamissense_tranx_domain_map",
            format!("transcript_pp_domain_mapping.{assembly}.pkl"),
        ),
        (
            "alphamissense_intolerant_domains",
            format!("domain_tolerance_analysis.{assembly}.tsv"),
        ),
        (
            "alphamissense_intolerant_motifs",
            format!("AlphaMissense_{assembly}.kde.pkl"),
        ),
        (
            "all_intolerant_domains",
            format!("all_intolerant_domains.{assembly}.pkl"),
        ),
        ("pm1_regions_pkl", format!("pm1_regions.{assembly}.pkl")),
    ] {
        pathfill(&mut effective, &mut updates, key, join(&alpha, name));
    }
    pathfill(
        &mut effective,
        &mut updates,
        "conservation_file",
        join(
            &plugin,
            format!("Conservation/gerp_conservation_scores.homo_sapiens.{genome}.bw"),
        ),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "splicevault_prescore",
        join(
            &plugin,
            format!("SpliceVault/SpliceVault_data_{genome}.tsv.gz"),
        ),
    );
    if assembly == "hg38" {
        pathfill(
            &mut effective,
            &mut updates,
            "mavedb_file",
            join(
                &plugin,
                format!("MaveDB/MaveDB_variants.{raw_assembly}.tsv.gz"),
            ),
        );
    }
    if config.get("pext_enabled").text().to_lowercase() == "true" {
        pathfill(
            &mut effective,
            &mut updates,
            "pext_mean_bw",
            join(
                &plugin,
                format!(
                    "pext/{assembly}/{}",
                    if assembly == "hg19" {
                        "mean_proportion.bw"
                    } else {
                        "exp_prop_mean.bw"
                    }
                ),
            ),
        );
        pathfill(
            &mut effective,
            &mut updates,
            "pext_tissue_bw_dir",
            join(&plugin, format!("pext/{assembly}")),
        );
    }
    let cadd = join(
        &effective.index("cadd_parent_dir")?.path()?,
        "CADD-scripts-1.7.2",
    );
    pathfill(&mut effective, &mut updates, "cadd_base_dir", cadd);
    pathfill(
        &mut effective,
        &mut updates,
        "hg38_hg19_chain",
        join(source_root, "data/liftover/hg38ToHg19.over.chain.gz"),
    );
    fill(
        &mut effective,
        &mut updates,
        "repeat_region_file_name",
        format!("ucsc.{assembly}.fasta.repeat.bed.gz"),
    );
    let interpro = join(source_root, "data/InterPro/Interpro_entry_mapping.pkl.gz");
    let interpro = if resolve(&interpro)?.starts_with(resolve(source_root)?.join("data/InterPro"))
        && is_file(&interpro)?
        && metadata(&interpro)?.len() > 0
    {
        interpro
    } else {
        join(&root, "InterPro/Interpro_entry_mapping.pkl.gz")
    };
    pathfill(
        &mut effective,
        &mut updates,
        "interpro_mapping_pickle",
        interpro,
    );
    for (key, suffix) in [
        ("clingen_map", "pkl.gz"),
        ("gene_dosage_sensitivity", "tsv"),
    ] {
        let relative = format!("data/clingen/{key}.{assembly}.{suffix}");
        pathfill(
            &mut effective,
            &mut updates,
            key,
            bundled_cache_default(
                key,
                &join(&root, format!("ClinGen/{key}.{assembly}.{suffix}")),
                source_root,
                Some(&relative),
            )?,
        );
    }
    pathfill(
        &mut effective,
        &mut updates,
        "cds_fasta_file",
        join(
            &plugin,
            format!("Ensembl_CDS/release-{version}/Homo_sapiens.{genome}.cds.all.fa.gz"),
        ),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "gene_mechanism_cache_dir",
        join(&root, "patho_mechanism"),
    );
    let mechanism = effective.index("gene_mechanism_cache_dir")?.path()?;
    pathfill(
        &mut effective,
        &mut updates,
        "gene_mechanism_raw_dir",
        join(&mechanism, "raw"),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "clinvar_vcv_match_audit",
        join(&mechanism, "clinvar_vcv_gofcards_matches.tsv"),
    );
    for (key, name) in [
        (
            "gene_nonlof_mechanism_json",
            "gene_nonlof_mechanism_curated_assertions.json.gz",
        ),
        (
            "hpo_condition_mechanism_json",
            "hpo_condition_mechanism_cache.json.gz",
        ),
        (
            "ddg2p_mechanism_evidence",
            "gene_pathogenic_mechanism_evidence.tsv",
        ),
    ] {
        pathfill(
            &mut effective,
            &mut updates,
            key,
            bundled_cache_default(key, &join(&mechanism, name), source_root, None)?,
        );
    }
    pathfill(
        &mut effective,
        &mut updates,
        "gofcards_exact_gof_cache",
        bundled_cache_default(
            "gofcards_exact_gof_cache",
            &join(&root, "gofcards/gofcards_exact_gof.json.gz"),
            source_root,
            None,
        )?,
    );
    pathfill(
        &mut effective,
        &mut updates,
        "gofcards_workdir",
        join(&root, "gofcards/build_work"),
    );
    let vcv = join(
        &effective.index("gene_mechanism_raw_dir")?.path()?,
        "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz",
    );
    pathfill(&mut effective, &mut updates, "clinvar_vcv_xml", vcv);
    let h = config.get("hpo_release");
    let h = if h.truth() {
        h
    } else {
        Value::Text("v2026-06-23")
    };
    let m = config.get("mondo_release");
    let m = if m.truth() {
        m
    } else {
        Value::Text("v2026-07-06")
    };
    let h = h.string().ok_or_else(|| {
        error(
            "TypeError",
            format!(
                "unsupported operand type(s) for /: 'PosixPath' and '{}'",
                h.kind()
            ),
        )
    })?;
    pathfill(
        &mut effective,
        &mut updates,
        "hpo_raw_dir",
        join(&root.join("hpo/raw"), h),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "hpo_assertions",
        join(&root, "hpo/genes_to_phenotype.assertions.tsv.gz"),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "mondo_obo",
        join(&root, format!("mondo/mondo-simple.{}.obo", m.text())),
    );
    pathfill(
        &mut effective,
        &mut updates,
        "mondo_disease_scope_registry",
        join(&root, "mondo/disease_scope.tsv.gz"),
    );
    Ok(updates)
}
pub(crate) fn portable_updates(ctx: &Context, config: &Config, root: &Path) -> Result<Paths> {
    let (_, mut updates) = plan::resolve_cache_paths(ctx, config, Some(root))?;
    for (k, v) in portable_defaults(&config.overlay(&updates), root)? {
        set(&mut updates, &k, v);
    }
    Ok(updates)
}
pub(crate) fn automatic_source_paths(
    ctx: &Context,
    config: &Config,
    updates: &Paths,
    root: &Path,
) -> Result<Vec<(String, Vec<PathBuf>)>> {
    if !matches!(
        config.get("assembly").string(),
        Some("hg19" | "GRCh37" | "hg38" | "GRCh38")
    ) {
        return Ok(vec![]);
    }
    let effective = config.overlay(updates);
    let mut default_config = effective.clone();
    let dirs = FOLDERS
        .iter()
        .map(|(_, _, k)| {
            Ok((
                k.to_string(),
                effective.index(k)?.path()?.display().to_string(),
            ))
        })
        .collect::<Result<Paths>>()?;
    for (key, _) in resource_groups(ctx, config, &dirs)? {
        if !generated_directory(&key) {
            default_config.remove(&key);
        }
    }
    let defaults = portable_defaults(&default_config, root)?;
    let mut paths: Vec<(String, Vec<PathBuf>)> = defaults
        .iter()
        .filter(|(k, _)| AUTO.contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), vec![path_value(v)]))
        .collect();
    let cache = checked_storage(config.index("cache_root")?, "cache_root")?;
    let old = matches!(config.get("assembly").string(), Some("hg19" | "GRCh37"));
    let ucsc = if old { "hg19" } else { "hg38" };
    paths.push((
        "ref_genome".into(),
        vec![join(&cache, format!("reference/ucsc.{ucsc}.fasta"))],
    ));
    let plugin = effective.index("vep_plugins_cachedir")?.path()?;
    let plugins = effective.index("vep_plugins_dir")?.path()?;
    let aliases = if old {
        ["hg19", "GRCh37"]
    } else {
        ["hg38", "GRCh38"]
    };
    let mut add = |key: &str, path: PathBuf| {
        if let Some((_, v)) = paths.iter_mut().find(|(k, _)| k == key) {
            if !v.contains(&path) {
                v.push(path);
            }
        } else {
            paths.push((key.into(), vec![path]));
        }
    };
    for (key, name) in [
        ("human_ancestor_fasta", "human_ancestor.fa.gz"),
        ("loftee_conservation_file", "phylocsf_gerp.sql"),
        ("gerp_bigwig", "gerp_conservation_scores.homo_sapiens.bw"),
    ] {
        if key != "gerp_bigwig" || !old {
            for alias in aliases {
                add(key, join(&plugin, format!("LoFtee/{alias}/{name}")));
            }
        }
    }
    for alias in aliases {
        add(
            "loftee_repo",
            join(&plugins, format!("LoFtee_repos/loftee-{alias}")),
        );
        if !old {
            add(
                "mavedb_file",
                join(&plugin, format!("MaveDB/MaveDB_variants.{alias}.tsv.gz")),
            );
        }
    }
    Ok(paths)
}
