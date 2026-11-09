//! Ports of the former tests/nextflow/test_adapter.py, test_cpu_contract.py
//! and test_task_boundary.py. Synthetic fixtures only.
use super::contract::{self, Config, PRIORITY_KEYS};
use super::tasks;
use super::{ContractError, Options};
use crate::json_text::Json;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// Tests that change the working directory run one at a time.
static CWD: Mutex<()> = Mutex::new(());

struct Folder(PathBuf);

impl Folder {
    fn new(name: &str) -> Self {
        // Only characters the task contract accepts in a path.
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let serial = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("priva_task_{name}_{}_{serial}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Folder(fs::canonicalize(path).unwrap())
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn in_dir<T>(path: &Path, body: impl FnOnce() -> T) -> T {
    let original = std::env::current_dir().unwrap();
    std::env::set_current_dir(path).unwrap();
    let result = body();
    std::env::set_current_dir(original).unwrap();
    result
}

fn message(error: Box<dyn std::error::Error>) -> String {
    assert!(
        error.is::<ContractError>(),
        "not a contract refusal: {error}"
    );
    error.to_string()
}

fn samples(items: &[&str]) -> Vec<String> {
    items.iter().map(|s| s.to_string()).collect()
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

#[test]
fn source_digest_uses_nextflow_path_string_order() {
    let folder = Folder::new("digest");
    let source = folder.0.join("source-order");
    fs::create_dir_all(source.join("conf/cadd")).unwrap();
    for (name, content) in [
        ("conf/cadd.config", "flat\n"),
        ("conf/cadd/profile.json", "nested\n"),
        ("main.nf", "main\n"),
        ("nextflow.config", "settings\n"),
    ] {
        fs::write(source.join(name), content).unwrap();
    }
    let expected = crate::sha256::sha256_hex(
        b"conf/cadd.config\0flat\nconf/cadd/profile.json\0nested\nmain.nf\0main\nnextflow.config\0settings\n",
    );
    assert_eq!(contract::source_digest(&source).unwrap(), expected);
}

#[test]
fn ped_preserves_cohort_context_and_selects_family() {
    let folder = Folder::new("ped1");
    let path = folder.0.join("cohort.ped");
    let text = "OTHER\tother\t0\t0\t2\t2\nFAM\tchild\tdad\tmom\t1\t2\tHP:0000001\nFAM\tmom\t0\t0\t2\t1\nFAM\tdad\t0\t0\t1\t1\n";
    fs::write(&path, text).unwrap();
    let all = samples(&["child", "dad", "mom", "other"]);
    assert_eq!(
        tasks::validate_ped(&path, Some("FAM"), &all).unwrap(),
        vec![("FAM".to_owned(), samples(&["child", "mom", "dad"]))]
    );
    // Without a selection every family runs, in order of first appearance.
    assert_eq!(
        tasks::validate_ped(&path, None, &all).unwrap(),
        vec![
            ("OTHER".to_owned(), samples(&["other"])),
            ("FAM".to_owned(), samples(&["child", "mom", "dad"])),
        ]
    );
    assert_eq!(fs::read_to_string(&path).unwrap(), text);
}

#[test]
fn every_family_needs_its_affected_proband_first() {
    let folder = Folder::new("ped4");
    let path = folder.0.join("cohort.ped");
    fs::write(&path, "FAM\tchild\t0\t0\t1\t2\nOTHER\tother\t0\t0\t2\t1\n").unwrap();
    let all = samples(&["child", "other"]);
    // Selecting the valid family alone is accepted; the whole PED is not.
    assert!(tasks::validate_ped(&path, Some("FAM"), &all).is_ok());
    let error = tasks::validate_ped(&path, None, &all).unwrap_err();
    assert!(message(error).contains("affected proband first (family OTHER)"));
}

#[test]
fn ped_accepts_missing_parent_phenotypes_in_family_order() {
    let folder = Folder::new("ped2");
    let path = folder.0.join("cohort.ped");
    let all = samples(&["dad", "other", "child", "mom"]);
    for (mother, father) in [("-9", "-9"), ("0", "0"), ("-9", "0"), ("0", "-9")] {
        let text = format!(
            "OTHER\tother\t0\t0\t2\t2\nFAM\tchild\tdad\tmom\t0\t2\nFAM\tmom\t0\t0\t2\t{mother}\nFAM\tdad\t0\t0\t1\t{father}\n"
        );
        fs::write(&path, &text).unwrap();
        assert_eq!(
            tasks::validate_ped(&path, Some("FAM"), &all).unwrap(),
            vec![("FAM".to_owned(), samples(&["child", "mom", "dad"]))]
        );
    }
}

#[test]
fn ped_refusals() {
    let folder = Folder::new("ped3");
    let path = folder.0.join("family.ped");
    for (text, members, expected) in [
        ("FAM\tchild\t0\t0\t1\t2\n", vec![], "Every PED sample"),
        (
            "OTHER\tchild\t0\t0\t1\t2\n",
            vec!["child"],
            "absent from the PED",
        ),
        (
            "FAM\tchild\t0\t0\t1\t3\n",
            vec!["child"],
            "PED requires sex",
        ),
        (
            "FAM\tchild\t0\t0\t1\t1\n",
            vec!["child"],
            "affected proband first",
        ),
        ("FAM child 0 0 1 2\n", vec!["child"], "tab-delimited"),
        (
            "FAM\tmom\t0\t0\t2\t-9\nFAM\tchild\tdad\tmom\t0\t2\nFAM\tdad\t0\t0\t1\t-9\n",
            vec!["child", "mom", "dad"],
            "affected proband first",
        ),
        (
            "FAM\tchild\t0\tmom\t0\t2\nFAM\tmom\t0\t0\t3\t-9\n",
            vec!["child", "mom"],
            "PED requires sex",
        ),
        (
            "FAM\tchild\t0\tmom\t0\t2\nFAM\tmom\t0\t0\t2\t.\n",
            vec!["child", "mom"],
            "PED requires sex",
        ),
    ] {
        fs::write(&path, text).unwrap();
        let error = tasks::validate_ped(&path, Some("FAM"), &samples(&members)).unwrap_err();
        assert!(message(error).contains(expected), "{text:?}");
    }
}

#[test]
fn duplicate_yaml_key_is_refused() {
    let folder = Folder::new("yaml");
    let path = folder.0.join("bad.yaml");
    fs::write(&path, "assembly: hg19\nassembly: hg38\n").unwrap();
    assert!(message(contract::read_config(&path).err().unwrap()).contains("Duplicate YAML key"));
}

#[test]
fn spliceai_versions_are_read_from_headers_and_must_match() {
    let folder = Folder::new("splice");
    let path = folder.0.join("profile.yaml");
    let cache = folder.0.join("scores.any_extension");
    let header = b"##fileformat=VCFv4.2\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    for key in ["spliceai_snv_prescore", "spliceai_indel_prescore"] {
        for compressed in [false, true] {
            fs::write(
                &cache,
                if compressed {
                    gzip(header)
                } else {
                    header.to_vec()
                },
            )
            .unwrap();
            fs::write(
                &path,
                format!("assembly: hg38\n{key}: {}\n", cache.display()),
            )
            .unwrap();
            assert!(contract::read_config(&path).is_ok());
        }
    }
    let v1 = folder.0.join("v1.vcf");
    fs::write(
        &v1,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    fs::write(
        &path,
        format!(
            "assembly: hg38\nspliceai_snv_prescore: {}\nspliceai_indel_prescore: {}\n",
            cache.display(),
            v1.display()
        ),
    )
    .unwrap();
    let error = contract::read_config(&path).err().unwrap();
    assert!(message(error).contains("must be the same version"));
    fs::write(
        &cache,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    fs::write(
        &path,
        format!(
            "assembly: hg38\nspliceai_snv_prescore: {0}\nspliceai_indel_prescore: {0}\n",
            cache.display()
        ),
    )
    .unwrap();
    assert!(contract::read_config(&path).is_ok());
    fs::write(&cache, "not a VCF\n").unwrap();
    let error = contract::read_config(&path).err().unwrap();
    assert!(message(error).contains("Cannot verify spliceai_snv_prescore header"));
}

#[test]
fn missing_relative_and_lfs_paths_are_refused() {
    let folder = Folder::new("paths");
    let missing = folder.0.join("missing").display().to_string();
    for value in [
        "missing/cache",
        missing.as_str(),
        "",
        "/path with space/cache",
    ] {
        assert!(
            contract::checked_path(Some(value), "cache").is_err(),
            "{value}"
        );
    }
    let pointer = folder.0.join("pointer.tsv");
    fs::write(
        &pointer,
        "version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 123\n",
    )
    .unwrap();
    let error = contract::checked_path(Some(&pointer.display().to_string()), "cache").unwrap_err();
    assert!(message(error).contains("Git LFS pointer"));
}

#[test]
fn copy_is_private_even_when_original_is_read_only() {
    use std::os::unix::fs::MetadataExt;
    let folder = Folder::new("copy");
    let source = folder.0.join("source.tsv");
    fs::write(&source, "original\n").unwrap();
    fs::set_permissions(&source, fs::Permissions::from_mode(0o444)).unwrap();
    let copied = contract::copy_writable(&source, &folder.0.join("task/copy.tsv")).unwrap();
    fs::write(&copied, "changed\n").unwrap();
    assert_eq!(fs::read_to_string(&source).unwrap(), "original\n");
    assert_ne!(
        fs::metadata(&source).unwrap().ino(),
        fs::metadata(&copied).unwrap().ino()
    );
    assert_eq!(
        fs::metadata(&source).unwrap().mtime(),
        fs::metadata(&copied).unwrap().mtime()
    );
}

#[test]
fn inventory_tracks_nested_changes() {
    let folder = Folder::new("inventory");
    let root = folder.0.join("cache");
    fs::create_dir_all(root.join("release")).unwrap();
    let path = root.join("release/data.tsv");
    fs::write(&path, "before").unwrap();
    let paths = vec![("vep".to_owned(), root.clone())];
    let before = contract::inventory(&paths, false).unwrap();
    fs::write(&path, "after and different size").unwrap();
    let after = contract::inventory(&paths, false).unwrap();
    let key = "vep/release/data.tsv";
    assert!(before.get(key).is_some());
    assert!(!crate::json_text::py_equal(
        before.get(key).unwrap(),
        after.get(key).unwrap()
    ));
    // A link to a folder inside the resource is recorded, not followed.
    std::os::unix::fs::symlink(root.join("release"), root.join("inside")).unwrap();
    let records = contract::inventory(&paths, false).unwrap();
    let link = records.get("vep/inside").unwrap();
    assert!(matches!(link.get("type"), Some(Json::Str(kind)) if kind == "directory_link"));
    assert!(records.get("vep/inside/data.tsv").is_none());
    // A link leaving the resource is refused.
    let outside = folder.0.join("outside");
    fs::create_dir_all(&outside).unwrap();
    std::os::unix::fs::symlink(&outside, root.join("away")).unwrap();
    let error = contract::inventory(&paths, false).unwrap_err();
    assert!(message(error).contains("points outside"));
}

fn options(pairs: &[(&str, &str)], flags: &[&str]) -> Options {
    let mut args = Vec::new();
    for (key, value) in pairs {
        args.push(format!("--{key}"));
        args.push(value.to_string());
    }
    for flag in flags {
        args.push(format!("--{flag}"));
    }
    let names: Vec<&str> = pairs.iter().map(|(k, _)| *k).collect();
    Options::parse(&args, &names, flags).unwrap()
}

fn yaml_text(config: &Config, key: &str) -> String {
    config.text_of(key, "<absent>")
}

/// The private task configuration with the given extra settings.
fn private(
    folder: &Path,
    extra: &str,
    phase: &str,
    cpus: i64,
    mode: &str,
) -> Result<Config, String> {
    let bundle = folder.join("bundle");
    let task = folder.join("task");
    let _ = fs::remove_dir_all(&task);
    fs::create_dir_all(&bundle).unwrap();
    fs::create_dir_all(&task).unwrap();
    fs::write(bundle.join("family.ped"), "FAM\tchild\t0\t0\t0\t2\n").unwrap();
    let original = format!("assembly: hg38\n{extra}");
    fs::write(bundle.join("legacy.yaml"), &original).unwrap();
    let record = crate::json_text::parse(
        &format!(
            "{{\"base\": \"input\", \"assembly\": \"hg38\", \"mode\": \"{mode}\", \"vcf_hub_mode\": \"off\", \"cadd_hub_mode\": \"off\"}}"
        ),
        false,
    )
    .unwrap();
    let outcome = in_dir(&task, || {
        tasks::private_config_for_test(&bundle, &record, phase, cpus)
    });
    assert_eq!(
        fs::read_to_string(bundle.join("legacy.yaml")).unwrap(),
        original
    );
    match outcome {
        Ok(path) => Ok(contract::read_config(&path).unwrap()),
        Err(error) => {
            assert_eq!(
                fs::read_dir(&task).unwrap().count(),
                0,
                "task written before refusal"
            );
            Err(message(error))
        }
    }
}

#[test]
fn cpu_allocations_keep_the_scientific_worker_counts() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let folder = Folder::new("cpu");
    let counts = "threads: 22\nthreads_per_fam: 4\n";
    let config = private(&folder.0, counts, "annotation", 24, "raw").unwrap();
    assert_eq!(
        (
            yaml_text(&config, "threads"),
            yaml_text(&config, "threads_per_fam")
        ),
        ("22".into(), "4".into())
    );
    assert!(matches!(
        config.get("threads"),
        Some(crate::yaml_config::ConfigValue::Scalar(
            crate::yaml_config::Scalar::Int(_)
        ))
    ));
    let config = private(
        &folder.0,
        "threads: '22'\nthreads_per_fam: '4'\n",
        "annotation",
        24,
        "raw",
    )
    .unwrap();
    assert_eq!(config.str("threads"), Some("22"));
    assert_eq!(config.str("threads_per_fam"), Some("4"));
    let error = private(&folder.0, counts, "annotation", 23, "raw").unwrap_err();
    assert!(error.contains("at least 24 allocated CPUs"), "{error}");
    let config = private(&folder.0, counts, "prioritization", 8, "raw").unwrap();
    assert_eq!(yaml_text(&config, "threads_per_fam"), "4");
    let error = private(&folder.0, counts, "filtration", 3, "raw").unwrap_err();
    assert!(error.contains("at least 4 allocated CPUs"), "{error}");
    let config = private(&folder.0, counts, "annotation", 1, "annotated").unwrap();
    assert_eq!(yaml_text(&config, "threads"), "22");
    let error = private(&folder.0, "threads_per_fam: 4\n", "annotation", 24, "raw").unwrap_err();
    assert!(error.contains("explicit legacy-config threads"), "{error}");
    for value in [
        "true", "false", "0", "-1", "null", "''", "'2.0'", "auto", "2.0",
    ] {
        let error = private(
            &folder.0,
            &format!("threads: {value}\n"),
            "annotation",
            24,
            "raw",
        )
        .unwrap_err();
        assert!(
            error.contains("must be a positive integer"),
            "{value}: {error}"
        );
    }
    for cpus in [0, -1] {
        let error = private(&folder.0, "threads: 22\n", "annotation", cpus, "raw").unwrap_err();
        assert!(error.contains("Allocated cpus"), "{error}");
    }
    // The family phase fills an absent count from the allocation.
    let config = private(&folder.0, "threads: 22\n", "filtration", 3, "raw").unwrap();
    assert_eq!(yaml_text(&config, "threads_per_fam"), "3");
}

#[test]
fn private_config_keeps_scientific_values_and_isolates_mutable_resources() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let folder = Folder::new("private");
    let controls = folder.0.join("controls.vcf.gz");
    fs::write(&controls, b"immutable-control-bytes").unwrap();
    fs::write(format!("{}.tbi", controls.display()), b"index").unwrap();
    let cds = folder.0.join("coding.fa.gz");
    fs::write(&cds, b"fasta").unwrap();
    fs::write(folder.0.join("coding.cds_cache.pkl.gz"), b"pickle").unwrap();
    let extra = format!(
        "af_cutoff: 0.05\nextreme_rare_PAF: 0.0001\napply_clingen_override: no  # PyYAML false\ncontrol_vcf: {}\ncds_fasta_file: {}\nclinvar_aa_stat: /shared/same/cache.pkl\nhub_vcf_file: /shared/hub.vcf.gz\n",
        controls.display(),
        cds.display()
    );
    let config = private(&folder.0, &extra, "prioritization", 2, "raw").unwrap();
    let task = folder.0.join("task");
    assert_eq!(yaml_text(&config, "af_cutoff"), "0.05");
    assert!(
        config
            .text
            .contains("apply_clingen_override: false  # PyYAML false")
    );
    assert_eq!(
        config.str("clinvar_aa_stat"),
        Some("/shared/same/cache.pkl")
    );
    assert_eq!(yaml_text(&config, "threads_per_fam"), "2");
    assert_eq!(config.str("hub_vcf_file"), Some(""));
    assert_eq!(config.str("assembly"), Some("hg38"));
    let control = config.str("control_vcf").unwrap().to_owned();
    assert!(control.starts_with(&task.display().to_string()));
    let copied_cds = config
        .str("cds_fasta_file")
        .unwrap()
        .replace(".fa.gz", ".cds_cache.pkl.gz");
    assert!(Path::new(&copied_cds).exists());
    fs::write(&control, b"modified").unwrap();
    assert_eq!(fs::read(&controls).unwrap(), b"immutable-control-bytes");
}

struct LaunchFixture {
    folder: Folder,
    args: Vec<(String, String)>,
    resource: PathBuf,
}

impl LaunchFixture {
    fn new(name: &str) -> Self {
        let folder = Folder::new(name);
        let source = folder.0.join("source");
        fs::create_dir_all(&source).unwrap();
        for name in ["main.nf", "nextflow.config"] {
            fs::write(source.join(name), "frozen source\n").unwrap();
        }
        for relative in [
            "data/hpo/genes_to_phenotype.assertions.tsv.gz",
            "data/alphamissense/alphamissense_mean_score.tsv",
            "data/ClinVar/gene.pkl",
            "data/repeats/repeats.bed.gz",
            "data/hgnc/non_alt_loci_set.tsv",
            "data/dispensable_genes/dispensable_gene_list.txt",
        ] {
            let path = source.join(relative);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"synthetic resource fixture").unwrap();
        }
        let resource = folder.0.join("resource.pkl");
        fs::write(&resource, b"same fixed resource").unwrap();
        let priva = folder.0.join("priva");
        fs::write(&priva, "#!/bin/sh\necho 'priva 0.1.0 (fixture)'\n").unwrap();
        fs::set_permissions(&priva, fs::Permissions::from_mode(0o755)).unwrap();
        let mut config = String::new();
        for key in PRIORITY_KEYS {
            config.push_str(&format!("{key}: {}\n", resource.display()));
        }
        config.push_str(&format!(
            "assembly: hg19\nclinvar_gene_stat: gene.pkl\nrepeat_region_file_name: repeats.bed.gz\npriva_binary: {}\n",
            priva.display()
        ));
        let config_path = folder.0.join("legacy.yaml");
        fs::write(&config_path, config).unwrap();
        let manifest = folder.0.join("resources.json");
        let s = |p: &Path| p.display().to_string();
        tasks::build_manifest(&options(
            &[
                ("config", &s(&config_path)),
                ("legacy-root", &s(&source)),
                ("mode", "annotated"),
                ("output", &s(&manifest)),
            ],
            &[],
        ))
        .unwrap();
        fs::write(folder.0.join("family.ped"), "FAM\tchild\t0\t0\t1\t2\n").unwrap();
        fs::write(folder.0.join("input.vcf.gz"), b"fixture only").unwrap();
        fs::write(folder.0.join("input.vcf.gz.tbi"), b"index").unwrap();
        let digest = contract::source_digest(&source).unwrap();
        let args = [
            ("config", s(&config_path)),
            ("legacy-root", s(&source)),
            ("mode", "annotated".into()),
            ("vcf-hub-mode", "off".into()),
            ("cadd-hub-mode", "copy".into()),
            ("resources", s(&manifest)),
            ("source-digest", digest),
            ("software-id", "frozen-test-environment".into()),
            ("family", "FAM".into()),
            ("assembly", "hg19".into()),
            ("vcf", s(&folder.0.join("input.vcf.gz"))),
            ("vcf-index", s(&folder.0.join("input.vcf.gz.tbi"))),
            ("ped", s(&folder.0.join("family.ped"))),
            ("annotation-cpus", "1".into()),
            ("family-cpus", "1".into()),
            ("outdir", s(&folder.0.join("published"))),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_owned(), v))
        .collect();
        LaunchFixture {
            folder,
            args,
            resource,
        }
    }

    fn launch(&self) -> super::Result<()> {
        let pairs: Vec<(&str, &str)> = self
            .args
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        tasks::launch(&options(&pairs, &[]))
    }

    fn set(&mut self, key: &str, value: &str) {
        self.args.iter_mut().find(|(k, _)| k == key).unwrap().1 = value.to_owned();
    }

    fn value(&self, key: &str) -> String {
        self.args.iter().find(|(k, _)| k == key).unwrap().1.clone()
    }
}

#[test]
fn launch_rechecks_resources_even_for_identical_resume() {
    let fixture = LaunchFixture::new("launch1");
    fixture.launch().unwrap();
    fixture.launch().unwrap();
    fs::write(&fixture.resource, b"changed after completed cached run").unwrap();
    let refusal = message(fixture.launch().unwrap_err());
    assert!(refusal.contains("Resource changed"));
    // The refusal names the changed fields with both values.
    assert!(refusal.contains("size 19 -> 34"), "{refusal}");
}

#[test]
fn folder_size_is_not_compared_but_other_fields_are() {
    let record = |kind: &str, size: i128, mtime: i128| {
        contract::object(vec![
            ("path", contract::text("/r")),
            ("resolved_path", contract::text("/r")),
            ("type", contract::text(kind)),
            ("size", Json::Int(size)),
            ("mtime_ns", Json::Int(mtime)),
        ])
    };
    let change = |a: &Json, b: &Json| contract::record_change(a, b).unwrap();
    assert_eq!(
        change(&record("directory", 4096, 1), &record("directory", 0, 1)),
        None
    );
    assert_eq!(
        change(&record("directory", 0, 2), &record("directory", 0, 1)).as_deref(),
        Some(" (mtime_ns 1 -> 2)")
    );
    assert_eq!(
        change(&record("file", 5, 1), &record("file", 4, 1)).as_deref(),
        Some(" (size 4 -> 5)")
    );
}

#[test]
fn launch_rejects_changed_family_and_changed_source() {
    let mut fixture = LaunchFixture::new("launch2");
    fixture.launch().unwrap();
    fixture.set("family", "SECOND_FAMILY");
    assert!(message(fixture.launch().unwrap_err()).contains("different inputs"));
    let source = PathBuf::from(fixture.value("legacy-root"));
    fs::write(source.join("main.nf"), "changed source").unwrap();
    assert!(message(fixture.launch().unwrap_err()).contains("Source changed"));
}

#[test]
fn mixed_spliceai_versions_rejected_before_output_creation() {
    let fixture = LaunchFixture::new("launch3");
    let config = PathBuf::from(fixture.value("config"));
    let cache = fixture.folder.0.join("v2.vcf");
    fs::write(
        &cache,
        "##fileformat=VCFv4.2\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    let v1 = fixture.folder.0.join("v1.vcf");
    fs::write(
        &v1,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    let text = fs::read_to_string(&config).unwrap();
    fs::write(
        &config,
        format!(
            "{text}spliceai_snv_prescore: {}\nspliceai_indel_prescore: {}\n",
            cache.display(),
            v1.display()
        ),
    )
    .unwrap();
    assert!(message(fixture.launch().unwrap_err()).contains("must be the same version"));
    assert!(!Path::new(&fixture.value("outdir")).exists());
}

#[test]
fn mixed_spliceai_versions_task_rejected_before_any_output() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let fixture = LaunchFixture::new("launch4");
    let bundle = fixture.folder.0.join("v2-bundle");
    fs::create_dir_all(&bundle).unwrap();
    let cache = fixture.folder.0.join("v2.vcf");
    fs::write(
        &cache,
        "##fileformat=VCFv4.2\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    let v1 = fixture.folder.0.join("v1.vcf");
    fs::write(
        &v1,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    let config = fs::read_to_string(fixture.value("config")).unwrap();
    fs::write(
        bundle.join("legacy.yaml"),
        format!(
            "{config}spliceai_snv_prescore: {}\nspliceai_indel_prescore: {}\n",
            v1.display(),
            cache.display()
        ),
    )
    .unwrap();
    fs::copy(fixture.value("resources"), bundle.join("resources.json")).unwrap();
    let record = format!(
        "{{\"mode\": \"annotated\", \"stub\": false, \"legacy_root\": \"{}\", \"source_digest\": \"{}\", \"vcf_hub_mode\": \"off\", \"cadd_hub_mode\": \"copy\"}}",
        fixture.value("legacy-root"),
        fixture.value("source-digest")
    );
    fs::write(bundle.join("contract.json"), record).unwrap();
    for phase in ["annotation", "filtration", "prioritization"] {
        let bundle_text = bundle.display().to_string();
        let outcome = in_dir(&fixture.folder.0, || {
            tasks::run(&options(
                &[
                    ("phase", phase),
                    ("bundle", &bundle_text),
                    ("vcf", "input.vcf.gz"),
                    ("vcf-index", "input.vcf.gz.tbi"),
                    ("cpus", "1"),
                ],
                &[],
            ))
        });
        assert!(message(outcome.unwrap_err()).contains("must be the same version"));
        assert!(!fixture.folder.0.join("result").exists());
    }
}

#[test]
fn family_tasks_name_a_validated_family() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let fixture = LaunchFixture::new("family_run");
    let bundle = fixture.folder.0.join("family-bundle");
    fs::create_dir_all(&bundle).unwrap();
    let record = format!(
        "{{\"mode\": \"annotated\", \"stub\": true, \"base\": \"input\", \"families\": [\"FAM\", \"OTHER\"], \"legacy_root\": \"{}\", \"source_digest\": \"{}\", \"vcf_hub_mode\": \"off\", \"cadd_hub_mode\": \"copy\"}}",
        fixture.value("legacy-root"),
        fixture.value("source-digest")
    );
    fs::write(bundle.join("contract.json"), record).unwrap();
    let bundle_text = bundle.display().to_string();
    let run = |phase: &str, family: Option<&str>| {
        let mut pairs = vec![
            ("phase", phase),
            ("bundle", bundle_text.as_str()),
            ("vcf", "input.vcf.gz"),
            ("vcf-index", "input.vcf.gz.tbi"),
            ("cpus", "1"),
        ];
        if let Some(family) = family {
            pairs.push(("family", family));
        }
        in_dir(&fixture.folder.0, || tasks::run(&options(&pairs, &["stub"])))
    };
    let refusal = |phase, family| message(run(phase, family).unwrap_err());
    assert!(refusal("filtration", None).contains("name one with --family"));
    assert!(refusal("prioritization", Some("NOPE")).contains("not in this run's validated families"));
    assert!(refusal("annotation", Some("FAM")).contains("takes no --family"));
    // Refusals come before the task writes anything.
    assert!(!fixture.folder.0.join("result").exists());
}

#[test]
fn manifest_bundle_identity_is_checked() {
    let fixture = LaunchFixture::new("manifest");
    let path = PathBuf::from(fixture.value("resources"));
    let text = fs::read_to_string(&path).unwrap();
    let parsed = crate::json_text::parse(&text, false).unwrap();
    assert!(matches!(parsed.get("bundle_id"), Some(Json::Str(_))));
    assert!(text.ends_with("}\n") && text.starts_with("{\n  \"assembly\": \"hg19\""));
    fs::write(
        &path,
        text.replace("\"mode\": \"annotated\"", "\"mode\": \"raw\""),
    )
    .unwrap();
    assert!(message(fixture.launch().unwrap_err()).contains("bundle_id does not match"));
}

#[test]
fn raw_mode_checks_the_annotation_resources() {
    let fixture = LaunchFixture::new("gate");
    let mut raw = fixture.args.clone();
    raw.iter_mut().find(|(k, _)| k == "mode").unwrap().1 = "raw".into();
    let pairs: Vec<(&str, &str)> = raw.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
    // Raw mode needs explicit VEP workers and every annotation resource.
    let error = tasks::launch(&options(&pairs, &[])).unwrap_err();
    assert!(message(error).contains("explicit legacy-config threads"));
}

#[test]
fn boundary_checks_task_controls_and_aliases() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let folder = Folder::new("boundary");
    let task = folder.0.join("task");
    fs::create_dir_all(&task).unwrap();
    for name in [".command.sh", ".command.run"] {
        fs::write(task.join(name), "true\n").unwrap();
    }
    let check = |dir: &Path, cwd: &Path| {
        let args = vec![
            "check".to_owned(),
            "--task-dir".into(),
            dir.display().to_string(),
        ];
        in_dir(cwd, || super::boundary::main(&args))
    };
    assert_eq!(check(&task, &task), 0);
    assert_eq!(check(&task, &folder.0), 125);
    let link = folder.0.join("link");
    std::os::unix::fs::symlink(&task, &link).unwrap();
    assert_eq!(check(&link, &task), 125);
    let outside = folder.0.join("outside");
    fs::create_dir_all(&outside).unwrap();
    std::os::unix::fs::symlink(&outside, task.join("input")).unwrap();
    assert_eq!(check(&task, &task), 0);
    let cache = folder.0.join("cache");
    fs::write(&cache, "unchanged\n").unwrap();
    fs::hard_link(&cache, task.join("aliased")).unwrap();
    assert_eq!(check(&task, &task), 125);
    fs::remove_file(task.join("aliased")).unwrap();
    fs::remove_file(task.join(".command.run")).unwrap();
    assert_eq!(check(&task, &task), 125);
    std::os::unix::fs::symlink(task.join(".command.sh"), task.join(".command.run")).unwrap();
    assert_eq!(check(&task, &task), 125);
}

#[test]
fn cds_cache_copies_follow_the_detector_names() {
    let _guard = CWD.lock().unwrap_or_else(|e| e.into_inner());
    let folder = Folder::new("cds");
    // Names whose cache sits beside the FASTA, as the Python detector names it.
    for relative in [
        "standard.fa.gz",
        "plain.fa",
        "generic.other.gz",
        "refs.fa.gz/selected.fa",
        "refs.fa.gz/selected.other.gz",
        "refs.fa.gz/standard.fa.gz",
    ] {
        let fasta = folder.0.join("inputs").join(relative);
        fs::create_dir_all(fasta.parent().unwrap()).unwrap();
        fs::write(&fasta, b"fasta bytes").unwrap();
        let cache = contract::cds_cache_path(&fasta.display().to_string());
        fs::write(&cache, b"parsed cache").unwrap();
        let config = private(
            &folder.0,
            &format!("cds_fasta_file: {}\n", fasta.display()),
            "prioritization",
            1,
            "annotated",
        )
        .unwrap();
        let copied = config.str("cds_fasta_file").unwrap().to_owned();
        assert!(copied.starts_with(&folder.0.join("task").display().to_string()));
        assert_eq!(fs::read(&copied).unwrap(), b"fasta bytes");
        assert_eq!(
            fs::read(contract::cds_cache_path(&copied)).unwrap(),
            b"parsed cache",
            "{relative}"
        );
        fs::remove_file(&cache).unwrap();
    }
    // An old whole-path sidecar is never copied over the private FASTA.
    let fasta = folder.0.join("inputs/refs.fa.gz/selected.fa");
    let old = PathBuf::from(
        fasta
            .display()
            .to_string()
            .replace(".fa.gz", ".cds_cache.pkl.gz"),
    );
    fs::create_dir_all(old.parent().unwrap()).unwrap();
    fs::write(&old, b"old misplaced pickle").unwrap();
    let config = private(
        &folder.0,
        &format!("cds_fasta_file: {}\n", fasta.display()),
        "prioritization",
        1,
        "annotated",
    )
    .unwrap();
    let copied = config.str("cds_fasta_file").unwrap().to_owned();
    assert_eq!(fs::read(&copied).unwrap(), b"fasta bytes");
    assert!(!contract::cds_cache_path(&copied).exists());
    assert!(old.exists());
    // Other phases leave the configured FASTA in place.
    let config = private(
        &folder.0,
        &format!("cds_fasta_file: {}\n", fasta.display()),
        "filtration",
        1,
        "annotated",
    )
    .unwrap();
    assert_eq!(
        config.str("cds_fasta_file"),
        Some(fasta.display().to_string().as_str())
    );
}

#[test]
fn the_pipeline_program_is_required_and_content_pinned() {
    let fixture = LaunchFixture::new("program");
    let manifest = crate::json_text::parse(
        &fs::read_to_string(fixture.value("resources")).unwrap(),
        false,
    )
    .unwrap();
    let record = manifest
        .get("resources")
        .unwrap()
        .get("priva_binary")
        .unwrap();
    let program = fixture.folder.0.join("priva");
    assert!(matches!(record.get("sha256"), Some(Json::Str(digest))
        if *digest == crate::sha256::sha256_file(&program).unwrap()));
    assert!(matches!(
        record.get("executable_mode"),
        Some(Json::Int(0o755))
    ));
    // The same size and time but other bytes: the launch is refused.
    let original = fs::metadata(&program).unwrap();
    fs::write(&program, "#!/bin/sh\necho 'priva 0.1.1 (fixture)'\n").unwrap();
    let times = fs::FileTimes::new().set_modified(original.modified().unwrap());
    fs::OpenOptions::new()
        .write(true)
        .open(&program)
        .unwrap()
        .set_times(times)
        .unwrap();
    assert!(message(fixture.launch().unwrap_err()).contains("Resource changed"));
    // A configuration without the program is refused.
    let config = PathBuf::from(fixture.value("config"));
    let text = fs::read_to_string(&config).unwrap();
    let without: String = text
        .lines()
        .filter(|line| !line.starts_with("priva_binary:"))
        .map(|line| format!("{line}\n"))
        .collect();
    fs::write(&config, without).unwrap();
    let error = contract::runtime_resources(
        &contract::read_config(&config).unwrap(),
        Path::new(&fixture.value("legacy-root")),
        "annotated",
        "off",
        "copy",
    )
    .unwrap_err();
    assert!(message(error).contains("priva_binary"));
}
