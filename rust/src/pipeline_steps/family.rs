//! The family steps of the analysis without the shell scripts: filtration
//! (`main_filtration` of `filtration_vcf_per_fam.sh`) and the prioritization
//! steps (`prepare_combined_tab`, `interpret_splicing_annotations`,
//! `interpret_utr_annotations` and `acmg_assign_table` of
//! `prioritization_vcf_per_fam.sh`, as the task's step wrapper ran them).
//!
//! Each runs the same tools' commands with the same arguments: bcftools and
//! tabix, and priva's `pedigree-filter`, `combine-annotations`,
//! `splicing-analysis`, `utr-anno-interpret` and `acmg-assign`. Settings are
//! read from the task's private configuration as the scripts' `read_yaml`
//! read them (`yq '.KEY'` with the quotes taken off: a missing key or null is
//! `null`, a boolean `true` or `false`). A value the scripts passed unquoted
//! is split at white space as the shell split it, so an empty value gives no
//! word. The scripts' checks of the files they read are kept, and so are their
//! answers for an empty family: an empty table, or no table when the CADD
//! table has no row.
//!
//! Left out: the checks that skipped a step whose output was newer than its
//! inputs (a task starts in a new folder, so its output never exists before
//! it), and the scripts' table displays (log lines only).
use super::tools::{
    Context, Result, args, grep_matches, has_column, is_file, line_count, log, newer, not_empty,
    option, quoted_path, random_tag, refuse, replace_vcf, unquoted_path, words,
};
use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Stdio;

/// `check_path VALUE file NAME`, VALUE unquoted, collecting the message.
fn unquoted_file(value: &str, name: &str, errors: &mut Vec<String>) {
    errors.extend(unquoted_path(value, "file", name));
}

/// `check_path "VALUE" file NAME`, collecting the message.
fn quoted_file(value: &str, name: &str, errors: &mut Vec<String>) {
    errors.extend(quoted_path(value, "file", name));
}

// ---------------------------------------------------------------------------
// Prioritization

/// `combine`: the family table from the filtered VCF and CADD table.
pub(crate) fn combine(context: &Context, vcf: &str, cadd: &str) -> Result<()> {
    let table = replace_vcf(vcf, ".tsv");
    if context.count_records(vcf) < 1 && context.vcf_valid(vcf, 0, None) {
        fs::write(&table, "")?;
        return Ok(());
    }
    if line_count(cadd) <= 1 {
        log(&format!("The CADD table {cadd} is empty"));
        return Ok(());
    }
    if fs::File::open(context.config_path).is_err() {
        return refuse(format!(
            "The configuration file is not readable: {}",
            context.config_path.display()
        ));
    }
    let configured = context.value("hpo_assertions");
    let hpo = if configured.is_empty() || configured == "null" {
        context.data("hpo/genes_to_phenotype.assertions.tsv.gz")
    } else {
        configured
    };
    if !(not_empty(&hpo) && is_file(&hpo)) {
        return refuse(format!(
            "The selected HPO assertions file is missing or empty: {hpo}"
        ));
    }
    if !context.vcf_valid(vcf, 1, None) {
        return refuse(format!(
            "The input VCF file {vcf} is not valid, please check the VCF file"
        ));
    }
    // prepare_combined_tab quoted each value.
    let threads = context.value("threads_per_fam");
    context.priva(
        "combine-annotations",
        args(&[
            "-i", vcf, "-c", cadd, "-o", &table, "-p", &hpo, "-t", &threads,
        ]),
    )?;
    if !is_file(&table) {
        return refuse(format!("combine-annotations wrote no table {table}"));
    }
    Ok(())
}

/// The splicing and 5' UTR steps kept a table that already had their columns
/// and was newer than what they read.
fn already_done(table: &str, inputs: &[&str], binary: &str, columns: &[&str]) -> bool {
    is_file(table)
        && inputs.iter().all(|input| newer(table, input))
        && newer(table, binary)
        && columns.iter().all(|column| has_column(table, column))
}

/// `splicing`: priva splicing-analysis adds the splicing columns.
pub(crate) fn splicing(context: &Context, table: &str) -> Result<()> {
    if !not_empty(table) {
        return Ok(());
    }
    let domain_map = context.value("alphamissense_tranx_domain_map");
    let intolerant = context.value("all_intolerant_domains");
    let interpro = context.value("interpro_mapping_pickle");
    let columns = [
        "splicing_lof",
        "splicing_len_changing",
        "splicing_canonical_offset",
    ];
    if already_done(
        table,
        &[&interpro, &domain_map, &intolerant],
        context.binary,
        &columns,
    ) {
        log(&format!(
            "The input table {table} is up to date, skip the splicing interpretation"
        ));
        return Ok(());
    }
    let mut args = Vec::new();
    args.extend(option("--anno_table", table));
    args.extend(option("--transcript_domain_map", &domain_map));
    args.extend(option("--intolerant_domains", &intolerant));
    args.extend(option("--interpro_entry_map_pkl", &interpro));
    args.extend(option("--config", &context.config_text()));
    args.extend(option("--threads", &context.value("threads_per_fam")));
    context.priva("splicing-analysis", args)
}

/// `utr`: priva utr-anno-interpret adds the 5' UTR columns.
pub(crate) fn utr(context: &Context, table: &str) -> Result<()> {
    if !not_empty(table) {
        return Ok(());
    }
    let domain_map = context.value("alphamissense_tranx_domain_map");
    let intolerant = context.value("all_intolerant_domains");
    if already_done(
        table,
        &[&domain_map, &intolerant],
        context.binary,
        &["5UTR_lof", "5UTR_length_changing"],
    ) {
        log(&format!(
            "The input table {table} is up to date, skip the UTR interpretation"
        ));
        return Ok(());
    }
    let mut args = Vec::new();
    args.extend(option("--variants_table", table));
    args.extend(option("--domain_map", &domain_map));
    args.extend(option("--intolerant_domains", &intolerant));
    context.priva("utr-anno-interpret", args)
}

/// `acmg`: priva acmg-assign assigns the ACMG criteria and writes the
/// criteria matrix TABLE.acmg.tsv beside the table.
/// The review workbook of a family: the rows of the top-ranked variants of the
/// ACMG step's table as one Excel sheet (`crate::review_workbook`). The
/// setting `review_workbook_top_variants` says how many variants (1,000).
pub(crate) fn workbook(
    context: &Context,
    table: &str,
    family: Option<&str>,
    output: &std::path::Path,
) -> Result<()> {
    let setting = |key: &str| {
        let value = context.value(key);
        if value == "null" {
            String::new()
        } else {
            value
        }
    };
    let top_variants = match setting("review_workbook_top_variants").as_str() {
        "" => crate::review_workbook::DEFAULT_TOP_VARIANTS,
        count => count.parse().map_err(|_| {
            format!("review_workbook_top_variants: {count} is not a number of variants")
        })?,
    };
    let expected_incidence = match setting("exp_disease_incidence").as_str() {
        "" => 0.0001,
        value => value
            .parse()
            .map_err(|_| format!("exp_disease_incidence: {value} is not a number"))?,
    };
    let (ped, ranking) = (setting("ped_file"), setting("phenotype_gene_ranking"));
    let data = std::path::PathBuf::from(context.data(""));
    let summary = crate::review_workbook::write(&crate::review_workbook::Inputs {
        table: std::path::Path::new(table),
        output,
        top_variants,
        ped: is_file(&ped).then(|| std::path::Path::new(&ped)),
        family: family.filter(|name| !name.is_empty()),
        gene_ranking: is_file(&ranking).then(|| std::path::Path::new(&ranking)),
        expected_incidence,
        pext_tissues: &setting("pext_tissues"),
        data_dir: &data,
    })?;
    log(&format!(
        "Review workbook: {} rows of {} variants, {} columns",
        summary.rows, summary.variants, summary.columns
    ));
    Ok(())
}

pub(crate) fn acmg(context: &Context, table: &str, family: Option<&str>) -> Result<()> {
    if !not_empty(table) {
        // ${table%.tsv}.acmg.tsv
        fs::write(
            format!("{}.acmg.tsv", table.strip_suffix(".tsv").unwrap_or(table)),
            "",
        )?;
        return Ok(());
    }
    let v = |key: &str| context.value(key);
    let data = |relative: &str| context.data(relative);
    let mean_am_score_table = data("alphamissense/alphamissense_mean_score.tsv");
    let ped_table = v("ped_file");
    let clinvar_aa = v("clinvar_aa_stat");
    let clinvar_splice = v("clinvar_splice_stat");
    let clinvar_pm5 = v("clinvar_pm5_index");
    let intolerant = v("all_intolerant_domains");
    let am_intolerant = v("alphamissense_intolerant_domains");
    let mut clinvar_gene_stat = v("clinvar_gene_stat");
    if !clinvar_gene_stat.contains('/') {
        clinvar_gene_stat = data(&format!("ClinVar/{clinvar_gene_stat}"));
    }
    let domain_map = v("alphamissense_tranx_domain_map");
    let pm1_regions = v("pm1_regions_pkl");
    let alt_disease_vcf = v("alt_disease_vcf");
    let extreme_rare = v("extreme_rare_PAF");
    let mut incidence = v("exp_disease_incidence");
    if incidence.is_empty() || incidence == "null" {
        incidence = "0.0001".to_owned();
    }
    let patho_af = v("clinvar_patho_af_stat");
    let patho_exon_af = v("clinvar_patho_exon_af_stat");
    let interpro = v("interpro_mapping_pickle");
    let repeat_region = data(&format!("repeats/{}", v("repeat_region_file_name")));
    let pp1_vcf = v("pp1_vcf");
    let pp1_ped = v("pp1_ped");
    let clingen_map = v("clingen_map");
    let apply_clingen_override = v("apply_clingen_override");
    let dosage = v("gene_dosage_sensitivity");
    let relevant = v("relevant_gene_list");
    let dispensable = v("dispensable_gene_list");
    let cds_fasta = v("cds_fasta_file");
    let pext_tissues = v("pext_tissues");
    let condition_mechanism = v("hpo_condition_mechanism_json");
    let nonlof_mechanism = v("gene_nonlof_mechanism_json");

    let mut errors = Vec::new();
    for (value, name) in [
        (&clinvar_aa, "clinvar_aa_stat"),
        (&clinvar_splice, "clinvar_splice_stat"),
        (&clinvar_pm5, "clinvar_pm5_index"),
        (&intolerant, "all_intolerant_domains"),
        (&am_intolerant, "alphamissense_intolerant_domains"),
    ] {
        unquoted_file(value, name, &mut errors);
    }
    quoted_file(&clinvar_gene_stat, "clinvar_gene_stat", &mut errors);
    for (value, name) in [
        (&pm1_regions, "pm1_regions_pkl"),
        (&domain_map, "alphamissense_tranx_domain_map"),
        (&patho_af, "clinvar_patho_af_stat"),
        (&patho_exon_af, "clinvar_patho_exon_af_stat"),
        (&interpro, "interpro_mapping_pickle"),
        (&repeat_region, "repeat_region_file"),
        (&clingen_map, "clingen_map"),
        (&dosage, "gene_dosage_sensitivity"),
        (&condition_mechanism, "hpo_condition_mechanism_json"),
        (&nonlof_mechanism, "gene_nonlof_mechanism_json"),
    ] {
        unquoted_file(value, name, &mut errors);
    }
    let cds_fasta_arg = if is_file(&cds_fasta) {
        vec!["--cds_fasta_path".to_owned(), cds_fasta.clone()]
    } else {
        log(&format!(
            "WARNING: CDS FASTA file not found at {cds_fasta}. Alternative start codon detection for start_lost variants will use deprecated transcript-based approach."
        ));
        Vec::new()
    };
    let clingen_override_arg = if apply_clingen_override == "false" {
        log(
            "ClinGen override DISABLED by config. PriVA will report the criteria it derived itself.",
        );
        vec!["--ignore_clingen_override".to_owned()]
    } else {
        Vec::new()
    };
    let assembly = v("assembly");
    let mavedb_arg = if assembly == "hg38" || assembly == "GRCh38" {
        let metadata = data("MaveDB/assay_mapping_result.neutral_bands.tsv");
        unquoted_file(&metadata, "mavedb_assay_metadata", &mut errors);
        option("--mavedb_metadata_tsv", &metadata)
    } else {
        Vec::new()
    };
    if !errors.is_empty() {
        for error in &errors {
            log(&format!("Error: {error}"));
        }
        return refuse(format!(
            "Failed to offer the valid required files for the ACMG criteria assignment: {}",
            errors.join("; ")
        ));
    }

    let pp1_arg = if !pp1_vcf.is_empty()
        && words(&pp1_vcf)
            .first()
            .is_some_and(|path| context.vcf_valid(path, 1, None))
        && is_file(&pp1_ped)
    {
        let mut arg = option("--pp1_vcf", &pp1_vcf);
        arg.extend(option("--pp1_ped", &pp1_ped));
        arg
    } else {
        Vec::new()
    };
    let ped_arg = if is_file(&ped_table) {
        option("--ped_table", &ped_table)
    } else {
        Vec::new()
    };
    let fam_arg = match family.filter(|f| !f.is_empty()) {
        Some(name) => option("--fam_name", name),
        None => Vec::new(),
    };
    let alt_disease_arg = if !alt_disease_vcf.is_empty()
        && words(&alt_disease_vcf)
            .first()
            .is_some_and(|path| context.vcf_valid(path, 1, None))
    {
        option("--alt_disease_vcf", &alt_disease_vcf)
    } else {
        Vec::new()
    };
    let relevant_arg = if is_file(&relevant) {
        option("--relevant_gene_list", &relevant)
    } else {
        Vec::new()
    };
    let dispensable_arg = if is_file(&dispensable) {
        option("--dispensable_gene_list", &dispensable)
    } else {
        Vec::new()
    };
    let pext_arg = if !pext_tissues.is_empty() && pext_tissues != "null" {
        option("--pext_tissues", &pext_tissues)
    } else {
        Vec::new()
    };

    let mut args = Vec::new();
    for (name, value) in [
        ("--anno_table", table),
        ("--am_score_table", mean_am_score_table.as_str()),
        ("--clinvar_patho_af_stat", patho_af.as_str()),
        ("--clinvar_patho_exon_af_stat", patho_exon_af.as_str()),
        ("--clinvar_aa_dict_pkl", clinvar_aa.as_str()),
        ("--clinvar_splice_dict_pkl", clinvar_splice.as_str()),
        ("--clinvar_pm5_index", clinvar_pm5.as_str()),
        ("--interpro_entry_map_pkl", interpro.as_str()),
        ("--intolerant_domains_pkl", intolerant.as_str()),
        ("--gene_dosage_sensitivity", dosage.as_str()),
        (
            "--hpo_condition_mechanism_json",
            condition_mechanism.as_str(),
        ),
        ("--gene_nonlof_mechanism_json", nonlof_mechanism.as_str()),
        ("--am_intol_domains_tsv", am_intolerant.as_str()),
        ("--pm1_regions_pkl", pm1_regions.as_str()),
        ("--repeat_region_file", repeat_region.as_str()),
    ] {
        args.extend(option(name, value));
    }
    // The one value the script quoted.
    args.extend(["--clinvar_gene_stat_pkl".to_owned(), clinvar_gene_stat]);
    let config_text = context.config_path.display().to_string();
    let threads = v("threads_per_fam");
    for (name, value) in [
        ("--clingen_map_pkl", clingen_map.as_str()),
        ("--gnomAD_extreme_rare_threshold", extreme_rare.as_str()),
        ("--expected_incidence", incidence.as_str()),
        ("--config", config_text.as_str()),
        ("--threads", threads.as_str()),
        ("--tranx_exon_domain_map_pkl", domain_map.as_str()),
    ] {
        args.extend(option(name, value));
    }
    for extra in [
        ped_arg,
        fam_arg,
        alt_disease_arg,
        mavedb_arg,
        pp1_arg,
        relevant_arg,
        dispensable_arg,
        cds_fasta_arg,
        pext_arg,
        clingen_override_arg,
    ] {
        args.extend(extra);
    }
    context.priva("acmg-assign", args)
}

/// The whole prioritization of one family (`main_prioritization`): the
/// table, then splicing, 5' UTR and ACMG on that same table.
pub(crate) fn prioritize(
    context: &Context,
    vcf: &str,
    cadd: &str,
    family: Option<&str>,
) -> Result<()> {
    let table = replace_vcf(vcf, ".tsv");
    if context.count_records(vcf) < 1 && context.vcf_valid(vcf, 0, None) {
        log(&format!(
            "WARNING: The input VCF file {vcf} is empty, skip the prioritization with warnings"
        ));
        fs::write(&table, "")?;
        fs::write(format!("{}.acmg.tsv", &table[..table.len() - 4]), "")?;
        return Ok(());
    }
    if line_count(cadd) <= 1 {
        log(&format!(
            "WARNING: The CADD annotation file {cadd} is empty, just quit with warning"
        ));
        return Ok(());
    }
    combine(context, vcf, cadd)?;
    if !is_file(&table) {
        return refuse("Failed to prepare the combined annotation table");
    }
    splicing(context, &table)?;
    utr(context, &table)?;
    acmg(context, &table, family)
}

// ---------------------------------------------------------------------------
// Filtration

/// `main_filtration`: with a PED file, the family's samples
/// (`bcftools view -s`), the records priva pedigree-filter keeps, then the
/// allele-frequency filter on that file; without one, the allele-frequency
/// filter on the cohort VCF. Then the CADD rows of the records that are left.
pub(crate) fn filtration(context: &Context, vcf: &str, family: Option<&str>) -> Result<()> {
    let ped = context.value("ped_file");
    let threads = context.value("threads_per_fam");
    let family = family.unwrap_or("");
    let filtered = if is_file(&ped) {
        let family_vcf = extract_family(context, vcf, &ped, family, &threads)?;
        let output = replace_vcf(&family_vcf, ".filtered.vcf.gz");
        let mut args = Vec::new();
        args.extend(option("-v", &family_vcf));
        args.extend(option("-p", &ped));
        args.extend(option("-f", family));
        args.extend(option("-o", &output));
        context.priva("pedigree-filter", args)?;
        context.run(
            "tabix",
            &["-p".into(), "vcf".into(), "-f".into(), output.clone()],
        )?;
        filter_af(context, &output, None, &threads)?;
        output
    } else {
        log(
            "No pedigree table and family name provided, directly filter the variants based on the allele frequency to get the filtered VCF file",
        );
        let output = replace_vcf(vcf, ".filtered.vcf.gz");
        filter_af(context, vcf, Some(&output), &threads)?;
        output
    };
    if !is_file(&filtered) {
        log(&format!(
            "Skipping CADD filtering: Filtered VCF ({filtered}) not found"
        ));
        return Ok(());
    }
    let name = Path::new(&filtered)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let base = name.strip_suffix(".filtered.vcf.gz").unwrap_or(&name);
    let parent = Path::new(&filtered)
        .parent()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| ".".into());
    let family_cadd = format!("{parent}/{base}.filtered.cadd.tsv");
    let input_vcf = context.value("input_vcf");
    let input_name = Path::new(&replace_vcf(&input_vcf, ".anno.cadd.tsv"))
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let original_cadd = format!("{}/{input_name}", context.value("output_dir"));
    if !is_file(&original_cadd) {
        log(&format!(
            "Skipping CADD filtering: Original CADD file ({original_cadd}) not found"
        ));
        return Ok(());
    }
    filter_cadd(context, &filtered, &original_cadd, &family_cadd)
}

/// `extract_fam_vcf`: the family's samples from the PED file (column 2 of
/// the rows whose column 1 is the family), with `bcftools view --force-samples`.
fn extract_family(
    context: &Context,
    vcf: &str,
    ped: &str,
    family: &str,
    threads: &str,
) -> Result<String> {
    let threads = if threads.is_empty() { "1" } else { threads };
    let mut members = String::new();
    for line in fs::read_to_string(ped)?.lines() {
        let mut fields = line.split_whitespace();
        if fields.next() == Some(family) {
            if let Some(sample) = fields.next() {
                members.push_str(sample);
                members.push(',');
            }
        }
    }
    // ${fam_members::-1}: the list without its last comma; the shell refused
    // the expansion of an empty list.
    if members.pop().is_none() {
        return refuse(format!(
            "The PED file {ped} lists no member of family {family}"
        ));
    }
    if !context.vcf_valid(vcf, 1, None) {
        return refuse(format!("The input VCF {vcf} is not valid"));
    }
    let staged = replace_vcf(vcf, &format!(".{}.vcf.gz", random_tag()));
    let output = replace_vcf(vcf, &format!(".{family}.vcf.gz"));
    let mut args = vec!["view".to_owned()];
    args.extend(option("--threads", threads));
    args.push("--force-samples".into());
    args.extend(option("-s", &members));
    args.extend(["-Oz".into(), "-o".into(), staged.clone(), vcf.to_owned()]);
    context.run("bcftools", &args)?;
    fs::rename(&staged, &output)?;
    context.run(
        "tabix",
        &["-p".into(), "vcf".into(), "-f".into(), output.clone()],
    )?;
    Ok(output)
}

/// `filter_af`: the records whose AF_grpmax_joint is not above af_cutoff,
/// `bcftools view -e | bcftools sort`. Without OUTPUT the input is replaced.
/// A VCF whose header already records this filter is kept as it is.
fn filter_af(context: &Context, vcf: &str, output: Option<&str>, threads: &str) -> Result<()> {
    let threads = if threads.is_empty() { "1" } else { threads };
    if !context.config_path.is_file() {
        return refuse(format!(
            "The config file {} does not exist",
            context.config_path.display()
        ));
    }
    let cutoff = context.value("af_cutoff");
    let expression = format!("AF_grpmax_joint > {cutoff}");
    if context.vcf_valid(vcf, 1, None)
        && context
            .output("bcftools", &["view", "-h", vcf])
            .lines()
            .any(|line| grep_matches(line, &expression))
    {
        log(&format!(
            "The input VCF {vcf} already been filtered on allele frequency at the cutoff {cutoff}"
        ));
        return Ok(());
    }
    let target = match output {
        Some(path) => path.to_owned(),
        None => replace_vcf(vcf, &format!(".{}.vcf.gz", random_tag())),
    };
    let mut view_args = vec!["view".to_owned()];
    view_args.extend(option("--threads", threads));
    view_args.extend(["-e".into(), expression, "-Ou".into(), vcf.to_owned()]);
    log(&format!(
        "Running: bcftools {} | bcftools sort -Oz -o {target}",
        view_args.join(" ")
    ));
    let mut view = context
        .command("bcftools")
        .args(&view_args)
        .stdout(Stdio::piped())
        .spawn()?;
    let piped = view.stdout.take().expect("bcftools view output is piped");
    let sort = context
        .command("bcftools")
        .args(["sort", "-Oz", "-o", &target])
        .stdin(piped)
        .status()?;
    // The script's status was the sort's: no pipefail.
    let _ = view.wait();
    if !sort.success() {
        return refuse("Failed to filter the variants based on the allele frequency");
    }
    let indexed = match output {
        Some(_) => target,
        None => {
            fs::rename(&target, vcf)?;
            vcf.to_owned()
        }
    };
    context.run("tabix", &["-p".into(), "vcf".into(), "-f".into(), indexed])
}

/// `filter_cadd_based_on_vcf`: the CADD table's first line, then its rows
/// (from the second line on) whose chromosome (without `chr`), position,
/// reference and alternative allele are a record of the VCF.
fn filter_cadd(context: &Context, vcf: &str, cadd: &str, output: &str) -> Result<()> {
    let listing = context.output(
        "bcftools",
        &["query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", vcf],
    );
    let mut wanted: HashSet<String> = HashSet::new();
    for line in listing.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let chrom = get(0).strip_prefix("chr").unwrap_or(get(0));
        wanted.insert(format!("{chrom}\t{}\t{}\t{}", get(1), get(2), get(3)));
    }
    log(&format!("Loaded {} variants from VCF", wanted.len()));
    let staged = format!("{output}.{}.partial", random_tag());
    {
        let mut reader = BufReader::new(fs::File::open(cadd)?);
        let mut writer = BufWriter::new(fs::File::create(&staged)?);
        let mut line = Vec::new();
        let mut first = true;
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if first {
                // head -n 1
                first = false;
                writer.write_all(&line)?;
                continue;
            }
            let text = String::from_utf8_lossy(&line);
            let row = text.trim_end_matches('\n');
            let fields: Vec<&str> = row.split('\t').collect();
            let get = |i: usize| fields.get(i).copied().unwrap_or("");
            let chrom = get(0).strip_prefix("chr").unwrap_or(get(0));
            if wanted.contains(&format!("{chrom}\t{}\t{}\t{}", get(1), get(2), get(3))) {
                writer.write_all(row.as_bytes())?;
                writer.write_all(b"\n")?;
            }
        }
        writer.flush()?;
    }
    fs::rename(&staged, output)?;
    log(&format!(
        "CADD table filtering complete: {} of {} lines kept in {output}",
        line_count(output),
        line_count(cadd)
    ));
    Ok(())
}
