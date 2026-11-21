//! The annotation's first step without `preprocess_vcf` of annotation_vcf.sh:
//! the records of the primary chromosomes (`bcftools view -r | bcftools
//! sort`), UCSC chromosome names (`bcftools annotate --rename-chrs` when the
//! VCF uses GRC names), biallelic, normalized records (`bcftools norm`
//! twice, `bcftools filter`, `bcftools sort`), then priva's variant IDs.
//! The intermediate files are named and removed as the function did; it
//! went on after a failed rename or normalization, and so does this, so the
//! failure shows at the variant IDs.
use super::annotation::{fasta_assembly, rename_chromosomes, vcf_assembly};
use super::tools::{
    Context, Result, args, log, option, refuse, replace_first, replace_vcf, unquoted_path, words,
};
use std::fs;
use std::path::Path;

/// `announce_remove_tmps FILE`: removes a file (a folder is not given here).
fn remove_announced(path: &str) {
    if Path::new(path).is_file() {
        match fs::remove_file(path) {
            Ok(()) => log(&format!("{path} has been succesfully removed")),
            Err(_) => log(&format!(
                "File {path} failed to be deleted (either its not existed or occupied)"
            )),
        }
    }
}

/// `announce_remove_tmps ${PREFIX}*tmp*vcf*`: the files of PREFIX's folder
/// whose name starts with PREFIX's name and holds `tmp`, then `vcf`, after it.
fn remove_temporary(prefix: &str) {
    let path = Path::new(prefix);
    let folder = path.parent().filter(|p| !p.as_os_str().is_empty());
    let start = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let Ok(entries) = fs::read_dir(folder.unwrap_or(Path::new("."))) else {
        return;
    };
    let mut names: Vec<String> = entries
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    for name in names {
        let Some(rest) = name.strip_prefix(&start) else {
            continue;
        };
        let matches = rest
            .find("tmp")
            .is_some_and(|at| rest[at + 3..].contains("vcf"));
        if matches && !name.starts_with('.') {
            let full = match folder {
                Some(folder) => folder.join(&name).display().to_string(),
                None => name.clone(),
            };
            remove_announced(&full);
        }
    }
}

/// `normalize_vcf IN OUT REF THREADS`: `bcftools norm -m -both -c s -f REF
/// --atom-overlaps . [--keep-sum AD] -a | bcftools norm -d exact |
/// bcftools filter -i EXPR | bcftools sort`, then `tabix`. AD is summed
/// when the header defines it and GT too; without GT only `*` alleles go.
fn normalize(
    context: &Context,
    input: &str,
    output: &str,
    reference: &str,
    threads: &str,
) -> Result<bool> {
    if !output.ends_with(".vcf.gz") {
        log(&format!(
            "Output vcf file {output} is not ending with .vcf.gz. Only support .vcf.gz format now."
        ));
        return Ok(false);
    }
    let header = context.header(input);
    let mut keep_sum = if header.contains("##FORMAT=<ID=AD,") {
        args(&["--keep-sum", "AD"])
    } else {
        Vec::new()
    };
    let filter = if header.contains("##FORMAT=<ID=GT,") {
        "ALT!=\"*\" && GT=\"alt\""
    } else {
        keep_sum.clear();
        "ALT != \"*\""
    };
    let mut first = vec!["norm".to_owned()];
    first.extend(option("--threads", threads));
    first.extend(args(&["-m", "-both", "-c", "s"]));
    first.extend(option("-f", reference));
    first.extend(args(&["--atom-overlaps", "."]));
    first.extend(keep_sum);
    first.extend(args(&["-Ou", "-a"]));
    first.extend(words(input));
    let mut second = vec!["norm".to_owned()];
    second.extend(option("--threads", threads));
    second.extend(args(&["-d", "exact", "-Ou", "-"]));
    let mut third = vec!["filter".to_owned()];
    third.extend(option("--threads", threads));
    third.extend(args(&["-i", filter, "-Ou", "-"]));
    let sorted = context.pipeline(
        &[
            ("bcftools", first),
            ("bcftools", second),
            ("bcftools", third),
            ("bcftools", args(&["sort", "-Oz", "-o", output])),
        ],
        false,
    )?;
    Ok(sorted && context.succeeds("tabix", &args(&["-f", "-p", "vcf", output])))
}

/// `preprocess_vcf -i VCF -o OUTPUT -r ref_genome -t threads`.
pub(crate) fn run(context: &Context, input: &str, output: &str) -> Result<()> {
    let reference = context.value("ref_genome");
    let configured_threads = context.value("threads");
    for (value, name) in [(input, "input_vcf"), (reference.as_str(), "ref_genome")] {
        if let Some(error) = unquoted_path(value, "file", name) {
            log(&format!("Error: {error}"));
            return refuse(error);
        }
    }
    let assembly = vcf_assembly(context, input);
    let threads = if configured_threads.is_empty() {
        "1"
    } else {
        configured_threads.as_str()
    };
    if reference.is_empty() {
        return refuse(
            "User does not specify the ref genome fasta file. Cant proceed now. Quit with Error!",
        );
    }
    if !assembly.is_empty() && fasta_assembly(&reference) != assembly {
        return refuse(format!(
            "User specified ref_genome fasta file {reference} seems not match with the VCF used assembly {assembly}. Quit with Error"
        ));
    }
    log(&format!("Expect the output vcf file at {output}"));

    // The records of the primary chromosomes (chrM to chrY, MT to Y).
    let primary = replace_vcf(input, ".primary.vcf.gz");
    let contigs = fs::read_to_string(context.data("liftover/ucsc_GRC.primary.contigs.tsv"))?
        .replace('\n', ",");
    let mut view = vec!["view".to_owned()];
    view.extend(option("--threads", threads));
    view.extend(args(&["-r", &contigs, "-Ou", input]));
    let kept = context.pipeline(
        &[
            ("bcftools", view),
            ("bcftools", args(&["sort", "-Oz", "-o", &primary])),
        ],
        false,
    )? && context.succeeds("tabix", &args(&["-f", "-p", "vcf", &primary]));
    if !kept {
        return refuse(
            "Failed to generate a VCF file that only contains records in primary chromosomes",
        );
    }

    // check_vcf_contig_version: UCSC names when the first record's starts
    // with chr, GRC names otherwise (also for a VCF without records).
    let ucsc = replace_vcf(input, ".ucsc.vcf.gz");
    let first = context.first_line("bcftools", &["query", "-f", "%CHROM\\n", &primary]);
    if first.starts_with("chr") {
        if fs::copy(&primary, &ucsc).is_err() {
            log(&format!("Failed to copy {primary} to {ucsc}"));
        }
    } else {
        log(
            "The input vcf is detected to map variants to GRC assemblies instead of UCSC assemblies",
        );
        let map = context.data("liftover/GRC_to_ucsc.contig.map.txt");
        if !Path::new(&map).is_file() {
            log(&format!(
                "The contig map file {map} is not found, please check the file"
            ));
        } else if !rename_chromosomes(context, &primary, &map, &ucsc)? {
            log(&format!("Failed to rename the chromosomes of {primary}"));
        }
    }

    let normalized = replace_vcf(input, ".norm.vcf.gz");
    if !normalize(context, &ucsc, &normalized, &reference, threads)? {
        log(&format!("Failed to normalize {ucsc}"));
    }
    remove_temporary(&replace_vcf(input, ""));

    // prepare_vcf_add_varID ${input_vcf/.vcf/.norm.vcf} OUTPUT.
    let ids_input = replace_first(input, ".vcf", ".norm.vcf");
    context.priva(
        "variant-ids",
        args(&["-m", "encode", "-i", &ids_input, "-o", output]),
    )?;
    if !context.succeeds("tabix", &args(&["-f", "-p", "vcf", output])) {
        return refuse(format!("Failed to index {output}"));
    }
    for path in [&ids_input, &ucsc, &primary] {
        remove_announced(path);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temporary_files_are_found_as_the_glob_found_them() {
        let folder = std::env::temp_dir().join(format!("priva_preprocess_{}", std::process::id()));
        let _ = fs::remove_dir_all(&folder);
        fs::create_dir_all(&folder).unwrap();
        for name in [
            "c.vcf.gz",
            "c.x.tmp.vcf.gz",
            "c.tmp.tsv",
            "d.tmp.vcf",
            "c.anno.vcf.gz",
        ] {
            fs::write(folder.join(name), "").unwrap();
        }
        remove_temporary(&folder.join("c").display().to_string());
        let mut left: Vec<String> = fs::read_dir(&folder)
            .unwrap()
            .flatten()
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        left.sort();
        assert_eq!(
            left,
            ["c.anno.vcf.gz", "c.tmp.tsv", "c.vcf.gz", "d.tmp.vcf"]
        );
        fs::remove_dir_all(&folder).unwrap();
    }
}
