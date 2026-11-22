//! The annotated ClinVar VCF a configuration selects (`deployment_context`)
//! and the bounded check of its VEP header (`deployment_header`).
use super::*;
use crate::installer_prerequisites::paths::{abspath, is_relative_to, is_symlink};
use crate::python_path::absolute;

fn resolve(path: &str) -> Result<String> {
    crate::installer_prerequisites::paths::resolve(path)
}

/// `str(path).split(".vcf", 1)[0]`
fn first_vcf_stem(path: &str) -> &str {
    path.split_once(".vcf").map_or(path, |(stem, _)| stem)
}

/// The annotated ClinVar VCF a configuration selects, with the assembly and
/// VEP major its header must name.
pub(super) struct Context {
    target: String,
    assembly: &'static str,
    major: String,
}

/// `deployment_context(config_path)`: the selected annotation, refused when
/// it is a link, or when the step's own output lies outside cache_root.
pub(super) fn deployment_context(config_path: &str) -> Result<Context> {
    let config_path = absolute(Path::new(config_path))?;
    let config = load_config(&config_path)?;
    let assembly = assembly_of(&get_or_null(&config, "assembly")?)?;
    let major = vep_major(&config);
    let Some(assembly) = assembly.filter(|_| positive(&major)) else {
        return Err(value_error(
            "ClinVar deployment requires a configured assembly and numeric VEP major",
        ));
    };
    let folder = path_of(item(&config, "clinvar_vcf_dir")?)?;
    require(
        folder.starts_with('/'),
        "ClinVar deployment requires an absolute clinvar_vcf_dir",
    )?;
    let folder = abspath(&folder)?;
    let ucsc = ucsc(assembly);
    let target = abspath(&selected_vcf(&config, &folder, ucsc)?)?;
    let owned = target == join(&folder, &format!("clinvar.{ucsc}.vep.vcf.gz"));
    require(
        !is_symlink(&target)?,
        format!("Linked ClinVar annotation preserved without replacement: {target}"),
    )?;
    let root = match config.get("cache_root") {
        Some(Json::Str(text)) if text.starts_with('/') => Some(pathlib(text)),
        _ => None,
    };
    if owned {
        let inside = match &root {
            Some(root) => is_relative_to(&resolve(&folder)?, &resolve(root)?),
            None => false,
        };
        require(
            inside,
            "Owned ClinVar output must remain inside the configured cache_root; prior paths preserved",
        )?;
        require(
            format!("{}.vcf.gz", first_vcf_stem(&target)) == target,
            "Owned ClinVar output cannot have an earlier .vcf path component",
        )?;
    }
    Ok(Context {
        target,
        assembly,
        major,
    })
}

/// `deployment_header(context, input_file)`
pub(super) fn deployment_header(context: &Context, input: Option<&str>) -> Result<Json> {
    crate::alpha_annotation_base::header::check_resource(
        Path::new(input.unwrap_or(&context.target)),
        &string(context.assembly),
        &string(context.major.clone()),
        "ClinVar",
    )
}
