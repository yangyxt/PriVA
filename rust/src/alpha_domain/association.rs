//! Metadata-cache publication identity. As in the Python helper, VCF identity
//! uses metadata, implementation identity uses content, and a changed output
//! invalidates reuse before its payload is deserialized.
use super::{
    Json, Options, Result, builder_identity, exception, file_identity, json_text, parquet,
    verify_input,
};
use std::io::Write;
use std::path::{Path, PathBuf};

fn sidecar(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".source.json");
    name.into()
}
fn regular(path: &Path) -> Result<()> {
    let value = std::fs::metadata(path)?;
    if !value.is_file() || value.len() == 0 {
        return Err(exception(
            "ValueError",
            "metadata-cache input/output must be a nonempty regular file",
        ));
    }
    Ok(())
}
pub(super) fn source(path: &Path) -> Result<()> {
    regular(path)
}

fn encoded(output: &Path, input: &str, options: &Options) -> Result<String> {
    regular(output)?;
    let value = Json::Object(vec![
        (
            "schema".into(),
            Json::Str("priva-domain-metadata-association-v1".into()),
        ),
        (
            "input".into(),
            json_text::parse(input, true).map_err(super::invalid)?,
        ),
        (
            "builder_source_sha256".into(),
            Json::Str(builder_identity()?),
        ),
        (
            "assembly".into(),
            options
                .assembly
                .clone()
                .map(Json::Str)
                .unwrap_or(Json::Null),
        ),
        ("output".into(), file_identity(output)?),
    ]);
    json_text::dumps(&value, true).map_err(super::invalid)
}
pub(super) fn matches(output: &Path, input: &str, options: &Options) -> bool {
    (|| -> Result<bool> {
        let expected = encoded(output, input, options)?;
        let text = std::fs::read_to_string(sidecar(output))?;
        let value = json_text::parse(&text, true).map_err(super::invalid)?;
        Ok(json_text::dumps(&value, true).map_err(super::invalid)? == expected)
    })()
    .unwrap_or(false)
}
pub(super) fn publish(output: &Path, input: &str, options: &Options) -> Result<()> {
    let receipt = sidecar(output);
    let source = json_text::parse(input, true).map_err(super::invalid)?;
    let source = source
        .get("vcf")
        .and_then(|v| v.get("path"))
        .and_then(Json::as_str)
        .ok_or_else(|| exception("ValueError", "metadata association has no VCF path"))?;
    super::protect(Path::new(source), &receipt)?;
    let expected = encoded(output, input, options)?;
    parquet::atomic_file(&receipt, |mut file, _| {
        file.write_all(expected.as_bytes())?;
        file.write_all(b"\n")?;
        file.sync_all()?;
        verify_input(input)?;
        if encoded(output, input, options)? != expected {
            return Err(exception(
                "ValueError",
                "metadata output changed before association",
            ));
        }
        Ok(())
    })
}
