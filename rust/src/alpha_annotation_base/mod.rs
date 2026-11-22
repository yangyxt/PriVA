//! The AlphaMissense annotation checks that other commands share: the
//! annotated VCF's header against the configured VEP release and assembly
//! (`header`), and the file and JSON helpers they use (`io`). The
//! `alpha-annotation-base` command, which only the shell installer's
//! AlphaMissense functions ran, went with that installer.
pub(crate) mod header;
pub(crate) mod io;

use crate::cli::{Result, invalid};
use crate::json_text::Json;

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
