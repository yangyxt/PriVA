//! HTTP source policy on the shared CPython urlsplit syntax checks.
use super::*;
use crate::python_text::is_python_space;
pub(crate) fn check(value: &str) -> Result<()> {
    let parsed = crate::python_url::split(value)?;
    if !matches!(parsed.scheme.as_str(), "http" | "https")
        || parsed.netloc.is_empty()
        || value.chars().any(is_python_space)
    {
        return Err(error(
            "ValueError",
            "Source URL must be an absolute HTTP(S) URL",
        ));
    }
    Ok(())
}
