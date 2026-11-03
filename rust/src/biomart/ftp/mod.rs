//! Optional original FTP parser/projection procedures; no caller activation.
use super::*;
mod metadata;
pub(super) use metadata::{named_pin, publisher, verify_raw};
mod projection;
mod raw;
mod sqlite;
pub(super) fn sqlite_information() -> Json {
    sqlite::information()
}
type Projection = (String, String, Vec<String>);
fn projections() -> &'static [Projection] {
    static PROJECTIONS: OnceLock<Vec<Projection>> = OnceLock::new();
    PROJECTIONS.get_or_init(|| {
        let value = json_text::parse(include_str!("projections.json"), false)
            .expect("source projection constants");
        value
            .as_object()
            .unwrap()
            .iter()
            .map(|(alias, value)| {
                let pair = array(value);
                (
                    alias.clone(),
                    text(&pair[0]).to_owned(),
                    array(&pair[1]).iter().map(|v| text(v).to_owned()).collect(),
                )
            })
            .collect()
    })
}
pub(super) fn schema(path: &Path) -> Result<Json> {
    Ok(raw::schema_json(&raw::read_schema(path)?))
}
pub(super) fn load(
    db: &Path,
    raw: &Path,
    alias: &str,
    table: &str,
    columns: &Json,
    schema: &Json,
) -> Result<Json> {
    let db = sqlite::Connection::open(db)?;
    let columns = array(columns)
        .iter()
        .map(|v| {
            v.as_str()
                .map(str::to_owned)
                .ok_or_else(|| error("TypeError", "projection column must be a string"))
        })
        .collect::<Result<Vec<_>>>()?;
    let schema = array(schema)
        .iter()
        .map(|v| {
            let fields = array(v);
            if fields.len() != 2 {
                return Err(error(
                    "ValueError",
                    "projection schema entries require two values",
                ));
            }
            Ok((py_str(&fields[0]), py_str(&fields[1])))
        })
        .collect::<Result<Vec<_>>>()?;
    projection::load(&db, raw, alias, table, &columns, &schema)
}
pub(super) fn audit(db: &Path) -> Result<Json> {
    projection::audit(&sqlite::Connection::open(db)?)
}
pub(super) fn export(db: &Path, output: &Path, audit: &Json) -> Result<Json> {
    projection::export(&sqlite::Connection::open(db)?, output, audit)
}
/// The prepared-input projection block in FTP.install. Acquisition and final
/// scientific publication are separate procedures; this never selects inputs.
pub(super) fn project(raw: &Path, schema: &Path, db: &Path, output: &Path) -> Result<Json> {
    let schema = raw::read_schema(schema)?;
    let db = sqlite::Connection::open(db)?;
    db.execute("PRAGMA cache_size=-131072")?;
    let mut imports = Vec::new();
    for (alias, table, columns) in projections() {
        let fields = &schema.iter().find(|(name, _)| name == table).unwrap().1;
        imports.push((
            alias.clone(),
            projection::load(&db, raw, alias, table, columns, fields)?,
        ));
    }
    let audit = projection::audit(&db)?;
    let exported = projection::export(&db, output, &audit)?;
    require(
        db.scalar("PRAGMA integrity_check")? == sqlite::Value::Text("ok".into()),
        "FTP projection SQLite integrity failed",
    )?;
    Ok(object([
        ("imports", Json::Object(imports)),
        ("audit", audit),
        ("export", exported),
    ]))
}
