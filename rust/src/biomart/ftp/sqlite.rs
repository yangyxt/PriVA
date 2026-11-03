//! Safe ownership around SQLite's stable C API, with CPython sqlite3's default
//! transaction and scalar-conversion boundaries. SQL remains in projection.rs.
//! Error classification follows CPython 3.12.3 Modules/_sqlite/util.c; this is
//! an independent narrow Rust adapter, not a replacement query engine.
use super::super::*;
use std::ffi::{CStr, CString, c_char, c_int, c_void};
use std::marker::PhantomData;
use std::os::unix::ffi::OsStrExt;
use std::ptr;

type Handle = *mut c_void;
unsafe extern "C" {
    fn sqlite3_open_v2(
        path: *const c_char,
        db: *mut Handle,
        flags: c_int,
        vfs: *const c_char,
    ) -> c_int;
    fn sqlite3_close_v2(db: Handle) -> c_int;
    fn sqlite3_busy_timeout(db: Handle, milliseconds: c_int) -> c_int;
    fn sqlite3_errcode(db: Handle) -> c_int;
    fn sqlite3_errmsg(db: Handle) -> *const c_char;
    fn sqlite3_prepare_v2(
        db: Handle,
        sql: *const c_char,
        bytes: c_int,
        stmt: *mut Handle,
        tail: *mut *const c_char,
    ) -> c_int;
    fn sqlite3_step(stmt: Handle) -> c_int;
    fn sqlite3_finalize(stmt: Handle) -> c_int;
    fn sqlite3_reset(stmt: Handle) -> c_int;
    fn sqlite3_clear_bindings(stmt: Handle) -> c_int;
    fn sqlite3_bind_null(stmt: Handle, index: c_int) -> c_int;
    fn sqlite3_bind_int64(stmt: Handle, index: c_int, value: i64) -> c_int;
    fn sqlite3_bind_text(
        stmt: Handle,
        index: c_int,
        value: *const c_char,
        bytes: c_int,
        destructor: Option<unsafe extern "C" fn(*mut c_void)>,
    ) -> c_int;
    fn sqlite3_bind_parameter_count(stmt: Handle) -> c_int;
    fn sqlite3_column_count(stmt: Handle) -> c_int;
    fn sqlite3_column_type(stmt: Handle, index: c_int) -> c_int;
    fn sqlite3_column_int64(stmt: Handle, index: c_int) -> i64;
    fn sqlite3_column_double(stmt: Handle, index: c_int) -> f64;
    fn sqlite3_column_text(stmt: Handle, index: c_int) -> *const u8;
    fn sqlite3_column_blob(stmt: Handle, index: c_int) -> *const c_void;
    fn sqlite3_column_name(stmt: Handle, index: c_int) -> *const c_char;
    fn sqlite3_column_bytes(stmt: Handle, index: c_int) -> c_int;
    fn sqlite3_get_autocommit(db: Handle) -> c_int;
    fn sqlite3_libversion() -> *const c_char;
    fn sqlite3_sourceid() -> *const c_char;
    fn sqlite3_compileoption_get(index: c_int) -> *const c_char;
}

#[derive(Debug, Clone, PartialEq)]
pub(super) enum Value {
    Null,
    Int(i64),
    Float(f64),
    Text(String),
    Blob(Vec<u8>),
    // Python first constructs arbitrary-size integer parameters. Overflow is
    // raised by executemany binding, not while reading a still-unflushed batch.
    OversizeInteger,
}
impl Value {
    pub(super) fn json(&self) -> Json {
        match self {
            Self::Null => Json::Null,
            Self::Int(n) => Json::Int(*n as i128),
            Self::Float(n) => Json::Float(*n),
            Self::Text(s) => string(s),
            Self::Blob(bytes) => Json::Array(bytes.iter().map(|b| Json::Int(*b as i128)).collect()),
            Self::OversizeInteger => unreachable!("only a parameter marker"),
        }
    }
}

pub(super) struct Connection(Handle);
pub(super) struct Statement<'a> {
    db: &'a Connection,
    handle: Handle,
    done: bool,
    // Explicit lifetime ownership; a statement is always finalized before DB.
    _owner: PhantomData<&'a Connection>,
}
unsafe fn c_text(value: *const c_char) -> String {
    // SAFETY: SQLite supplies a valid NUL-terminated diagnostic/static string.
    unsafe { CStr::from_ptr(value) }
        .to_string_lossy()
        .into_owned()
}
fn strip_sql(mut text: &str) -> &str {
    loop {
        text = text.trim_start_matches([' ', '\t', '\x0c', '\n', '\r']);
        if let Some(tail) = text.strip_prefix("--") {
            text = tail.find('\n').map_or("", |at| &tail[at + 1..]);
        } else if let Some(tail) = text.strip_prefix("/*") {
            text = tail.find("*/").map_or("", |at| &tail[at + 2..]);
        } else {
            return text;
        }
    }
}
impl Connection {
    pub(super) fn open(path: &Path) -> Result<Self> {
        let path = CString::new(path.as_os_str().as_bytes())
            .map_err(|_| error("ValueError", "embedded null byte"))?;
        let mut handle = ptr::null_mut();
        // SAFETY: valid string/out pointer; READWRITE|CREATE matches connect().
        let code = unsafe { sqlite3_open_v2(path.as_ptr(), &mut handle, 2 | 4, ptr::null()) };
        if handle.is_null() {
            return Err(error("MemoryError", ""));
        }
        let db = Self(handle);
        db.check(code)?;
        // Original sqlite3.connect defaults to a five-second busy timeout.
        db.check(unsafe { sqlite3_busy_timeout(db.0, 5000) })?;
        Ok(db)
    }
    fn check(&self, code: c_int) -> Result<()> {
        if matches!(code, 0 | 100 | 101) {
            return Ok(());
        }
        // SAFETY: connection owns a live handle; error read follows failing call.
        let (kind, message) = unsafe { (sqlite3_errcode(self.0), c_text(sqlite3_errmsg(self.0))) };
        let category = match kind {
            2 | 12 => "InternalError",
            7 => "MemoryError",
            1 | 3 | 4 | 5 | 6 | 8 | 9 | 10 | 13 | 14 | 15 | 16 | 17 => "OperationalError",
            18 => "DataError",
            19 | 20 => "IntegrityError",
            21 | 25 => "InterfaceError",
            _ => "DatabaseError",
        };
        Err(error(category, message))
    }
    pub(super) fn prepare(&self, sql: &str) -> Result<Statement<'_>> {
        let bytes = CString::new(sql)
            .map_err(|_| error("ProgrammingError", "the query contains a null character"))?;
        let size = c_int::try_from(bytes.as_bytes_with_nul().len())
            .map_err(|_| error("DataError", "query string is too large"))?;
        let mut handle = ptr::null_mut();
        let mut tail = ptr::null();
        // SAFETY: SQL remains alive during preparation; SQLite compiles a copy.
        self.check(unsafe {
            sqlite3_prepare_v2(self.0, bytes.as_ptr(), size, &mut handle, &mut tail)
        })?;
        let statement = Statement {
            db: self,
            handle,
            done: handle.is_null(),
            _owner: PhantomData,
        };
        // SAFETY: tail points inside this still-live SQL buffer.
        if !strip_sql(unsafe { CStr::from_ptr(tail) }.to_str().expect("UTF-8 SQL")).is_empty() {
            return Err(error(
                "ProgrammingError",
                "You can only execute one statement at a time.",
            ));
        }
        Ok(statement)
    }
    pub(super) fn execute(&self, sql: &str) -> Result<()> {
        let mut statement = self.prepare(sql)?;
        while statement.next()?.is_some() {}
        Ok(())
    }
    pub(super) fn scalar(&self, sql: &str) -> Result<Value> {
        let mut statement = self.prepare(sql)?;
        statement
            .next()?
            .and_then(|row| row.into_iter().next())
            .ok_or_else(|| error("TypeError", "'NoneType' object is not subscriptable"))
    }
    pub(super) fn executemany(&self, sql: &str, rows: &[Vec<Value>]) -> Result<()> {
        // Prepare errors precede BEGIN or binding failures in the original.
        let mut statement = self.prepare(sql)?;
        if unsafe { sqlite3_get_autocommit(self.0) } != 0 {
            self.execute("BEGIN ")?;
        }
        for row in rows {
            statement.bind(row)?;
            while statement.next()?.is_some() {}
            self.check(unsafe { sqlite3_reset(statement.handle) })?;
            statement.done = false;
        }
        Ok(())
    }
    pub(super) fn commit(&self) -> Result<()> {
        if unsafe { sqlite3_get_autocommit(self.0) } == 0 {
            self.execute("COMMIT")?;
        }
        Ok(())
    }
}
impl Drop for Connection {
    fn drop(&mut self) {
        // SAFETY: every statement borrows this owner; close rolls back an open
        // default transaction, matching the source's finally: db.close().
        unsafe {
            sqlite3_close_v2(self.0);
        }
    }
}
impl Statement<'_> {
    fn bind(&mut self, values: &[Value]) -> Result<()> {
        let required = unsafe { sqlite3_bind_parameter_count(self.handle) } as usize;
        if values.len() != required {
            return Err(error(
                "ProgrammingError",
                format!(
                    "Incorrect number of bindings supplied. The current statement uses {required}, and there are {} supplied.",
                    values.len()
                ),
            ));
        }
        self.db
            .check(unsafe { sqlite3_clear_bindings(self.handle) })?;
        for (index, value) in values.iter().enumerate() {
            let index = index as c_int + 1;
            // SQLITE_TRANSIENT's public C API sentinel asks SQLite to copy text
            // before this call returns, so no borrowed Rust allocation escapes.
            let transient = unsafe {
                std::mem::transmute::<isize, Option<unsafe extern "C" fn(*mut c_void)>>(-1)
            };
            let code = match value {
                Value::Null => unsafe { sqlite3_bind_null(self.handle, index) },
                Value::Int(n) => unsafe { sqlite3_bind_int64(self.handle, index, *n) },
                Value::Text(text) => {
                    let count = c_int::try_from(text.len())
                        .map_err(|_| error("OverflowError", "string longer than INT_MAX bytes"))?;
                    unsafe {
                        sqlite3_bind_text(
                            self.handle,
                            index,
                            text.as_ptr().cast(),
                            count,
                            transient,
                        )
                    }
                }
                Value::OversizeInteger => {
                    return Err(error(
                        "OverflowError",
                        "Python int too large to convert to SQLite INTEGER",
                    ));
                }
                _ => {
                    return Err(error(
                        "RuntimeError",
                        "Unsupported internal projection parameter",
                    ));
                }
            };
            self.db.check(code)?;
        }
        Ok(())
    }
    pub(super) fn next(&mut self) -> Result<Option<Vec<Value>>> {
        if self.done {
            return Ok(None);
        }
        let code = unsafe { sqlite3_step(self.handle) };
        self.db.check(code)?;
        if code == 101 {
            self.done = true;
            return Ok(None);
        }
        let count = unsafe { sqlite3_column_count(self.handle) };
        let mut row = Vec::with_capacity(count as usize);
        for index in 0..count {
            let kind = unsafe { sqlite3_column_type(self.handle, index) };
            let value = match kind {
                1 => Value::Int(unsafe { sqlite3_column_int64(self.handle, index) }),
                2 => Value::Float(unsafe { sqlite3_column_double(self.handle, index) }),
                5 => Value::Null,
                3 | 4 => {
                    // Column bytes must be asked AFTER text/blob conversion.
                    let pointer = if kind == 3 {
                        unsafe { sqlite3_column_text(self.handle, index) }
                    } else {
                        unsafe { sqlite3_column_blob(self.handle, index).cast() }
                    };
                    let count = unsafe { sqlite3_column_bytes(self.handle, index) } as usize;
                    if pointer.is_null()
                        && (count != 0 || unsafe { sqlite3_errcode(self.db.0) } == 7)
                    {
                        return Err(error("MemoryError", ""));
                    }
                    let bytes = if count == 0 {
                        &[][..]
                    } else {
                        unsafe { std::slice::from_raw_parts(pointer, count) }
                    };
                    if kind == 4 {
                        Value::Blob(bytes.to_vec())
                    } else {
                        let text = std::str::from_utf8(bytes).map_err(|_| {
                            let name = unsafe { c_text(sqlite3_column_name(self.handle, index)) };
                            let mut message =
                                format!("Could not decode to UTF-8 column '{name}' with text '")
                                    .into_bytes();
                            message.extend(bytes.split(|c| *c == 0).next().unwrap_or_default());
                            message.push(b'\'');
                            message.truncate(198);
                            error(
                                "OperationalError",
                                message
                                    .iter()
                                    .map(|c| if c.is_ascii() { *c as char } else { '\u{fffd}' })
                                    .collect::<String>(),
                            )
                        })?;
                        Value::Text(text.to_owned())
                    }
                }
                _ => return Err(error("RuntimeError", "Unknown SQLite scalar type")),
            };
            row.push(value);
        }
        Ok(Some(row))
    }
}
impl Drop for Statement<'_> {
    fn drop(&mut self) {
        unsafe {
            sqlite3_finalize(self.handle);
        }
    }
}
pub(super) fn information() -> Json {
    let mut options = Vec::new();
    let mut index = 0;
    loop {
        let value = unsafe { sqlite3_compileoption_get(index) };
        if value.is_null() {
            break;
        }
        options.push(string(unsafe { c_text(value) }));
        index += 1;
    }
    object([
        ("version", string(unsafe { c_text(sqlite3_libversion()) })),
        ("source_id", string(unsafe { c_text(sqlite3_sourceid()) })),
        ("compile_options", Json::Array(options)),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn engine_keeps_null_empty_multiplicity_and_exact_integer_scalars() {
        let db = Connection::open(Path::new(":memory:")).unwrap();
        db.execute("CREATE TABLE x(k INTEGER,v TEXT)").unwrap();
        db.executemany(
            "INSERT INTO x VALUES (?,?)",
            &[
                vec![Value::Int(1), Value::Null],
                vec![Value::Int(1), Value::Text(String::new())],
                vec![Value::Int(i64::MAX), Value::Text("A".into())],
            ],
        )
        .unwrap();
        db.commit().unwrap();
        let mut query = db.prepare("SELECT k,v FROM x ORDER BY rowid").unwrap();
        assert_eq!(
            query.next().unwrap().unwrap(),
            vec![Value::Int(1), Value::Null]
        );
        assert_eq!(
            query.next().unwrap().unwrap(),
            vec![Value::Int(1), Value::Text(String::new())]
        );
        assert_eq!(
            query.next().unwrap().unwrap(),
            vec![Value::Int(i64::MAX), Value::Text("A".into())]
        );
        assert!(query.next().unwrap().is_none());
        assert_eq!(get(&information(), "version").as_str(), Some("3.46.0"));
    }
    #[test]
    fn binding_overflow_keeps_earlier_pending_rows_until_connection_closes() {
        let db = Connection::open(Path::new(":memory:")).unwrap();
        db.execute("CREATE TABLE x(k INTEGER)").unwrap();
        let cause = db
            .executemany(
                "INSERT INTO x VALUES (?)",
                &[vec![Value::Int(1)], vec![Value::OversizeInteger]],
            )
            .unwrap_err();
        assert_eq!(
            cause.to_string(),
            "OverflowError: Python int too large to convert to SQLite INTEGER"
        );
        assert_eq!(db.scalar("SELECT COUNT(*) FROM x").unwrap(), Value::Int(1));
        assert_eq!(
            db.execute("SELECT 1; SELECT 2").unwrap_err().to_string(),
            "ProgrammingError: You can only execute one statement at a time."
        );
    }
    #[test]
    fn effective_defaults_and_new_file_layout_match_python_engine() {
        let path = std::env::temp_dir().join(format!(
            "priva-sqlite-defaults-{}-{}.sqlite",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        assert!(!path.exists());
        {
            let db = Connection::open(&path).unwrap();
            assert_eq!(db.scalar("PRAGMA auto_vacuum").unwrap(), Value::Int(0));
            assert_eq!(
                db.scalar("PRAGMA recursive_triggers").unwrap(),
                Value::Int(0)
            );
            assert_eq!(db.scalar("PRAGMA secure_delete").unwrap(), Value::Int(1));
            db.execute("CREATE TABLE x(k INTEGER)").unwrap();
            assert_eq!(db.scalar("PRAGMA page_count").unwrap(), Value::Int(2));
        }
        fs::remove_file(path).unwrap();
    }
}
