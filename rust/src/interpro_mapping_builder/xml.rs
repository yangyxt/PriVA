//! Narrow safe owner of the stable narrow-character Expat C API.
//!
//! ElementTree feeds 16 KiB chunks and drains its queued end events AFTER
//! each feed. That timing matters: elem.clear() can remove a nested entry's
//! already parsed tail. Retain it, rather than clearing in the C end callback.
//! No external entity or DTD loader is installed. Expat owns XML grammar;
//! the Rust callbacks only build the selected tree and the source mapping.
use super::{Element, Mapping, record};
use crate::builder_error::category;
use crate::cli::{Result, invalid};
use crate::python_io::GzipFailure;
use std::cell::RefCell;
use std::ffi::{CStr, c_char, c_int, c_void};
use std::io::Read;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::ptr;
use std::rc::Rc;

type Parser = *mut c_void;
type Node = Option<Rc<RefCell<Element>>>;

// ABI declarations from the public expat.h API. We use the default narrow
// XML_Char build (UTF-8 callbacks), never XML_UNICODE or a layout-private type.
#[repr(C)]
struct Encoding {
    map: [c_int; 256],
    data: *mut c_void,
    convert: Option<unsafe extern "C" fn(*mut c_void, *const c_char) -> c_int>,
    release: Option<unsafe extern "C" fn(*mut c_void)>,
}

unsafe extern "C" {
    fn XML_ParserCreateNS(encoding: *const c_char, separator: c_char) -> Parser;
    fn XML_ParserFree(parser: Parser);
    fn XML_SetUserData(parser: Parser, data: *mut c_void);
    fn XML_SetElementHandler(
        parser: Parser,
        start: Option<unsafe extern "C" fn(*mut c_void, *const c_char, *const *const c_char)>,
        end: Option<unsafe extern "C" fn(*mut c_void, *const c_char)>,
    );
    fn XML_SetCharacterDataHandler(
        parser: Parser,
        handler: Option<unsafe extern "C" fn(*mut c_void, *const c_char, c_int)>,
    );
    fn XML_SetDefaultHandlerExpand(
        parser: Parser,
        handler: Option<unsafe extern "C" fn(*mut c_void, *const c_char, c_int)>,
    );
    fn XML_SetExternalEntityRefHandler(
        parser: Parser,
        handler: Option<
            unsafe extern "C" fn(
                Parser,
                *const c_char,
                *const c_char,
                *const c_char,
                *const c_char,
            ) -> c_int,
        >,
    );
    fn XML_SetUnknownEncodingHandler(
        parser: Parser,
        handler: Option<unsafe extern "C" fn(*mut c_void, *const c_char, *mut Encoding) -> c_int>,
        data: *mut c_void,
    );
    fn XML_Parse(parser: Parser, data: *const c_char, length: c_int, final_chunk: c_int) -> c_int;
    fn XML_StopParser(parser: Parser, resumable: u8) -> c_int;
    fn XML_GetErrorCode(parser: Parser) -> c_int;
    fn XML_GetCurrentLineNumber(parser: Parser) -> libc::c_ulong;
    fn XML_GetCurrentColumnNumber(parser: Parser) -> libc::c_ulong;
    fn XML_ErrorString(code: c_int) -> *const c_char;
    fn XML_ExpatVersion() -> *const c_char;
}

struct Owner(Parser);
impl Drop for Owner {
    fn drop(&mut self) {
        // SAFETY: Owner is created only for a nonnull unique parser, and no
        // callback/parse can run after this synchronous owner's destruction.
        unsafe { XML_ParserFree(self.0) };
    }
}

struct State {
    parser: Parser,
    stack: Vec<Node>,
    ends: Vec<Node>,
    first_end: bool,
    result: Mapping,
    failure: Option<Box<dyn std::error::Error>>,
    product6: Option<Product6>,
}

#[derive(Default)]
struct Product6 {
    root: Option<String>,
    has_disorders: bool,
}

impl State {
    fn drain(&mut self) {
        for node in std::mem::take(&mut self.ends) {
            if self.first_end {
                if let Some(element) = node {
                    if element.borrow().name == "interpro" {
                        record(&element.borrow(), &mut self.result);
                        element.borrow_mut().clear();
                    }
                }
            }
            self.first_end = true;
        }
    }
}

/// Expat strings are valid NUL-terminated UTF-8 during the callback only.
unsafe fn text(pointer: *const c_char) -> Result<String> {
    // SAFETY: callers pass an Expat callback name/value or static error text.
    Ok(unsafe { CStr::from_ptr(pointer) }.to_str()?.to_owned())
}

fn expanded(name: String) -> String {
    // The same namespace separator used by ElementTree. Namespaced tags never
    // compare equal to the source's literal unqualified names.
    if name.contains('}') {
        format!("{{{name}")
    } else {
        name
    }
}

unsafe fn callback(data: *mut c_void, operation: impl FnOnce(&mut State) -> Result<()>) {
    // SAFETY: data points to the boxed State, which outlives Owner; callbacks
    // are synchronous and Expat never calls two callbacks simultaneously.
    let state = unsafe { &mut *data.cast::<State>() };
    if state.failure.is_some() {
        return;
    }
    let result = catch_unwind(AssertUnwindSafe(|| operation(state)));
    let error = match result {
        Ok(Ok(())) => return,
        Ok(Err(error)) => error,
        Err(_) => category("RuntimeError", invalid("XML callback panicked")),
    };
    state.failure = Some(error);
    // SAFETY: this callback runs inside this live parser's XML_Parse call.
    unsafe { XML_StopParser(state.parser, 0) };
}

unsafe extern "C" fn start(
    data: *mut c_void,
    name: *const c_char,
    attributes: *const *const c_char,
) {
    unsafe {
        callback(data, |state| {
            let name = expanded(text(name)?);
            if let Some(product) = &mut state.product6 {
                if name.contains('}') {
                    return Err(category(
                        "ValueError",
                        invalid(
                            "The original Orphadata reader requires unnamespaced Product 6 XML",
                        ),
                    ));
                }
                if product.root.is_none() {
                    product.root = Some(name.clone());
                }
                product.has_disorders |= name == "DisorderList";
                return Ok(());
            }
            let capture = name == "interpro" || state.stack.last().is_some_and(Option::is_some);
            let node = if capture {
                let mut element = Element {
                    name,
                    ..Element::default()
                };
                let mut at = 0;
                // Expat supplies NUL-terminated (name,value) pairs followed by a
                // null name pointer. All pointed-to text lives for this callback.
                while !(*attributes.add(at)).is_null() {
                    element.attributes.insert(
                        expanded(text(*attributes.add(at))?),
                        text(*attributes.add(at + 1))?,
                    );
                    at += 2;
                }
                let node = Rc::new(RefCell::new(element));
                if let Some(Some(parent)) = state.stack.last() {
                    parent.borrow_mut().children.push(node.clone());
                }
                Some(node)
            } else {
                None
            };
            state.stack.push(node);
            Ok(())
        })
    }
}

unsafe extern "C" fn end(data: *mut c_void, _name: *const c_char) {
    unsafe {
        callback(data, |state| {
            if state.product6.is_some() {
                return Ok(());
            }
            let node = state
                .stack
                .pop()
                .ok_or_else(|| invalid("XML element stack underflow"))?;
            state.ends.push(node);
            Ok(())
        })
    }
}

unsafe extern "C" fn characters(data: *mut c_void, value: *const c_char, length: c_int) {
    unsafe {
        callback(data, |state| {
            if state.product6.is_some() {
                return Ok(());
            }
            // Expat's character payload is length-delimited, not NUL-terminated.
            // The supplied nonnegative byte count covers valid callback memory.
            let bytes = std::slice::from_raw_parts(value.cast::<u8>(), usize::try_from(length)?);
            if let Some(Some(node)) = state.stack.last() {
                node.borrow_mut().append_text(std::str::from_utf8(bytes)?);
            }
            Ok(())
        })
    }
}

unsafe extern "C" fn external_entity(
    _: Parser,
    _: *const c_char,
    _: *const c_char,
    _: *const c_char,
    _: *const c_char,
) -> c_int {
    // ElementTree rejects an external entity reference; no filesystem/network
    // resolver is provided. External DTD declarations alone are not loaded.
    0
}

unsafe extern "C" fn product6_default(data: *mut c_void, value: *const c_char, length: c_int) {
    unsafe {
        callback(data, |state| {
            // ElementTree's expanding default handler rejects an external
            // entity by its literal reference; no external resolver is used.
            let bytes = std::slice::from_raw_parts(value.cast::<u8>(), usize::try_from(length)?);
            if bytes.starts_with(b"&") {
                let value = std::str::from_utf8(bytes)?;
                let line = XML_GetCurrentLineNumber(state.parser);
                let column = XML_GetCurrentColumnNumber(state.parser);
                return Err(category(
                    "ParseError",
                    invalid(format!(
                        "undefined entity {value}: line {line}, column {column}"
                    )),
                ));
            }
            Ok(())
        })
    }
}

unsafe extern "C" fn unknown_encoding(
    data: *mut c_void,
    name: *const c_char,
    info: *mut Encoding,
) -> c_int {
    let mut accepted = false;
    unsafe {
        callback(data, |state| {
            let name = text(name)?;
            // Python normalizes ASCII codec labels before standard-library
            // lookup. WHATWG aliases differ (notably Latin-1 and UTF-16-LE), so
            // retain the pinned CPython whole-byte-map fallback instead.
            let mut normalized = String::new();
            let mut separator = false;
            for byte in name.bytes() {
                if byte.is_ascii_alphanumeric() || byte == b'.' {
                    if separator && !normalized.is_empty() {
                        normalized.push('_');
                    }
                    normalized.push(byte.to_ascii_lowercase() as char);
                    separator = false;
                } else {
                    separator = true;
                }
            }
            let Some(codec) = super::encodings::lookup(&normalized) else {
                state.failure = Some(category(
                    "LookupError",
                    invalid(format!("unknown encoding: {name}")),
                ));
                return Ok(());
            };
            match codec {
                super::encodings::Codec::Map(map) => (*info).map = *map,
                super::encodings::Codec::Error(kind, message) => {
                    state.failure = Some(category(kind, invalid(message)));
                    return Ok(());
                }
            }
            (*info).data = ptr::null_mut();
            (*info).convert = None;
            (*info).release = None;
            accepted = true;
            Ok(())
        })
    }
    c_int::from(accepted)
}

pub(super) fn version() -> String {
    // SAFETY: Expat returns a process-lifetime static C string.
    unsafe { CStr::from_ptr(XML_ExpatVersion()) }
        .to_string_lossy()
        .into_owned()
}

pub(super) fn parse(source: impl Read) -> Result<Mapping> {
    parse_with(source, false)
}

pub(super) fn validate_product6(source: impl Read) -> Result<()> {
    parse_with(source, true).map(|_| ())
}

fn parse_with(mut source: impl Read, product6: bool) -> Result<Mapping> {
    // Declare state first so every early return drops the parser before its
    // callback state. The box keeps its address stable while parsing.
    let mut state = Box::new(State {
        parser: ptr::null_mut(),
        stack: vec![],
        ends: vec![],
        first_end: false,
        result: Mapping::new(),
        failure: None,
        product6: product6.then(Product6::default),
    });
    // SAFETY: null selects XML-declared/default encoding; namespace separator
    // is ElementTree's ASCII '}'. Null is checked before all other API calls.
    let parser = Owner(unsafe { XML_ParserCreateNS(ptr::null(), b'}' as c_char) });
    if parser.0.is_null() {
        std::mem::forget(parser);
        return Err(category(
            "MemoryError",
            invalid("could not create XML parser"),
        ));
    }
    state.parser = parser.0;
    let data = ptr::from_mut(state.as_mut()).cast::<c_void>();
    // SAFETY: the stable boxed data lives until all parsing is finished; these
    // callbacks borrow only during each synchronous XML_Parse invocation.
    unsafe {
        XML_SetUserData(parser.0, data);
        XML_SetElementHandler(parser.0, Some(start), Some(end));
        XML_SetCharacterDataHandler(parser.0, Some(characters));
        if product6 {
            XML_SetDefaultHandlerExpand(parser.0, Some(product6_default));
        } else {
            XML_SetExternalEntityRefHandler(parser.0, Some(external_entity));
        }
        XML_SetUnknownEncodingHandler(parser.0, Some(unknown_encoding), data);
    }
    let mut buffer = [0u8; 16 * 1024];
    loop {
        // gzip.GzipFile.read(16*1024), not a short underlying read: the queued
        // mutation boundary is defined by the decompressed byte count.
        let mut used = 0;
        while used < buffer.len() {
            let count = source
                .read(&mut buffer[used..])
                .map_err(|error| match GzipFailure::of(&error) {
                    Some(GzipFailure::EndOfFile) => category("EOFError", error.into()),
                    Some(GzipFailure::Bad(_)) => category("BadGzipFile", error.into()),
                    Some(GzipFailure::Zlib(_)) => category("error", error.into()),
                    None => category("OSError", error.into()),
                })?;
            if count == 0 {
                break;
            }
            used += count;
        }
        // SAFETY: byte buffer is valid for used <= 16384 bytes throughout the
        // synchronous parse; callbacks cannot retain its borrowed bytes.
        let status = unsafe {
            XML_Parse(
                parser.0,
                buffer.as_ptr().cast(),
                used as c_int,
                c_int::from(used == 0),
            )
        };
        state.drain();
        if let Some(error) = state.failure.take() {
            return Err(error);
        }
        if status != 1 {
            let reason = unsafe { text(XML_ErrorString(XML_GetErrorCode(parser.0))) }?;
            if product6 {
                // ElementTree includes Expat's original line/column for this
                // validator; retain the InterPro builder's existing diagnostic.
                let (line, column) = unsafe {
                    (
                        XML_GetCurrentLineNumber(parser.0),
                        XML_GetCurrentColumnNumber(parser.0),
                    )
                };
                return Err(category(
                    "ParseError",
                    invalid(format!("{reason}: line {line}, column {column}")),
                ));
            }
            return Err(category("ParseError", invalid(reason)));
        }
        if used == 0 {
            break;
        }
    }
    // Free the C owner before the callback's state, including on success.
    drop(parser);
    if let Some(product) = &state.product6 {
        if product.root.as_deref() != Some("JDBOR") || !product.has_disorders {
            return Err(category(
                "ValueError",
                invalid("Expected Orphadata Product 6 JDBOR/DisorderList XML"),
            ));
        }
    }
    Ok(std::mem::take(&mut state.result))
}
