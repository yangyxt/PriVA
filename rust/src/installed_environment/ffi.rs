//! Narrow owned-event boundary to the mature libyaml parser and emitter.
use super::*;
use std::ffi::{CStr, c_char, c_int, c_void};
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Mark {
    pub index: usize,
    pub line: usize,
    pub column: usize,
}
#[repr(C)]
#[derive(Default)]
struct EventView {
    kind: c_int,
    anchor: *const u8,
    tag: *const u8,
    value: *const u8,
    length: usize,
    style: c_int,
    implicit: c_int,
    quoted: c_int,
    index: usize,
    line: usize,
    column: usize,
}
#[repr(C)]
struct ErrorView {
    kind: c_int,
    context: *const c_char,
    problem: *const c_char,
    context_index: usize,
    context_line: usize,
    context_column: usize,
    problem_index: usize,
    problem_line: usize,
    problem_column: usize,
    offset: usize,
    value: c_int,
}
#[derive(Clone, Debug, Default)]
pub(super) struct Event {
    pub kind: i32,
    pub anchor: String,
    pub tag: String,
    pub value: String,
    pub implicit: bool,
    pub mark: Mark,
}
#[derive(Debug)]
pub(super) struct ParseError {
    pub kind: i32,
    pub context: String,
    pub problem: String,
    pub context_mark: Mark,
    pub problem_mark: Mark,
    pub offset: usize,
    pub value: i32,
}
unsafe extern "C" {
    fn priva_yaml_parse(
        input: *const u8,
        length: usize,
        callback: unsafe extern "C" fn(*mut c_void, *const EventView) -> c_int,
        failure: unsafe extern "C" fn(*mut c_void, *const ErrorView),
        context: *mut c_void,
    ) -> c_int;
    fn priva_yaml_emitter_new(
        callback: unsafe extern "C" fn(*mut c_void, *mut u8, usize) -> c_int,
        context: *mut c_void,
    ) -> *mut c_void;
    fn priva_yaml_emit(handle: *mut c_void, event: *const EventView) -> c_int;
    fn priva_yaml_emitter_delete(handle: *mut c_void);
}
fn c_text(value: *const u8) -> String {
    if value.is_null() {
        String::new()
    } else {
        // libyaml emits valid UTF-8; event storage lives until the callback returns.
        unsafe { CStr::from_ptr(value.cast()).to_string_lossy().into_owned() }
    }
}
struct Parser<'a> {
    callback: &'a mut dyn FnMut(Event) -> bool,
    failure: Option<ParseError>,
}
unsafe extern "C" fn event(context: *mut c_void, event: *const EventView) -> c_int {
    let state = unsafe { &mut *context.cast::<Parser<'_>>() };
    let event = unsafe { &*event };
    let value = if event.value.is_null() {
        String::new()
    } else {
        String::from_utf8_lossy(unsafe { std::slice::from_raw_parts(event.value, event.length) })
            .into_owned()
    };
    i32::from((state.callback)(Event {
        kind: event.kind,
        anchor: c_text(event.anchor),
        tag: c_text(event.tag),
        value,
        implicit: event.implicit != 0,
        mark: Mark {
            index: event.index,
            line: event.line,
            column: event.column,
        },
    }))
}
unsafe extern "C" fn failure(context: *mut c_void, error: *const ErrorView) {
    let state = unsafe { &mut *context.cast::<Parser<'_>>() };
    let e = unsafe { &*error };
    state.failure = Some(ParseError {
        kind: e.kind,
        context: c_text(e.context.cast()),
        problem: c_text(e.problem.cast()),
        context_mark: Mark {
            index: e.context_index,
            line: e.context_line,
            column: e.context_column,
        },
        problem_mark: Mark {
            index: e.problem_index,
            line: e.problem_line,
            column: e.problem_column,
        },
        offset: e.offset,
        value: e.value,
    });
}
pub(super) fn parse(
    text: &str,
    callback: &mut dyn FnMut(Event) -> bool,
) -> std::result::Result<(), ParseError> {
    let mut parser = Parser {
        callback,
        failure: None,
    };
    let result = unsafe {
        priva_yaml_parse(
            text.as_ptr(),
            text.len(),
            event,
            failure,
            (&mut parser as *mut Parser<'_>).cast(),
        )
    };
    if result == -1 {
        std::alloc::handle_alloc_error(std::alloc::Layout::new::<Parser<'_>>());
    }
    match parser.failure {
        Some(error) => Err(error),
        None => Ok(()),
    }
}
unsafe extern "C" fn output(context: *mut c_void, buffer: *mut u8, length: usize) -> c_int {
    let output = unsafe { &mut *context.cast::<Vec<u8>>() };
    output.extend_from_slice(unsafe { std::slice::from_raw_parts(buffer, length) });
    1
}
pub(super) fn emit(events: &[Event]) -> Result<String> {
    let mut bytes = Vec::new();
    let emitter = unsafe { priva_yaml_emitter_new(output, (&mut bytes as *mut Vec<u8>).cast()) };
    if emitter.is_null() {
        return Err(error("MemoryError", "libyaml emitter allocation failed"));
    }
    struct Guard(*mut c_void);
    impl Drop for Guard {
        fn drop(&mut self) {
            unsafe { priva_yaml_emitter_delete(self.0) };
        }
    }
    let _guard = Guard(emitter);
    for e in events {
        let anchor = std::ffi::CString::new(e.anchor.as_bytes())?;
        let tag = std::ffi::CString::new(e.tag.as_bytes())?;
        if e.value.len() > i32::MAX as usize {
            return Err(error("OverflowError", "YAML scalar is too large"));
        }
        let view = EventView {
            kind: e.kind,
            anchor: if e.anchor.is_empty() {
                std::ptr::null()
            } else {
                anchor.as_ptr().cast()
            },
            tag: if e.tag.is_empty() {
                std::ptr::null()
            } else {
                tag.as_ptr().cast()
            },
            value: e.value.as_ptr(),
            length: e.value.len(),
            implicit: i32::from(e.implicit),
            quoted: 1,
            ..EventView::default()
        };
        if unsafe { priva_yaml_emit(emitter, &view) } == 0 {
            return Err(error(
                "EmitterError",
                "libyaml could not emit the package export",
            ));
        }
    }
    String::from_utf8(bytes).map_err(Into::into)
}
