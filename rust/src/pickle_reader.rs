//! A restricted reader for the Python pickles PriVA stores (the InterPro
//! mapping, the transcript domain map and the intolerant domain set).
//!
//! A pickle is a small program for a stack machine, and Python's `pickle.load`
//! will run any importable function a file names. This reader instead handles
//! exactly the opcodes that plain data needs (None, booleans, integers up to
//! 64 bits, floats, text, tuples, lists, dicts, sets, frozen sets, and the
//! memo that lets one object appear twice) and stops with an error naming the
//! opcode and its byte offset on anything else: no GLOBAL, REDUCE, BUILD,
//! bytes or numpy objects. Protocols 2 to 5 as written by Python 3 are
//! covered; sets pickled with protocol 2 or 3 (which need REDUCE) are not.
//!
//! Containers keep Python's behaviour that matters here: a dict keeps the
//! first position of a key and the last value written to it, a set holds each
//! item once, and `1`, `1.0` and `True` are the same key. A list, dict or set
//! that appears twice in the pickle is the same object in both places.
use crate::cli::{Result, invalid};
use flate2::read::MultiGzDecoder;
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::rc::Rc;

#[derive(Clone, Debug)]
pub(crate) enum Value {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(Rc<str>),
    Tuple(Rc<Vec<Value>>),
    List(Rc<RefCell<Vec<Value>>>),
    Dict(Rc<RefCell<Dict>>),
    Set(Rc<RefCell<Set>>),
    FrozenSet(Rc<Set>),
}

#[derive(Debug, Default)]
pub(crate) struct Dict {
    pub(crate) entries: Vec<(Value, Value)>,
    positions: HashMap<String, usize>,
}

impl Dict {
    fn insert(&mut self, key: Value, value: Value) -> Result<()> {
        let identity = identity(&key)?;
        match self.positions.get(&identity) {
            Some(&position) => self.entries[position].1 = value,
            None => {
                self.positions.insert(identity, self.entries.len());
                self.entries.push((key, value));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Default)]
pub(crate) struct Set {
    pub(crate) items: Vec<Value>,
    seen: HashSet<String>,
}

impl Set {
    fn insert(&mut self, item: Value) -> Result<()> {
        if self.seen.insert(identity(&item)?) {
            self.items.push(item);
        }
        Ok(())
    }
}

/// A text that is equal exactly when Python considers two hashable values the
/// same dict key or set member.
fn identity(value: &Value) -> Result<String> {
    Ok(match value {
        Value::None => "N".to_owned(),
        Value::Bool(flag) => format!("n{}", i64::from(*flag)),
        Value::Int(number) => format!("n{number}"),
        Value::Float(number) => {
            if number.is_nan() {
                "nan".to_owned()
            } else if number.fract() == 0.0 && number.abs() < 9.0e18 {
                format!("n{}", *number as i64)
            } else {
                format!("f{:016x}", number.to_bits())
            }
        }
        Value::Str(text) => format!("s{text}"),
        Value::Tuple(items) => {
            let mut text = String::from("t");
            for item in items.iter() {
                let part = identity(item)?;
                text.push_str(&format!("{}:{part}", part.len()));
            }
            text
        }
        Value::FrozenSet(set) => {
            let mut parts = Vec::new();
            for item in &set.items {
                parts.push(identity(item)?);
            }
            parts.sort();
            let mut text = String::from("z");
            for part in parts {
                text.push_str(&format!("{}:{part}", part.len()));
            }
            text
        }
        Value::List(_) | Value::Dict(_) | Value::Set(_) => {
            return Err(invalid(
                "the pickle uses a list, dict or set as a dict key or set member",
            ));
        }
    })
}

impl Value {
    /// Python's `type(value).__name__`.
    pub(crate) fn type_name(&self) -> &'static str {
        match self {
            Value::None => "NoneType",
            Value::Bool(_) => "bool",
            Value::Int(_) => "int",
            Value::Float(_) => "float",
            Value::Str(_) => "str",
            Value::Tuple(_) => "tuple",
            Value::List(_) => "list",
            Value::Dict(_) => "dict",
            Value::Set(_) => "set",
            Value::FrozenSet(_) => "frozenset",
        }
    }

    /// Python's `len(value)` for the sized types, `None` for the others.
    pub(crate) fn length(&self) -> Option<usize> {
        match self {
            Value::Str(text) => Some(text.chars().count()),
            Value::Tuple(items) => Some(items.len()),
            Value::List(items) => Some(items.borrow().len()),
            Value::Dict(dict) => Some(dict.borrow().entries.len()),
            Value::Set(set) => Some(set.borrow().items.len()),
            Value::FrozenSet(set) => Some(set.items.len()),
            _ => None,
        }
    }

    pub(crate) fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(text) => Some(&**text),
            _ => None,
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8]> {
        let end = self
            .position
            .checked_add(count)
            .filter(|end| *end <= self.bytes.len())
            .ok_or_else(|| invalid("the pickle ends before its STOP opcode"))?;
        let all: &'a [u8] = self.bytes;
        let slice = &all[self.position..end];
        self.position = end;
        Ok(slice)
    }

    fn byte(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    /// A little-endian unsigned integer of up to eight bytes.
    fn unsigned(&mut self, width: usize) -> Result<u64> {
        let bytes = self.take(width)?;
        let mut value = 0_u64;
        for (index, byte) in bytes.iter().enumerate() {
            value |= u64::from(*byte) << (8 * index);
        }
        Ok(value)
    }
}

fn read_long(reader: &mut Reader, width: usize) -> Result<i64> {
    let length = reader.unsigned(width)? as usize;
    if length == 0 {
        return Ok(0);
    }
    if length > 8 {
        return Err(invalid("the pickle holds an integer wider than 64 bits"));
    }
    let bytes = reader.take(length)?;
    let mut value: i64 = if bytes[length - 1] & 0x80 != 0 { -1 } else { 0 };
    for byte in bytes.iter().rev() {
        value = (value << 8) | i64::from(*byte);
    }
    Ok(value)
}

fn read_text(reader: &mut Reader, length: u64) -> Result<Value> {
    let length = usize::try_from(length).map_err(|_| invalid("the pickle holds a huge string"))?;
    let bytes = reader.take(length)?;
    let text = std::str::from_utf8(bytes).map_err(|_| {
        invalid(format!(
            "the pickle holds a string that is not UTF-8 ending at byte {}",
            reader.position
        ))
    })?;
    Ok(Value::Str(Rc::from(text)))
}

fn pop(stack: &mut Vec<Value>) -> Result<Value> {
    stack
        .pop()
        .ok_or_else(|| invalid("the pickle takes a value from an empty stack"))
}

/// The values pushed since the last MARK, removed from the stack.
fn marked(stack: &mut Vec<Value>, marks: &mut Vec<usize>) -> Result<Vec<Value>> {
    let mark = marks
        .pop()
        .ok_or_else(|| invalid("the pickle closes a MARK it never opened"))?;
    if mark > stack.len() {
        return Err(invalid("the pickle's MARK lies above the stack"));
    }
    Ok(stack.split_off(mark))
}

fn rejected(opcode: u8, offset: usize) -> Box<dyn std::error::Error> {
    let name = match opcode {
        b'c' => "GLOBAL",
        0x93 => "STACK_GLOBAL",
        b'R' => "REDUCE",
        b'b' => "BUILD",
        0x81 => "NEWOBJ",
        0x92 => "NEWOBJ_EX",
        b'B' | b'C' | 0x8e | 0x96 => "a bytes opcode",
        b'I' | b'L' | b'F' | b'V' | b'S' => "a text-protocol opcode (protocol 0)",
        _ => "an unknown opcode",
    };
    invalid(format!(
        "the pickle uses opcode 0x{opcode:02x} ({name}) at byte {offset}, which the restricted \
         reader does not accept; it reads plain data only"
    ))
}

/// Reads one pickle from the start of `bytes`; anything after its STOP opcode
/// is ignored, as `pickle.load` ignores it.
pub(crate) fn load(bytes: &[u8]) -> Result<Value> {
    let mut reader = Reader { bytes, position: 0 };
    let mut stack: Vec<Value> = Vec::new();
    let mut marks: Vec<usize> = Vec::new();
    let mut memo: HashMap<u64, Value> = HashMap::new();
    loop {
        let offset = reader.position;
        let opcode = reader.byte()?;
        match opcode {
            0x80 => {
                let version = reader.byte()?;
                if version > 5 {
                    return Err(invalid(format!(
                        "pickle protocol {version} is not supported"
                    )));
                }
            }
            0x95 => {
                reader.take(8)?;
            }
            b'.' => return pop(&mut stack),
            b'(' => marks.push(stack.len()),
            b'N' => stack.push(Value::None),
            0x88 => stack.push(Value::Bool(true)),
            0x89 => stack.push(Value::Bool(false)),
            b'J' => {
                let number = reader.unsigned(4)? as u32 as i32;
                stack.push(Value::Int(i64::from(number)));
            }
            b'K' => {
                let number = reader.unsigned(1)?;
                stack.push(Value::Int(number as i64));
            }
            b'M' => {
                let number = reader.unsigned(2)?;
                stack.push(Value::Int(number as i64));
            }
            0x8a => {
                let number = read_long(&mut reader, 1)?;
                stack.push(Value::Int(number));
            }
            0x8b => {
                let number = read_long(&mut reader, 4)?;
                stack.push(Value::Int(number));
            }
            b'G' => {
                let mut array = [0_u8; 8];
                array.copy_from_slice(reader.take(8)?);
                stack.push(Value::Float(f64::from_be_bytes(array)));
            }
            0x8c => {
                let length = reader.unsigned(1)?;
                stack.push(read_text(&mut reader, length)?);
            }
            b'X' => {
                let length = reader.unsigned(4)?;
                stack.push(read_text(&mut reader, length)?);
            }
            0x8d => {
                let length = reader.unsigned(8)?;
                stack.push(read_text(&mut reader, length)?);
            }
            b']' => stack.push(Value::List(Rc::new(RefCell::new(Vec::new())))),
            b'a' => {
                let item = pop(&mut stack)?;
                match stack.last() {
                    Some(Value::List(list)) => list.borrow_mut().push(item),
                    _ => return Err(invalid(format!("APPEND at byte {offset} needs a list"))),
                }
            }
            b'e' => {
                let items = marked(&mut stack, &mut marks)?;
                match stack.last() {
                    Some(Value::List(list)) => list.borrow_mut().extend(items),
                    _ => return Err(invalid(format!("APPENDS at byte {offset} needs a list"))),
                }
            }
            b'}' => stack.push(Value::Dict(Rc::new(RefCell::new(Dict::default())))),
            b's' => {
                let value = pop(&mut stack)?;
                let key = pop(&mut stack)?;
                match stack.last() {
                    Some(Value::Dict(dict)) => dict.borrow_mut().insert(key, value)?,
                    _ => return Err(invalid(format!("SETITEM at byte {offset} needs a dict"))),
                }
            }
            b'u' => {
                let items = marked(&mut stack, &mut marks)?;
                if items.len() % 2 != 0 {
                    return Err(invalid(format!("SETITEMS at byte {offset} has a lone key")));
                }
                match stack.last() {
                    Some(Value::Dict(dict)) => {
                        let mut dict = dict.borrow_mut();
                        let mut pairs = items.into_iter();
                        while let (Some(key), Some(value)) = (pairs.next(), pairs.next()) {
                            dict.insert(key, value)?;
                        }
                    }
                    _ => return Err(invalid(format!("SETITEMS at byte {offset} needs a dict"))),
                }
            }
            b')' => stack.push(Value::Tuple(Rc::new(Vec::new()))),
            b't' => {
                let items = marked(&mut stack, &mut marks)?;
                stack.push(Value::Tuple(Rc::new(items)));
            }
            0x85 => {
                let first = pop(&mut stack)?;
                stack.push(Value::Tuple(Rc::new(vec![first])));
            }
            0x86 => {
                let second = pop(&mut stack)?;
                let first = pop(&mut stack)?;
                stack.push(Value::Tuple(Rc::new(vec![first, second])));
            }
            0x87 => {
                let third = pop(&mut stack)?;
                let second = pop(&mut stack)?;
                let first = pop(&mut stack)?;
                stack.push(Value::Tuple(Rc::new(vec![first, second, third])));
            }
            0x8f => stack.push(Value::Set(Rc::new(RefCell::new(Set::default())))),
            0x90 => {
                let items = marked(&mut stack, &mut marks)?;
                match stack.last() {
                    Some(Value::Set(set)) => {
                        let mut set = set.borrow_mut();
                        for item in items {
                            set.insert(item)?;
                        }
                    }
                    _ => return Err(invalid(format!("ADDITEMS at byte {offset} needs a set"))),
                }
            }
            0x91 => {
                let items = marked(&mut stack, &mut marks)?;
                let mut set = Set::default();
                for item in items {
                    set.insert(item)?;
                }
                stack.push(Value::FrozenSet(Rc::new(set)));
            }
            0x94 => {
                let top = stack
                    .last()
                    .cloned()
                    .ok_or_else(|| invalid("MEMOIZE with an empty stack"))?;
                let index = memo.len() as u64;
                memo.insert(index, top);
            }
            b'q' | b'r' => {
                let index = reader.unsigned(if opcode == b'q' { 1 } else { 4 })?;
                let top = stack
                    .last()
                    .cloned()
                    .ok_or_else(|| invalid("BINPUT with an empty stack"))?;
                memo.insert(index, top);
            }
            b'h' | b'j' => {
                let index = reader.unsigned(if opcode == b'h' { 1 } else { 4 })?;
                let item = memo.get(&index).cloned().ok_or_else(|| {
                    invalid(format!("BINGET at byte {offset} names an unset memo slot"))
                })?;
                stack.push(item);
            }
            _ => return Err(rejected(opcode, offset)),
        }
    }
}

/// The bytes of a pickle file, decompressed when `gzip` says so.
pub(crate) fn read_file(path: &Path, gzip: bool) -> Result<Vec<u8>> {
    let mut file = File::open(path)?;
    let mut bytes = Vec::new();
    if gzip {
        MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else {
        file.read_to_end(&mut bytes)?;
    }
    Ok(bytes)
}

/// The bytes of a pickle file that is decompressed when it starts with the
/// gzip magic bytes, whatever its name is.
pub(crate) fn read_file_by_magic(path: &Path) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    File::open(path)?.read_to_end(&mut bytes)?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::new();
        MultiGzDecoder::new(&bytes[..]).read_to_end(&mut plain)?;
        Ok(plain)
    } else {
        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_hex(text: &str) -> Vec<u8> {
        (0..text.len())
            .step_by(2)
            .map(|index| u8::from_str_radix(&text[index..index + 2], 16).unwrap())
            .collect()
    }

    // pickle.dumps(obj, protocol=4) for the object described in `plain_data_is_read`.
    const PLAIN: &str = "80049568010000000000007d94288c045046414d948c07504630303030319486945d94288c09495052303030303031948c06446f6d61696e948c074b72696e676c659468075d948c04474f3a31948c126d6f6c6563756c61725f66756e6374696f6e948c0762696e64696e67948794618c03616273947494618c04696e7473945d94284b004b074d2c014a701101008a060000000000014afbffffff4ad4feffff8a060000000000ff658c06666c6f61747394473fb999999999999a47c004000000000000477ff000000000000087948c05666c616773945d942888894e658c0474657874948c05c3a9e4b8ad948c05656d70747994298c027478947d948c05454e535431947d94288c0a70705f646f6d61696e73945d944b014b148c0164948794618c0c65786f6e5f646f6d61696e73947d948c0131945d948c0a673a5066616d3a50463194617375738c076d656d62657273948f94288c0161948c016294908c0666726f7a656e94288c0178949194752e";

    fn entry<'a>(dict: &'a Dict, name: &str) -> &'a Value {
        &dict
            .entries
            .iter()
            .find(|(key, _)| key.as_str() == Some(name))
            .unwrap()
            .1
    }

    #[test]
    fn plain_data_is_read() {
        // {("PFAM", "PF00001"): [("IPR000001", "Domain", "Kringle", "Kringle",
        //   [("GO:1", "molecular_function", "binding")], "abs")],
        //  "ints": [0, 7, 300, 70000, 2**40, -5, -300, -2**40],
        //  "floats": (0.1, -2.5, inf), "flags": [True, False, None],
        //  "text": "é中", "empty": (),
        //  "tx": {"ENST1": {"pp_domains": [(1, 20, "d")], "exon_domains": {"1": ["g:Pfam:PF1"]}}},
        //  "members": {"b", "a"}, "frozen": frozenset(["x"])}
        let value = load(&from_hex(PLAIN)).unwrap();
        assert_eq!(value.type_name(), "dict");
        assert_eq!(value.length(), Some(9));
        let Value::Dict(dict) = value else { panic!() };
        let dict = dict.borrow();
        let Value::List(ints) = entry(&dict, "ints") else {
            panic!()
        };
        let ints: Vec<i64> = ints
            .borrow()
            .iter()
            .map(|item| match item {
                Value::Int(number) => *number,
                _ => panic!(),
            })
            .collect();
        assert_eq!(ints, [0, 7, 300, 70000, 1 << 40, -5, -300, -(1 << 40)]);
        let Value::Tuple(floats) = entry(&dict, "floats") else {
            panic!()
        };
        assert!(matches!(floats[0], Value::Float(number) if number == 0.1));
        assert!(matches!(floats[2], Value::Float(number) if number == f64::INFINITY));
        assert_eq!(entry(&dict, "text").as_str(), Some("é中"));
        assert_eq!(entry(&dict, "empty").length(), Some(0));
        assert_eq!(entry(&dict, "members").length(), Some(2));
        assert_eq!(entry(&dict, "frozen").type_name(), "frozenset");
        assert_eq!(entry(&dict, "flags").length(), Some(3));
        let Value::Tuple(key) = &dict.entries[0].0 else {
            panic!()
        };
        assert_eq!(key[0].as_str(), Some("PFAM"));
        let Value::List(entries) = &dict.entries[0].1 else {
            panic!()
        };
        let entries = entries.borrow();
        let Value::Tuple(domain) = &entries[0] else {
            panic!()
        };
        assert_eq!(domain.len(), 6);
        let Value::List(go_terms) = &domain[4] else {
            panic!()
        };
        assert_eq!(go_terms.borrow().len(), 1);
    }

    #[test]
    fn a_shared_list_is_one_object() {
        // pickle.dumps({"first": shared, "second": shared, "t": (shared,)}, protocol=4)
        // with shared = [1, 2]
        let bytes = from_hex(
            "80049528000000000000007d94288c056669727374945d94284b014b02658c067365636f6e649468028c01749468028594752e",
        );
        let Value::Dict(dict) = load(&bytes).unwrap() else {
            panic!()
        };
        let dict = dict.borrow();
        let Value::List(first) = entry(&dict, "first") else {
            panic!()
        };
        let Value::List(second) = entry(&dict, "second") else {
            panic!()
        };
        assert!(Rc::ptr_eq(first, second));
        assert_eq!(first.borrow().len(), 2);
    }

    #[test]
    fn protocol_two_with_memo_puts_is_read() {
        // pickle.dumps({"a": [1]}, protocol=2)
        let value = load(&from_hex("80027d710058010000006171015d71024b0161732e")).unwrap();
        assert_eq!(value.length(), Some(1));
    }

    #[test]
    fn anything_else_is_refused_by_name() {
        // GLOBAL builtins set, as protocol 2 writes a set.
        let error = load(b"\x80\x02cbuiltins\nset\nq\x00.")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("GLOBAL") && error.contains("byte 2"),
            "{error}"
        );
        assert!(load(b"\x80\x04N").is_err(), "no STOP");
        assert!(load(b"\x80\x06N.").is_err(), "protocol 6");
        assert!(load(b"not a pickle").is_err());
        assert!(load(b"").is_err());
    }

    #[test]
    fn dict_keys_follow_python_equality() {
        // {1: "a", 1.0: "b", True: "c"} is one entry whose key stays 1 and value is "c".
        let mut dict = Dict::default();
        dict.insert(Value::Int(1), Value::Str(Rc::from("a")))
            .unwrap();
        dict.insert(Value::Float(1.0), Value::Str(Rc::from("b")))
            .unwrap();
        dict.insert(Value::Bool(true), Value::Str(Rc::from("c")))
            .unwrap();
        assert_eq!(dict.entries.len(), 1);
        assert!(matches!(dict.entries[0].0, Value::Int(1)));
        assert_eq!(dict.entries[0].1.as_str(), Some("c"));
        assert!(
            dict.insert(Value::List(Rc::new(RefCell::new(Vec::new()))), Value::None)
                .is_err()
        );
    }
}
