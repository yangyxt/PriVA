//! Read serialized exon-set insertion sequences from retained Alpha pickles.
//!
//! This is a metadata projection, not a numerical pickle decoder. Numerical
//! scalars, lists and NumPy buffers are discarded; the independently validated
//! Parquet intermediate remains the source of every scientific value. No
//! callable is executed. Only the NumPy constructors emitted by the pinned
//! producer and the protocol-4 bytearray constructor are recognized as opaque
//! values. Protocols 4 and 5 are supported.
//!
//! A set's serialized element order is needed because loading it constructs a
//! new CPython hash table. Canonical membership cannot recover that order.
//! Apply `python_set_order::StringSet` using an explicit seed after this read.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::rc::Rc;

const MAX_MEMO: usize = 250_000_000;
const MAX_RETAINED_MEMO: usize = 40_000_000;
const MAX_CONTAINERS: usize = 10_000_000;
const MAX_STACK: usize = 2_000_000;
const MAX_TEXT: u64 = 16 * 1024 * 1024;

#[derive(Clone)]
enum Value {
    Ignored,
    List,
    Text(Rc<str>),
    Dict(usize),
    Set(usize),
    Tuple(Rc<Vec<Value>>),
    Global(Global),
    Numpy,
}

#[derive(Clone, Copy)]
enum Global {
    Dtype,
    Reconstruct,
    FromBuffer,
    Scalar,
    Ndarray,
    Bytearray,
}

#[derive(Default)]
pub(crate) struct Statistics {
    pub(crate) bytes: u64,
    pub(crate) skipped_buffer_bytes: u64,
    pub(crate) opcodes: u64,
    pub(crate) memo_entries: usize,
    pub(crate) retained_memo_entries: usize,
    pub(crate) dictionaries: usize,
    pub(crate) sets: usize,
    pub(crate) set_items: usize,
    pub(crate) maximum_stack: usize,
    pub(crate) numpy_reductions: u64,
}

pub(crate) struct ExonSet {
    /// Dictionary path, including the final `exons` key for modern entries.
    pub(crate) path: Vec<String>,
    /// ADDITEMS operands in serialized order, before hash-table reconstruction.
    pub(crate) items: Vec<String>,
}

pub(crate) struct Projection {
    pub(crate) sets: Vec<ExonSet>,
    pub(crate) statistics: Statistics,
}

struct Reader<R> {
    input: BufReader<R>,
    length: u64,
    offset: u64,
    stack: Vec<Value>,
    marks: Vec<usize>,
    // Numeric residue lists dominate this pickle's memo. One tag preserves
    // their identity/type without a full Value allocation for every list.
    memo: Vec<u8>,
    retained_memo: HashMap<usize, Value>,
    dictionaries: Vec<Vec<(Rc<str>, Value)>>,
    sets: Vec<Vec<Rc<str>>>,
    statistics: Statistics,
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

impl<R: Read + Seek> Reader<R> {
    fn error(&self, message: &str) -> io::Error {
        invalid(format!("{message} at pickle byte {}", self.offset))
    }

    fn checked(&self, count: u64) -> io::Result<u64> {
        self.offset
            .checked_add(count)
            .filter(|end| *end <= self.length)
            .ok_or_else(|| self.error("truncated pickle operand"))
    }

    fn bytes<const N: usize>(&mut self) -> io::Result<[u8; N]> {
        let end = self.checked(N as u64)?;
        let mut value = [0; N];
        self.input.read_exact(&mut value)?;
        self.offset = end;
        Ok(value)
    }

    fn byte(&mut self) -> io::Result<u8> {
        Ok(self.bytes::<1>()?[0])
    }
    fn u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.bytes()?))
    }
    fn u64(&mut self) -> io::Result<u64> {
        Ok(u64::from_le_bytes(self.bytes()?))
    }

    fn skip(&mut self, count: u64) -> io::Result<()> {
        let end = self.checked(count)?;
        let relative = i64::try_from(count).map_err(|_| self.error("oversized operand"))?;
        self.input.seek_relative(relative)?;
        self.offset = end;
        Ok(())
    }

    fn text(&mut self, count: u64) -> io::Result<Value> {
        if count > MAX_TEXT {
            return Err(self.error("text exceeds the metadata projection bound"));
        }
        let end = self.checked(count)?;
        let mut bytes = vec![0; count as usize];
        self.input.read_exact(&mut bytes)?;
        self.offset = end;
        let text = String::from_utf8(bytes).map_err(|_| self.error("invalid UTF-8 text"))?;
        Ok(Value::Text(text.into()))
    }

    fn pop(&mut self) -> io::Result<Value> {
        self.stack
            .pop()
            .ok_or_else(|| self.error("operand stack underflow"))
    }

    fn marked(&mut self) -> io::Result<Vec<Value>> {
        let at = self.marks.pop().ok_or_else(|| self.error("missing MARK"))?;
        if at > self.stack.len() {
            return Err(self.error("invalid MARK stack"));
        }
        Ok(self.stack.split_off(at))
    }

    fn target(&self) -> io::Result<&Value> {
        self.stack
            .last()
            .ok_or_else(|| self.error("missing container target"))
    }

    fn assign(&mut self, items: Vec<Value>) -> io::Result<()> {
        let Value::Dict(id) = *self.target()? else {
            return Err(self.error("SETITEM target is not a dictionary"));
        };
        if !items.len().is_multiple_of(2) {
            return Err(self.error("unpaired dictionary items"));
        }
        let mut items = items.into_iter();
        while let Some(key) = items.next() {
            let value = items.next().unwrap();
            // Distribution dictionaries use non-text keys. Retain all text
            // assignments, including overwrites by ignored values, so first
            // key position and last assignment remain Python-compatible.
            if let Value::Text(key) = key {
                self.dictionaries[id].push((key, value));
            } else if matches!(value, Value::Dict(_) | Value::Set(_)) {
                return Err(self.error("non-text metadata dictionary key"));
            }
        }
        Ok(())
    }

    fn memoize(&mut self, index: usize) -> io::Result<()> {
        if index >= MAX_MEMO || index > self.memo.len() {
            return Err(self.error("memo exceeds the metadata projection bound"));
        }
        let value = self.target()?.clone();
        if index == self.memo.len() {
            self.memo.push(0);
        }
        self.retained_memo.remove(&index);
        let tag = match value {
            Value::Ignored => 1,
            Value::List => 2,
            Value::Numpy => 3,
            _ => {
                if self.retained_memo.len() >= MAX_RETAINED_MEMO {
                    return Err(self.error("retained memo exceeds projection bound"));
                }
                self.retained_memo.insert(index, value);
                4
            }
        };
        self.memo[index] = tag;
        Ok(())
    }

    fn run(&mut self) -> io::Result<usize> {
        if self.byte()? != 0x80 || !matches!(self.byte()?, 4 | 5) {
            return Err(self.error("metadata projection requires pickle protocol 4 or 5"));
        }
        loop {
            let opcode = self.byte()?;
            self.statistics.opcodes += 1;
            match opcode {
                b'.' => {
                    if self.offset != self.length || !self.marks.is_empty() || self.stack.len() != 1
                    {
                        return Err(self.error("invalid STOP or trailing bytes"));
                    }
                    let Value::Dict(root) = self.pop()? else {
                        return Err(self.error("root is not a dictionary"));
                    };
                    return Ok(root);
                }
                0x95 => {
                    let count = self.u64()?;
                    self.checked(count)?;
                } // FRAME
                b'(' => self.marks.push(self.stack.len()),
                b'N' | 0x88 | 0x89 => self.stack.push(Value::Ignored),
                b'K' => {
                    self.skip(1)?;
                    self.stack.push(Value::Ignored);
                }
                b'M' => {
                    self.skip(2)?;
                    self.stack.push(Value::Ignored);
                }
                b'J' => {
                    self.skip(4)?;
                    self.stack.push(Value::Ignored);
                }
                b'G' => {
                    self.skip(8)?;
                    self.stack.push(Value::Ignored);
                }
                0x8a => {
                    let count = u64::from(self.byte()?);
                    self.skip(count)?;
                    self.stack.push(Value::Ignored);
                }
                0x8b => {
                    let count = u64::from(self.u32()?);
                    self.skip(count)?;
                    self.stack.push(Value::Ignored);
                }
                0x8c | b'X' | 0x8d => {
                    let count = match opcode {
                        0x8c => u64::from(self.byte()?),
                        b'X' => u64::from(self.u32()?),
                        _ => self.u64()?,
                    };
                    let value = self.text(count)?;
                    self.stack.push(value);
                }
                b'C' | b'B' | 0x8e | 0x96 => {
                    let count = match opcode {
                        b'C' => u64::from(self.byte()?),
                        b'B' => u64::from(self.u32()?),
                        _ => self.u64()?,
                    };
                    self.skip(count)?;
                    self.statistics.skipped_buffer_bytes += count;
                    self.stack.push(Value::Ignored);
                }
                b'}' => {
                    if self.dictionaries.len() >= MAX_CONTAINERS {
                        return Err(self.error("too many dictionaries"));
                    }
                    let id = self.dictionaries.len();
                    self.dictionaries.push(Vec::new());
                    self.stack.push(Value::Dict(id));
                }
                b']' => self.stack.push(Value::List),
                0x8f => {
                    if self.sets.len() >= MAX_CONTAINERS {
                        return Err(self.error("too many sets"));
                    }
                    let id = self.sets.len();
                    self.sets.push(Vec::new());
                    self.stack.push(Value::Set(id));
                }
                b')' => self.stack.push(Value::Tuple(Rc::new(Vec::new()))),
                b't' => {
                    let items = self.marked()?;
                    self.stack.push(Value::Tuple(Rc::new(items)));
                }
                0x85..=0x87 => {
                    let count = usize::from(opcode - 0x84);
                    if self.stack.len() < count {
                        return Err(self.error("tuple stack underflow"));
                    }
                    let items = self.stack.split_off(self.stack.len() - count);
                    self.stack.push(Value::Tuple(Rc::new(items)));
                }
                b's' => {
                    let value = self.pop()?;
                    let key = self.pop()?;
                    self.assign(vec![key, value])?;
                }
                b'u' => {
                    let items = self.marked()?;
                    self.assign(items)?;
                }
                b'a' => {
                    self.pop()?;
                    if !matches!(self.target()?, Value::List) {
                        return Err(self.error("APPEND target is not a list"));
                    }
                }
                b'e' => {
                    self.marked()?;
                    if !matches!(self.target()?, Value::List) {
                        return Err(self.error("APPENDS target is not a list"));
                    }
                }
                0x90 => {
                    let items = self.marked()?;
                    let Value::Set(id) = *self.target()? else {
                        return Err(self.error("ADDITEMS target is not a set"));
                    };
                    for item in items {
                        let Value::Text(text) = item else {
                            return Err(self.error("exon set contains a non-text item"));
                        };
                        self.sets[id].push(text);
                        self.statistics.set_items += 1;
                    }
                }
                0x94 => self.memoize(self.memo.len())?,
                b'q' => {
                    let index = usize::from(self.byte()?);
                    self.memoize(index)?;
                }
                b'r' => {
                    let index = self.u32()? as usize;
                    self.memoize(index)?;
                }
                b'h' | b'j' => {
                    let index = if opcode == b'h' {
                        usize::from(self.byte()?)
                    } else {
                        self.u32()? as usize
                    };
                    let item = match self.memo.get(index) {
                        Some(1) => Value::Ignored,
                        Some(2) => Value::List,
                        Some(3) => Value::Numpy,
                        Some(4) => self.retained_memo.get(&index).unwrap().clone(),
                        _ => return Err(self.error("unset memo slot")),
                    };
                    self.stack.push(item);
                }
                0x93 => {
                    let name = self.pop()?;
                    let module = self.pop()?;
                    let (Value::Text(module), Value::Text(name)) = (module, name) else {
                        return Err(self.error("non-text STACK_GLOBAL"));
                    };
                    let global = match (&*module, &*name) {
                        ("numpy", "dtype") => Global::Dtype,
                        ("numpy", "ndarray") => Global::Ndarray,
                        ("numpy.core.multiarray", "_reconstruct") => Global::Reconstruct,
                        ("numpy.core.multiarray", "scalar") => Global::Scalar,
                        ("numpy.core.numeric", "_frombuffer") => Global::FromBuffer,
                        ("builtins", "bytearray") => Global::Bytearray,
                        _ => return Err(self.error("unsupported global in metadata projection")),
                    };
                    self.stack.push(Value::Global(global));
                }
                b'R' => {
                    let args = self.pop()?;
                    let callable = self.pop()?;
                    let (Value::Global(callable), Value::Tuple(args)) = (callable, args) else {
                        return Err(self.error("unsupported REDUCE operands"));
                    };
                    let count = match callable {
                        Global::Dtype | Global::Reconstruct => 3,
                        Global::Scalar => 2,
                        Global::FromBuffer => 4,
                        Global::Bytearray => 1,
                        Global::Ndarray => {
                            return Err(self.error("direct ndarray REDUCE is unsupported"));
                        }
                    };
                    if args.len() != count {
                        return Err(self.error("unexpected NumPy constructor arity"));
                    }
                    if matches!(callable, Global::Bytearray) {
                        self.stack.push(Value::Ignored);
                    } else {
                        self.statistics.numpy_reductions += 1;
                        self.stack.push(Value::Numpy);
                    }
                }
                b'b' => {
                    self.pop()?;
                    if !matches!(self.target()?, Value::Numpy) {
                        return Err(self.error("unsupported BUILD target"));
                    }
                }
                _ => return Err(self.error(&format!("unsupported opcode 0x{opcode:02x}"))),
            }
            self.statistics.maximum_stack = self.statistics.maximum_stack.max(self.stack.len());
            if self.stack.len() > MAX_STACK || self.marks.len() > 64 {
                return Err(self.error("operand stack exceeds projection bound"));
            }
        }
    }

    fn visit(
        &self,
        id: usize,
        path: &mut Vec<String>,
        active: &mut HashSet<usize>,
        output: &mut Vec<ExonSet>,
    ) -> io::Result<()> {
        if path.len() > 64 || !active.insert(id) {
            return Err(self.error("cyclic or excessively deep metadata dictionaries"));
        }
        let entries = &self.dictionaries[id];
        let mut positions = HashMap::new();
        let mut ordered = Vec::new();
        for (key, value) in entries {
            if let Some(&index) = positions.get(&**key) {
                ordered[index] = (key, value);
            } else {
                positions.insert(&**key, ordered.len());
                ordered.push((key, value));
            }
        }
        for (key, value) in ordered {
            path.push(key.to_string());
            match value {
                Value::Dict(child) => self.visit(*child, path, active, output)?,
                Value::Set(set) => output.push(ExonSet {
                    path: path.clone(),
                    items: self.sets[*set].iter().map(ToString::to_string).collect(),
                }),
                _ => {}
            }
            path.pop();
        }
        active.remove(&id);
        Ok(())
    }
}

/// Extract text-set opcode sequences without allocating numerical payloads.
pub(crate) fn read<R: Read + Seek>(mut input: R) -> io::Result<Projection> {
    let length = input.seek(SeekFrom::End(0))?;
    input.seek(SeekFrom::Start(0))?;
    let mut reader = Reader {
        input: BufReader::with_capacity(1024 * 1024, input),
        length,
        offset: 0,
        stack: Vec::new(),
        marks: Vec::new(),
        memo: Vec::new(),
        retained_memo: HashMap::new(),
        dictionaries: Vec::new(),
        sets: Vec::new(),
        statistics: Statistics::default(),
    };
    let root = reader.run()?;
    let mut sets = Vec::new();
    reader.visit(root, &mut Vec::new(), &mut HashSet::new(), &mut sets)?;
    reader.statistics.bytes = length;
    reader.statistics.memo_entries = reader.memo.len();
    reader.statistics.retained_memo_entries = reader.retained_memo.len();
    reader.statistics.dictionaries = reader.dictionaries.len();
    reader.statistics.sets = reader.sets.len();
    Ok(Projection {
        sets,
        statistics: reader.statistics,
    })
}
