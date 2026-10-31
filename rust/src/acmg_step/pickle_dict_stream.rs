//! A streaming reader for the large pickles of the ACMG step whose top level is
//! a dict: the ClinVar amino-acid and splice dictionaries (217 MB and 1.16 GB),
//! the ClinVar exon allele-frequency cache and the transcript domain map.
//!
//! `crate::pickle_reader` builds the whole object, and Python's own
//! unpickler keeps every object in its memo until the end; on the splice
//! dictionary that is 19.6 GB and 7.6 GB. The step needs only a few facts per
//! transcript, so this reader hands each top-level `(key, value)` pair to a
//! visitor as soon as the pickle has built it, and the visitor keeps what the
//! criteria need. Memory is then one batch of top-level entries (the pickler
//! writes at most 1,000 per SETITEMS) plus what the visitors keep.
//!
//! The memo is what would keep everything alive: protocol 4 memoizes every
//! string and container. The file is therefore read twice. The first pass only
//! notes which memo slots a BINGET reads back, and where each top-level key
//! starts; the second pass stores only those slots (on the ClinVar
//! dictionaries: repeated strings, never a container) and walks the values of
//! the keys the caller does not keep without building them (it still keeps a
//! string there that a later BINGET reads back).
//!
//! Like `pickle_reader`, it accepts the opcodes of plain data written by
//! Python 3 (protocols 2 to 5) and refuses everything else with the opcode and
//! its byte offset: no GLOBAL, REDUCE, BUILD, bytes, sets or numpy objects.
//! Dict keys below the top level must be text (all keys in these files are),
//! and a dict keeps the first position of a key and the last value written to
//! it, as Python's does.
use crate::cli::{Result, invalid};
use flate2::read::MultiGzDecoder;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::rc::Rc;

/// A value below the top-level dict.
#[derive(Clone, Debug)]
pub(crate) enum Item {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(Rc<str>),
    Tuple(Rc<Vec<Item>>),
    List(Rc<RefCell<Vec<Item>>>),
    Dict(Rc<RefCell<TextDict>>),
}

impl Item {
    /// Python's `type(value).__name__`.
    pub(crate) fn type_name(&self) -> &'static str {
        match self {
            Item::None => "NoneType",
            Item::Bool(_) => "bool",
            Item::Int(_) => "int",
            Item::Float(_) => "float",
            Item::Str(_) => "str",
            Item::Tuple(_) => "tuple",
            Item::List(_) => "list",
            Item::Dict(_) => "dict",
        }
    }

    pub(crate) fn as_str(&self) -> Option<&Rc<str>> {
        match self {
            Item::Str(text) => Some(text),
            _ => None,
        }
    }
}

/// A dict with text keys, in insertion order.
#[derive(Debug, Default)]
pub(crate) struct TextDict {
    pub(crate) entries: Vec<(Rc<str>, Item)>,
    /// Positions by key, built once the dict outgrows a linear search.
    index: Option<HashMap<Rc<str>, usize>>,
}

/// Dicts up to this size are searched linearly (the splice records have 11 keys).
const LINEAR_SEARCH_LIMIT: usize = 16;

impl TextDict {
    fn position(&self, key: &str) -> Option<usize> {
        match &self.index {
            Some(index) => index.get(key).copied(),
            None => self.entries.iter().position(|(name, _)| &**name == key),
        }
    }

    fn insert(&mut self, key: Rc<str>, value: Item) {
        if let Some(position) = self.position(&key) {
            self.entries[position].1 = value;
            return;
        }
        if let Some(index) = &mut self.index {
            index.insert(key.clone(), self.entries.len());
        }
        self.entries.push((key, value));
        if self.index.is_none() && self.entries.len() > LINEAR_SEARCH_LIMIT {
            self.index = Some(
                self.entries
                    .iter()
                    .enumerate()
                    .map(|(position, (name, _))| (name.clone(), position))
                    .collect(),
            );
        }
    }

    /// `dict.get(key)`.
    pub(crate) fn get(&self, key: &str) -> Option<&Item> {
        self.position(key).map(|position| &self.entries[position].1)
    }

    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }
}

/// What the stack holds: the top-level dict (never built), a value, or a part of the value
/// of a top-level key the caller does not keep (walked, not built).
#[derive(Clone, Debug)]
enum Slot {
    Root,
    Value(Item),
    Skipped,
}

/// The bytes of the file, decompressed when its name ends in `.gz`, as the step opens these
/// files (`gzip.open(path) if path.endswith(".gz") else open(path, "rb")`).
fn open(path: &str) -> Result<Stream> {
    let file =
        File::open(path).map_err(|error| invalid(format!("FileNotFoundError: {path}: {error}")))?;
    let source: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, file)))
    } else {
        Box::new(file)
    };
    Ok(Stream {
        source,
        chunk: vec![0; 1 << 22],
        start: 0,
        end: 0,
        offset: 0,
        buffer: Vec::new(),
    })
}

/// The opcode stream, read in large chunks so that taking one byte is an index.
struct Stream {
    source: Box<dyn Read>,
    chunk: Vec<u8>,
    /// The unread part of `chunk` is `start..end`.
    start: usize,
    end: usize,
    /// Bytes taken so far (for messages).
    offset: u64,
    /// Holds an item that crosses the end of a chunk.
    buffer: Vec<u8>,
}

fn ended() -> Box<dyn std::error::Error> {
    invalid("the pickle ends before its STOP opcode")
}

impl Stream {
    /// Reads the next chunk; false at the end of the file.
    fn fill(&mut self) -> Result<bool> {
        self.start = 0;
        self.end = 0;
        loop {
            match self.source.read(&mut self.chunk) {
                Ok(count) => {
                    self.end = count;
                    return Ok(count > 0);
                }
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(error) => return Err(error.into()),
            }
        }
    }

    fn byte(&mut self) -> Result<u8> {
        if self.start == self.end && !self.fill()? {
            return Err(ended());
        }
        let byte = self.chunk[self.start];
        self.start += 1;
        self.offset += 1;
        Ok(byte)
    }

    /// The next `count` bytes.
    fn take(&mut self, count: u64) -> Result<&[u8]> {
        let count = usize::try_from(count).map_err(|_| invalid("the pickle holds a huge item"))?;
        self.offset += count as u64;
        if self.end - self.start >= count {
            let slice = &self.chunk[self.start..self.start + count];
            self.start += count;
            return Ok(slice);
        }
        self.buffer.clear();
        while self.buffer.len() < count {
            if self.start == self.end && !self.fill()? {
                return Err(ended());
            }
            let part = (count - self.buffer.len()).min(self.end - self.start);
            self.buffer
                .extend_from_slice(&self.chunk[self.start..self.start + part]);
            self.start += part;
        }
        Ok(&self.buffer)
    }

    /// The next `count` bytes and the byte after them (the next opcode), which is not taken.
    fn take_and_peek(&mut self, count: u64) -> Result<(&[u8], Option<u8>)> {
        let count = usize::try_from(count).map_err(|_| invalid("the pickle holds a huge item"))?;
        self.offset += count as u64;
        if self.end - self.start > count {
            let next = self.chunk[self.start + count];
            let slice = &self.chunk[self.start..self.start + count];
            self.start += count;
            return Ok((slice, Some(next)));
        }
        self.buffer.clear();
        while self.buffer.len() < count {
            if self.start == self.end && !self.fill()? {
                return Err(ended());
            }
            let part = (count - self.buffer.len()).min(self.end - self.start);
            self.buffer
                .extend_from_slice(&self.chunk[self.start..self.start + part]);
            self.start += part;
        }
        let next = if self.start == self.end && !self.fill()? {
            None
        } else {
            Some(self.chunk[self.start])
        };
        Ok((&self.buffer, next))
    }

    /// A little-endian unsigned integer of up to eight bytes.
    fn unsigned(&mut self, width: usize) -> Result<u64> {
        let mut bytes = [0_u8; 8];
        bytes[..width].copy_from_slice(self.take(width as u64)?);
        Ok(u64::from_le_bytes(bytes))
    }

    /// Skips `count` bytes.
    fn skip(&mut self, count: u64) -> Result<()> {
        let mut left =
            usize::try_from(count).map_err(|_| invalid("the pickle holds a huge item"))?;
        self.offset += count;
        while left > 0 {
            if self.start == self.end && !self.fill()? {
                return Err(ended());
            }
            let part = left.min(self.end - self.start);
            self.start += part;
            left -= part;
        }
        Ok(())
    }

    fn text(&mut self, length: u64) -> Result<Rc<str>> {
        let offset = self.offset;
        let bytes = self.take(length)?;
        let text = std::str::from_utf8(bytes).map_err(|_| {
            invalid(format!(
                "the pickle holds a string that is not UTF-8 at byte {offset}"
            ))
        })?;
        Ok(Rc::from(text))
    }

    fn long(&mut self, width: usize) -> Result<i64> {
        let length = self.unsigned(width)?;
        if length == 0 {
            return Ok(0);
        }
        if length > 8 {
            return Err(invalid("the pickle holds an integer wider than 64 bits"));
        }
        let bytes = self.take(length)?;
        let mut value: i64 = if bytes[bytes.len() - 1] & 0x80 != 0 {
            -1
        } else {
            0
        };
        for byte in bytes.iter().rev() {
            value = (value << 8) | i64::from(*byte);
        }
        Ok(value)
    }
}

fn refused(opcode: u8, offset: u64) -> Box<dyn std::error::Error> {
    let name = match opcode {
        b'c' => "GLOBAL",
        0x93 => "STACK_GLOBAL",
        b'R' => "REDUCE",
        b'b' => "BUILD",
        0x81 => "NEWOBJ",
        0x92 => "NEWOBJ_EX",
        b'B' | b'C' | 0x8e | 0x96 => "a bytes opcode",
        0x8f | 0x90 | 0x91 => "a set opcode",
        b'I' | b'L' | b'F' | b'V' | b'S' | b'p' | b'g' => "a text-protocol opcode (protocol 0)",
        _ => "an opcode",
    };
    invalid(format!(
        "the pickle uses opcode 0x{opcode:02x} ({name}) at byte {offset}, which the streaming \
         reader does not accept; it reads plain data only"
    ))
}

/// A multiplicative hash for memo slot numbers: the memo is read back about a hundred million
/// times on the splice dictionary, where the standard hasher costs seconds.
#[derive(Default)]
struct SlotHasher(u64);

impl std::hash::Hasher for SlotHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 = (self.0 ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    fn write_u64(&mut self, number: u64) {
        self.0 = number.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    }
}

type Memo = HashMap<u64, Slot, std::hash::BuildHasherDefault<SlotHasher>>;

/// Which memo slots a BINGET or LONG_BINGET reads back, as a bit set.
struct FetchedSlots(Vec<u64>);

impl FetchedSlots {
    fn contains(&self, slot: u64) -> bool {
        let word = (slot / 64) as usize;
        self.0
            .get(word)
            .is_some_and(|bits| bits & (1 << (slot % 64)) != 0)
    }

    fn insert(&mut self, slot: u64) {
        let word = (slot / 64) as usize;
        if word >= self.0.len() {
            self.0.resize(word + 1, 0);
        }
        self.0[word] |= 1 << (slot % 64);
    }
}

/// What the first pass finds: the memo slots read back, and the byte offsets of the opcodes
/// that push the top-level keys, in order.
struct FirstPass {
    fetched: FetchedSlots,
    key_offsets: Vec<u64>,
}

fn underflow() -> Box<dyn std::error::Error> {
    invalid("the pickle takes a value from an empty stack")
}

/// The first pass: walks the opcodes without building anything, with only the offset of the
/// opcode that pushed each stack item. The memo slots must be given in order (MEMOIZE, or
/// BINPUT/LONG_BINPUT naming the next slot), as Python's pickler writes them.
fn first_pass(path: &str) -> Result<FirstPass> {
    let mut stream = open(path)?;
    let mut fetched = FetchedSlots(Vec::new());
    let mut key_offsets = Vec::new();
    let mut starts: Vec<u64> = Vec::new();
    let mut marks: Vec<usize> = Vec::new();
    let mut root_made = false;
    let mut next_slot: u64 = 0;
    loop {
        let offset = stream.offset;
        let opcode = stream.byte()?;
        match opcode {
            0x80 => {
                stream.byte()?;
            }
            0x95 => {
                stream.unsigned(8)?;
            }
            b'.' => {
                return Ok(FirstPass {
                    fetched,
                    key_offsets,
                });
            }
            b'(' => marks.push(starts.len()),
            b'N' | 0x88 | 0x89 | b']' | b')' => starts.push(offset),
            b'}' => {
                root_made |= starts.is_empty();
                starts.push(offset);
            }
            b'J' | b'K' | b'M' | b'G' => {
                stream.skip(match opcode {
                    b'J' => 4,
                    b'K' => 1,
                    b'M' => 2,
                    _ => 8,
                })?;
                starts.push(offset);
            }
            0x8a | 0x8b => {
                let length = stream.unsigned(if opcode == 0x8a { 1 } else { 4 })?;
                stream.skip(length)?;
                starts.push(offset);
            }
            0x8c | b'X' | 0x8d => {
                let length = stream.unsigned(match opcode {
                    0x8c => 1,
                    b'X' => 4,
                    _ => 8,
                })?;
                stream.skip(length)?;
                starts.push(offset);
            }
            b'a' => {
                starts.pop().ok_or_else(underflow)?;
            }
            b's' => {
                if starts.len() < 3 {
                    return Err(underflow());
                }
                if root_made && starts.len() == 3 {
                    key_offsets.push(starts[1]);
                }
                starts.truncate(starts.len() - 2);
            }
            b'e' | b'u' | b't' => {
                let mark = marks
                    .pop()
                    .ok_or_else(|| invalid("the pickle closes a MARK it never opened"))?;
                if mark > starts.len() {
                    return Err(invalid("the pickle's MARK lies above the stack"));
                }
                if opcode == b'u' && root_made && mark == 1 {
                    key_offsets.extend(starts[mark..].iter().step_by(2));
                }
                starts.truncate(mark);
                if opcode == b't' {
                    starts.push(offset);
                }
            }
            0x85 | 0x86 | 0x87 => {
                let count = usize::from(opcode - 0x84);
                if starts.len() < count {
                    return Err(underflow());
                }
                starts.truncate(starts.len() - count);
                starts.push(offset);
            }
            0x94 => next_slot += 1,
            b'q' | b'r' => {
                let slot = stream.unsigned(if opcode == b'q' { 1 } else { 4 })?;
                if slot != next_slot {
                    return Err(invalid(format!(
                        "the pickle's BINPUT at byte {offset} names memo slot {slot} out of order"
                    )));
                }
                next_slot += 1;
            }
            b'h' | b'j' => {
                let slot = stream.unsigned(if opcode == b'h' { 1 } else { 4 })?;
                if slot >= next_slot {
                    return Err(invalid(format!(
                        "BINGET at byte {offset} names an unset memo slot"
                    )));
                }
                fetched.insert(slot);
                starts.push(offset);
            }
            _ => return Err(refused(opcode, offset)),
        }
    }
}

fn pop(stack: &mut Vec<Slot>) -> Result<Slot> {
    stack
        .pop()
        .ok_or_else(|| invalid("the pickle takes a value from an empty stack"))
}

fn misplaced(slot: &Slot, offset: u64) -> Box<dyn std::error::Error> {
    invalid(match slot {
        Slot::Root => format!("the pickle uses its top-level dict as a value at byte {offset}"),
        _ => format!(
            "at byte {offset} the pickle reads back a part of a top-level entry that was not \
             kept; this reader cannot share it"
        ),
    })
}

fn pop_value(stack: &mut Vec<Slot>, offset: u64) -> Result<Item> {
    match pop(stack)? {
        Slot::Value(item) => Ok(item),
        other => Err(misplaced(&other, offset)),
    }
}

fn marked(stack: &mut Vec<Slot>, marks: &mut Vec<usize>) -> Result<Vec<Slot>> {
    let mark = marks
        .pop()
        .ok_or_else(|| invalid("the pickle closes a MARK it never opened"))?;
    if mark > stack.len() {
        return Err(invalid("the pickle's MARK lies above the stack"));
    }
    Ok(stack.split_off(mark))
}

fn values(slots: Vec<Slot>, offset: u64) -> Result<Vec<Item>> {
    slots
        .into_iter()
        .map(|slot| match slot {
            Slot::Value(item) => Ok(item),
            other => Err(misplaced(&other, offset)),
        })
        .collect()
}

fn text_key(key: Item, offset: u64) -> Result<Rc<str>> {
    match key {
        Item::Str(text) => Ok(text),
        other => Err(invalid(format!(
            "the pickle uses a {} as a dict key at byte {offset}; only text keys are expected",
            other.type_name()
        ))),
    }
}

/// Reads the pickle at `path`, whose top level must be a dict, and calls
/// `visit(key, value)` for each entry whose key `keep` accepts, in the order the
/// pickle sets them (a key set twice is visited twice; Python's dict keeps the
/// last value). The values of the other keys are walked without being built.
/// The top-level keys must be text. Returns `len(dict)`, the number of distinct
/// top-level keys.
pub(crate) fn visit_top_level_dict(
    path: &str,
    keep: impl Fn(&str) -> bool,
    mut visit: impl FnMut(Rc<str>, Item) -> Result<()>,
) -> Result<usize> {
    let FirstPass {
        fetched,
        key_offsets,
    } = first_pass(path)?;
    let mut stream = open(path)?;
    let mut stack: Vec<Slot> = Vec::new();
    let mut marks: Vec<usize> = Vec::new();
    let mut memo = Memo::default();
    let mut next_slot: u64 = 0;
    let mut root_made = false;
    let mut keys: std::collections::HashSet<Rc<str>> = std::collections::HashSet::new();
    let mut next_key = 0_usize;
    // Inside the value of a key that is not kept.
    let mut skipping = false;
    loop {
        let offset = stream.offset;
        let opcode = stream.byte()?;
        let pushes_key = key_offsets.get(next_key) == Some(&offset);
        if pushes_key {
            skipping = false;
        }
        let stack_before = stack.len();
        match opcode {
            0x80 => {
                let version = stream.byte()?;
                if version > 5 {
                    return Err(invalid(format!(
                        "pickle protocol {version} is not supported"
                    )));
                }
            }
            0x95 => {
                stream.unsigned(8)?;
            }
            b'.' => {
                return match (stack.len(), stack.first()) {
                    (1, Some(Slot::Root)) => Ok(keys.len()),
                    _ => Err(invalid("the pickle's top level is not a dict")),
                };
            }
            b'(' => marks.push(stack.len()),
            b'N' | 0x88 | 0x89 if skipping => stack.push(Slot::Skipped),
            b'N' => stack.push(Slot::Value(Item::None)),
            0x88 => stack.push(Slot::Value(Item::Bool(true))),
            0x89 => stack.push(Slot::Value(Item::Bool(false))),
            b'J' => {
                let number = stream.unsigned(4)? as u32 as i32;
                stack.push(Slot::Value(Item::Int(i64::from(number))));
            }
            b'K' => {
                let number = stream.unsigned(1)?;
                stack.push(Slot::Value(Item::Int(number as i64)));
            }
            b'M' => {
                let number = stream.unsigned(2)?;
                stack.push(Slot::Value(Item::Int(number as i64)));
            }
            0x8a => {
                let number = stream.long(1)?;
                stack.push(Slot::Value(Item::Int(number)));
            }
            0x8b => {
                let number = stream.long(4)?;
                stack.push(Slot::Value(Item::Int(number)));
            }
            b'G' => {
                let bits = stream.unsigned(8)?;
                stack.push(Slot::Value(Item::Float(f64::from_bits(bits.swap_bytes()))));
            }
            0x8c | b'X' | 0x8d => {
                let width = match opcode {
                    0x8c => 1,
                    b'X' => 4,
                    _ => 8,
                };
                let length = stream.unsigned(width)?;
                if skipping {
                    // Built only when the next opcode memoizes it into a slot read back later.
                    let (bytes, next) = stream.take_and_peek(length)?;
                    let memoized =
                        matches!(next, Some(0x94 | b'q' | b'r')) && fetched.contains(next_slot);
                    if memoized {
                        let text = std::str::from_utf8(bytes).map_err(|_| {
                            invalid(format!(
                                "the pickle holds a string that is not UTF-8 at byte {offset}"
                            ))
                        })?;
                        stack.push(Slot::Value(Item::Str(Rc::from(text))));
                    } else {
                        stack.push(Slot::Skipped);
                    }
                } else {
                    let text = stream.text(length)?;
                    stack.push(Slot::Value(Item::Str(text)));
                }
            }
            b']' | b'}' | b')' if skipping => stack.push(Slot::Skipped),
            b']' => stack.push(Slot::Value(Item::List(Rc::new(RefCell::new(Vec::new()))))),
            b'a' | b'e' => {
                let items = if opcode == b'a' {
                    vec![pop(&mut stack)?]
                } else {
                    marked(&mut stack, &mut marks)?
                };
                match stack.last() {
                    Some(Slot::Skipped) => {}
                    Some(Slot::Value(Item::List(list))) => {
                        list.borrow_mut().extend(values(items, offset)?)
                    }
                    _ => return Err(invalid(format!("APPEND at byte {offset} needs a list"))),
                }
            }
            b'}' => {
                if !root_made && stack.is_empty() {
                    root_made = true;
                    stack.push(Slot::Root);
                } else {
                    stack.push(Slot::Value(Item::Dict(Rc::new(RefCell::new(
                        TextDict::default(),
                    )))));
                }
            }
            b's' | b'u' => {
                let items = if opcode == b's' {
                    let value = pop(&mut stack)?;
                    let key = pop(&mut stack)?;
                    vec![key, value]
                } else {
                    marked(&mut stack, &mut marks)?
                };
                if items.len() % 2 != 0 {
                    return Err(invalid(format!("SETITEMS at byte {offset} has a lone key")));
                }
                let mut pairs = items.into_iter();
                match stack.last() {
                    Some(Slot::Skipped) => {}
                    Some(Slot::Root) => {
                        while let (Some(key), Some(value)) = (pairs.next(), pairs.next()) {
                            let key = match key {
                                Slot::Value(item) => text_key(item, offset)?,
                                other => return Err(misplaced(&other, offset)),
                            };
                            keys.insert(key.clone());
                            match value {
                                Slot::Skipped => {}
                                Slot::Value(item) => visit(key, item)?,
                                Slot::Root => return Err(misplaced(&value, offset)),
                            }
                        }
                        skipping = false;
                    }
                    Some(Slot::Value(Item::Dict(dict))) => {
                        let mut dict = dict.borrow_mut();
                        let items = values(pairs.collect(), offset)?;
                        let mut pairs = items.into_iter();
                        while let (Some(key), Some(value)) = (pairs.next(), pairs.next()) {
                            dict.insert(text_key(key, offset)?, value);
                        }
                    }
                    _ => return Err(invalid(format!("SETITEM at byte {offset} needs a dict"))),
                }
            }
            b')' => stack.push(Slot::Value(Item::Tuple(Rc::new(Vec::new())))),
            b't' | 0x85 | 0x86 | 0x87 => {
                let items = match opcode {
                    b't' => marked(&mut stack, &mut marks)?,
                    _ => {
                        let count = usize::from(opcode - 0x84);
                        if stack.len() < count {
                            return Err(underflow());
                        }
                        let start = stack.len() - count;
                        stack.split_off(start)
                    }
                };
                if skipping {
                    stack.push(Slot::Skipped);
                } else {
                    stack.push(Slot::Value(Item::Tuple(Rc::new(values(items, offset)?))));
                }
            }
            0x94 | b'q' | b'r' => {
                if opcode != 0x94 {
                    stream.unsigned(if opcode == b'q' { 1 } else { 4 })?;
                }
                let slot = next_slot;
                next_slot += 1;
                if fetched.contains(slot) {
                    let top = stack
                        .last()
                        .cloned()
                        .ok_or_else(|| invalid("MEMOIZE with an empty stack"))?;
                    if matches!(top, Slot::Root) {
                        return Err(invalid("the pickle reads its top-level dict back"));
                    }
                    memo.insert(slot, top);
                }
            }
            b'h' | b'j' => {
                let slot = stream.unsigned(if opcode == b'h' { 1 } else { 4 })?;
                let item = memo.get(&slot).cloned().ok_or_else(|| {
                    invalid(format!("BINGET at byte {offset} names an unset memo slot"))
                })?;
                if matches!(item, Slot::Skipped) && !skipping {
                    return Err(misplaced(&item, offset));
                }
                stack.push(item);
            }
            _ => return Err(refused(opcode, offset)),
        }
        if pushes_key {
            next_key += 1;
            match stack.last() {
                Some(Slot::Value(Item::Str(key))) if stack.len() == stack_before + 1 => {
                    skipping = !keep(key);
                }
                // Not text: SETITEMS refuses it with its type.
                _ => skipping = false,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn temporary(name: &str, bytes: &[u8]) -> String {
        let path = std::env::temp_dir().join(format!(
            "priva_pickle_dict_stream_{}_{name}",
            std::process::id()
        ));
        File::create(&path).unwrap().write_all(bytes).unwrap();
        path.to_str().unwrap().to_owned()
    }

    fn collect(path: &str) -> Result<Vec<(String, String)>> {
        collect_kept(path, |_| true).map(|(seen, _)| seen)
    }

    fn collect_kept(
        path: &str,
        keep: impl Fn(&str) -> bool,
    ) -> Result<(Vec<(String, String)>, usize)> {
        let mut seen = Vec::new();
        let count = visit_top_level_dict(path, keep, |key, value| {
            seen.push((key.to_string(), describe(&value)));
            Ok(())
        })?;
        Ok((seen, count))
    }

    fn describe(item: &Item) -> String {
        match item {
            Item::None => "None".into(),
            Item::Bool(flag) => format!("{flag}"),
            Item::Int(number) => format!("{number}"),
            Item::Float(number) => format!("{number:?}"),
            Item::Str(text) => format!("'{text}'"),
            Item::Tuple(items) => format!(
                "({})",
                items.iter().map(describe).collect::<Vec<_>>().join(", ")
            ),
            Item::List(items) => format!(
                "[{}]",
                items
                    .borrow()
                    .iter()
                    .map(describe)
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            Item::Dict(dict) => format!(
                "{{{}}}",
                dict.borrow()
                    .entries
                    .iter()
                    .map(|(key, value)| format!("'{key}': {}", describe(value)))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        }
    }

    #[test]
    fn reads_python_pickles_with_shared_strings() {
        // pickle.dumps({"T1": {"a": ["x", "x"], "b": (1, 2.5, None, True)},
        //               "T2": [{"a": "x"}]}, protocol=4), with "x" and "a" shared.
        let bytes = b"\x80\x04\x95?\x00\x00\x00\x00\x00\x00\x00}\x94(\x8c\x02T1\x94}\x94(\x8c\x01a\x94]\x94(\x8c\x01x\x94h\x05e\x8c\x01b\x94(K\x01G@\x04\x00\x00\x00\x00\x00\x00N\x88t\x94u\x8c\x02T2\x94]\x94}\x94h\x03h\x05sau.";
        let seen = collect(&temporary("protocol4", bytes)).unwrap();
        assert_eq!(
            seen,
            vec![
                (
                    "T1".to_owned(),
                    "{'a': ['x', 'x'], 'b': (1, 2.5, None, true)}".to_owned()
                ),
                ("T2".to_owned(), "[{'a': 'x'}]".to_owned()),
            ]
        );
        // The same shapes with protocol 2 (BINPUT, BININT, LONG1).
        let bytes = b"\x80\x02}q\x00(X\x02\x00\x00\x00T1q\x01}q\x02X\x01\x00\x00\x00aq\x03]q\x04(X\x01\x00\x00\x00xq\x05h\x05esX\x02\x00\x00\x00T2q\x06(J\xfd\xff\xff\xffM,\x01Jp\x11\x01\x00\x8a\x06\x00\x00\x00\x00\x00\x01\x8a\x06\x00\x00\x00\x00\x00\xfftq\x07u.";
        let seen = collect(&temporary("protocol2", bytes)).unwrap();
        assert_eq!(
            seen,
            vec![
                ("T1".to_owned(), "{'a': ['x', 'x']}".to_owned()),
                (
                    "T2".to_owned(),
                    "(-3, 300, 70000, 1099511627776, -1099511627776)".to_owned()
                ),
            ]
        );
    }

    #[test]
    fn skips_values_but_keeps_strings_read_back_later() {
        // The protocol 4 pickle above: "a" and "x" are first written inside T1's value and read
        // back (BINGET) in T2's.
        let bytes = b"\x80\x04\x95?\x00\x00\x00\x00\x00\x00\x00}\x94(\x8c\x02T1\x94}\x94(\x8c\x01a\x94]\x94(\x8c\x01x\x94h\x05e\x8c\x01b\x94(K\x01G@\x04\x00\x00\x00\x00\x00\x00N\x88t\x94u\x8c\x02T2\x94]\x94}\x94h\x03h\x05sau.";
        let path = temporary("skipping", bytes);
        let (seen, count) = collect_kept(&path, |key| key == "T2").unwrap();
        assert_eq!(seen, vec![("T2".to_owned(), "[{'a': 'x'}]".to_owned())]);
        assert_eq!(count, 2);
        let (seen, count) = collect_kept(&path, |_| false).unwrap();
        assert!(seen.is_empty());
        assert_eq!(count, 2);
    }

    #[test]
    fn refuses_globals_and_other_tops() {
        // pickle.dumps({"T": complex(1, 2)}, protocol=4) holds a STACK_GLOBAL.
        let bytes = b"\x80\x04\x955\x00\x00\x00\x00\x00\x00\x00}\x94\x8c\x01T\x94\x8c\x08builtins\x94\x8c\x07complex\x94\x93\x94G?\xf0\x00\x00\x00\x00\x00\x00G@\x00\x00\x00\x00\x00\x00\x00\x86\x94R\x94s.";
        let error = collect(&temporary("global", bytes)).unwrap_err();
        assert!(error.to_string().contains("STACK_GLOBAL"), "{error}");
        // pickle.dumps(["T"], protocol=4): a list at the top.
        let bytes = b"\x80\x04\x95\x08\x00\x00\x00\x00\x00\x00\x00]\x94\x8c\x01T\x94a.";
        let error = collect(&temporary("list", bytes)).unwrap_err();
        assert!(error.to_string().contains("not a dict"), "{error}");
        // pickle.dumps({"T": {1: 2}}, protocol=4): a number as a key below the top.
        let bytes =
            b"\x80\x04\x95\x0f\x00\x00\x00\x00\x00\x00\x00}\x94\x8c\x01T\x94}\x94K\x01K\x02ss.";
        let error = collect(&temporary("numberkey", bytes)).unwrap_err();
        assert!(error.to_string().contains("as a dict key"), "{error}");
    }

    #[test]
    fn a_repeated_key_keeps_its_place_and_takes_the_last_value() {
        let mut dict = TextDict::default();
        for number in 0..40 {
            dict.insert(
                Rc::from(format!("k{}", number % 20).as_str()),
                Item::Int(number),
            );
        }
        assert_eq!(dict.len(), 20);
        assert!(matches!(dict.get("k3"), Some(Item::Int(23))));
        assert_eq!(&*dict.entries[3].0, "k3");
        assert!(dict.get("k20").is_none());
    }
}
