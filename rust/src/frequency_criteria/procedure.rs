//! Native, source-ordered frequency procedures with explicit host effects.
//!
//! The host owns existing Python objects and mature package calls. It does not
//! choose a frequency branch, the next getter, an arithmetic operation, or the
//! order of publication. A reply advances exactly one outstanding effect.
//! Nothing here executes a Python procedure body or replays a completed prefix.
//!
//! This first complete procedure is `gnomAD_rare_AF` (source lines 256--269).
//! The other frequency entry points are deliberately not advertised by this
//! module until their native control programs exist. This module is currently
//! unregistered; the existing column-operation bridge remains a separate API.

use super::{calculate, Value};

pub const SOURCE_SHA256: &str =
    "ef2e70115ddfe7087b4d24002ce4af95c4f4b91fbf38a06e8a52275eecca4fb3";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Procedure {
    GnomadRareAf,
}

impl Procedure {
    pub fn named(name: &str) -> Result<Self, String> {
        match name {
            "gnomAD_rare_AF" => Ok(Self::GnomadRareAf),
            _ => Err("frequency procedure has no native control program".into()),
        }
    }

    pub fn name(self) -> &'static str {
        match self { Self::GnomadRareAf => "gnomAD_rare_AF" }
    }
}

/// These retain the existing bridge's explicit execution modes. Mixed mode
/// permits only the requested single mature operation, never a function retry.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mode { NativeMath, Mixed }

/// Object handles refer to the owning host's private registry. Registering an
/// argument must not inspect a frame, decode a cache, or admit its later fields.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Object {
    Handle(u64),
    None,
    Integer(i64),
    Float(u64),
    Text(&'static str),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Site { pub procedure: Procedure, pub line: u16, pub operation: u16 }

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TypeTest { Float, String }

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kernel { FillNa, Divide, GreaterEqual, LessEqual, Where }

impl Kernel {
    pub fn name(self) -> &'static str {
        match self {
            Self::FillNa => "fillna", Self::Divide => "div",
            Self::GreaterEqual => "ge", Self::LessEqual => "le", Self::Where => "where",
        }
    }
}

/// The original callable is captured before its arguments. An operator has no
/// mutable module callable to capture. Its exact operand objects are retained.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Invocation {
    Bound { callable: Object, args: Vec<Object> },
    Binary { operation: Kernel, left: Object, right: Object },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KernelCall {
    pub kernel: Kernel,
    pub invocation: Invocation,
    pub operands: Vec<Object>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Effect {
    IsInstance { object: Object, class: TypeTest },
    BindNumpyWhere,
    GetItem { owner: Object, key: Object },
    GetAttribute { owner: Object, name: &'static str },
    /// The host checks the captured callable and nonwarming operand transport
    /// at this source position. It returns typed values without doing math.
    PrepareKernel(KernelCall),
    /// The same typed Polars kernels used by the column API have already run
    /// in this Rust process. Host boxing preserves source axis/dtype/ownership.
    Materialize { call: KernelCall, value: Value },
    /// Inspect only exact built-in float identity/bits. A subclass or NumPy
    /// scalar is not relabeled as a Python float by a conversion callback.
    InspectBuiltinFloat(Object),
    InvokeMature { invocation: Invocation, reason: String },
}

#[derive(Clone, Debug, PartialEq)]
pub enum Action {
    Effect { token: u64, site: Site, effect: Effect },
    Return(Object),
    Raise { class: &'static str, message: &'static str },
    /// Preserve the actual host exception object, including args and cause.
    HostRaised { exception: u64 },
    Unsupported { site: Site, reason: String },
}

#[derive(Clone, Debug, PartialEq)]
pub enum Reply {
    Truth(bool),
    Object(u64),
    Prepared(Vec<Value>),
    NotRepresented(String),
    BuiltinFloat(u64),
    OtherScalar,
    Raised(u64),
}

#[derive(Clone, Copy)]
enum Arg { Slot(usize), Int(i64), Float(u64), Text(&'static str) }

#[derive(Clone, Copy)]
enum Op {
    Where(usize),
    Item { owner: usize, key: Arg, dest: usize },
    Attribute { owner: usize, name: &'static str, dest: usize },
    Fill { owner: usize, callable: usize, value: Arg, dest: usize },
    Reciprocal { operand: usize, dest: usize },
    Binary { kernel: Kernel, left: Arg, right: Arg, dest: usize },
    Select { callable: usize, condition: usize, yes: usize, no: usize, dest: usize },
    Return(usize),
}

#[derive(Clone, Copy)]
struct Instruction { line: u16, op: Op }

const fn at(line: u16, op: Op) -> Instruction { Instruction { line, op } }

// Slots 0 and 1 are df and cutoff. Every repeated source getter remains a
// separate instruction. np.where is bound before ANY of its eager arguments.
const FLOAT: &[Instruction] = &[
    at(263, Op::Where(2)),
    at(263, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AN_max"), dest: 3 }),
    at(263, Op::Attribute { owner: 3, name: "fillna", dest: 4 }),
    at(263, Op::Fill { owner: 3, callable: 4, value: Arg::Int(1_000_000), dest: 5 }),
    at(263, Op::Reciprocal { operand: 1, dest: 6 }),
    at(263, Op::Binary { kernel: Kernel::GreaterEqual, left: Arg::Slot(5), right: Arg::Slot(6), dest: 7 }),
    at(264, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AF_max"), dest: 8 }),
    at(264, Op::Attribute { owner: 8, name: "fillna", dest: 9 }),
    at(264, Op::Fill { owner: 8, callable: 9, value: Arg::Int(0), dest: 10 }),
    at(264, Op::Binary { kernel: Kernel::LessEqual, left: Arg::Slot(10), right: Arg::Slot(1), dest: 11 }),
    at(265, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AF"), dest: 12 }),
    at(265, Op::Attribute { owner: 12, name: "fillna", dest: 13 }),
    at(265, Op::Fill { owner: 12, callable: 13, value: Arg::Int(0), dest: 14 }),
    at(265, Op::Binary { kernel: Kernel::LessEqual, left: Arg::Slot(14), right: Arg::Slot(1), dest: 15 }),
    at(263, Op::Select { callable: 2, condition: 7, yes: 11, no: 15, dest: 16 }),
    at(263, Op::Return(16)),
];

const STRING: &[Instruction] = &[
    at(267, Op::Where(2)),
    at(267, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AN_max"), dest: 3 }),
    at(267, Op::Attribute { owner: 3, name: "fillna", dest: 4 }),
    at(267, Op::Fill { owner: 3, callable: 4, value: Arg::Int(1_000_000), dest: 5 }),
    at(267, Op::Item { owner: 0, key: Arg::Slot(1), dest: 6 }),
    at(267, Op::Attribute { owner: 6, name: "fillna", dest: 7 }),
    // The source spells 10e-6, a float, rather than an integer zero default.
    at(267, Op::Fill { owner: 6, callable: 7, value: Arg::Float(0x3ee4_f8b5_88e3_68f1), dest: 8 }),
    at(267, Op::Binary { kernel: Kernel::Divide, left: Arg::Int(1), right: Arg::Slot(8), dest: 9 }),
    at(267, Op::Binary { kernel: Kernel::GreaterEqual, left: Arg::Slot(5), right: Arg::Slot(9), dest: 10 }),
    at(268, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AF_max"), dest: 11 }),
    at(268, Op::Attribute { owner: 11, name: "fillna", dest: 12 }),
    at(268, Op::Fill { owner: 11, callable: 12, value: Arg::Int(0), dest: 13 }),
    at(268, Op::Item { owner: 0, key: Arg::Slot(1), dest: 14 }),
    at(268, Op::Attribute { owner: 14, name: "fillna", dest: 15 }),
    at(268, Op::Fill { owner: 14, callable: 15, value: Arg::Int(0), dest: 16 }),
    at(268, Op::Binary { kernel: Kernel::LessEqual, left: Arg::Slot(13), right: Arg::Slot(16), dest: 17 }),
    at(269, Op::Item { owner: 0, key: Arg::Text("gnomAD_joint_AF"), dest: 18 }),
    at(269, Op::Attribute { owner: 18, name: "fillna", dest: 19 }),
    at(269, Op::Fill { owner: 18, callable: 19, value: Arg::Int(0), dest: 20 }),
    at(269, Op::Item { owner: 0, key: Arg::Slot(1), dest: 21 }),
    at(269, Op::Attribute { owner: 21, name: "fillna", dest: 22 }),
    at(269, Op::Fill { owner: 21, callable: 22, value: Arg::Int(0), dest: 23 }),
    at(269, Op::Binary { kernel: Kernel::LessEqual, left: Arg::Slot(20), right: Arg::Slot(23), dest: 24 }),
    at(267, Op::Select { callable: 2, condition: 10, yes: 17, no: 24, dest: 25 }),
    at(267, Op::Return(25)),
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Phase { Ready, FloatCheck, StringCheck, Float, String, Finished, Poisoned }

#[derive(Clone, Debug)]
enum Continuation {
    FloatCheck,
    StringCheck,
    Store(usize),
    Kernel { call: KernelCall, dest: usize },
    Reciprocal { operand: Object, dest: usize },
}

#[derive(Clone, Debug)]
struct Pending { token: u64, site: Site, continuation: Continuation }

/// Only this native state determines the next operation. Python argument
/// registration is identity-only, and no later input is inspected up front.
pub struct Program {
    procedure: Procedure,
    mode: Mode,
    phase: Phase,
    pc: usize,
    slots: Vec<Option<Object>>,
    pending: Option<Pending>,
    next_token: u64,
    pub native_operations: u64,
    pub mature_operations: u64,
}

impl Program {
    pub fn new(procedure: Procedure, mode: Mode, df: u64, cutoff: u64) -> Result<Self, String> {
        if df == 0 || cutoff == 0 { return Err("zero host object handle".into()); }
        let mut slots = vec![None; 26];
        slots[0] = Some(Object::Handle(df));
        slots[1] = Some(Object::Handle(cutoff));
        Ok(Self { procedure, mode, phase: Phase::Ready, pc: 0, slots,
            pending: None, next_token: 0, native_operations: 0, mature_operations: 0 })
    }

    fn slot(&self, slot: usize) -> Result<Object, String> {
        self.slots.get(slot).and_then(Clone::clone).ok_or_else(|| "unassigned native procedure slot".into())
    }

    fn arg(&self, arg: Arg) -> Result<Object, String> {
        match arg { Arg::Slot(slot) => self.slot(slot), Arg::Int(value) => Ok(Object::Integer(value)),
            Arg::Float(bits) => Ok(Object::Float(bits)), Arg::Text(value) => Ok(Object::Text(value)) }
    }

    fn issue(&mut self, site: Site, effect: Effect, continuation: Continuation) -> Result<Action, String> {
        if self.pending.is_some() { return Err("native procedure already awaits an effect".into()); }
        let token = self.next_token;
        self.next_token = token.checked_add(1).ok_or("effect token overflow")?;
        self.pending = Some(Pending { token, site, continuation });
        Ok(Action::Effect { token, site, effect })
    }

    pub fn start(&mut self) -> Result<Action, String> {
        if self.phase != Phase::Ready { self.phase = Phase::Poisoned; return Err("procedure already started".into()); }
        self.phase = Phase::FloatCheck;
        self.issue(Site { procedure: self.procedure, line: 262, operation: 0 },
            Effect::IsInstance { object: self.slot(1)?, class: TypeTest::Float }, Continuation::FloatCheck)
    }

    fn kernel(&mut self, site: Site, call: KernelCall, dest: usize) -> Result<Action, String> {
        self.issue(site, Effect::PrepareKernel(call.clone()), Continuation::Kernel { call, dest })
    }

    fn drive(&mut self) -> Result<Action, String> {
        let tape = match self.phase { Phase::Float => FLOAT, Phase::String => STRING,
            _ => return Err("native procedure is not at a running instruction".into()) };
        let instruction = *tape.get(self.pc).ok_or("native procedure counter outside body")?;
        let site = Site { procedure: self.procedure, line: instruction.line, operation: self.pc as u16 + 2 };
        match instruction.op {
            Op::Where(dest) => self.issue(site, Effect::BindNumpyWhere, Continuation::Store(dest)),
            Op::Item { owner, key, dest } => self.issue(site,
                Effect::GetItem { owner: self.slot(owner)?, key: self.arg(key)? }, Continuation::Store(dest)),
            Op::Attribute { owner, name, dest } => self.issue(site,
                Effect::GetAttribute { owner: self.slot(owner)?, name }, Continuation::Store(dest)),
            Op::Fill { owner, callable, value, dest } => {
                let fill = self.arg(value)?;
                self.kernel(site, KernelCall { kernel: Kernel::FillNa,
                    invocation: Invocation::Bound { callable: self.slot(callable)?, args: vec![fill.clone()] },
                    operands: vec![self.slot(owner)?, fill] }, dest)
            }
            Op::Reciprocal { operand, dest } => {
                let operand = self.slot(operand)?;
                self.issue(site, Effect::InspectBuiltinFloat(operand.clone()), Continuation::Reciprocal { operand, dest })
            }
            Op::Binary { kernel, left, right, dest } => {
                let left = self.arg(left)?; let right = self.arg(right)?;
                self.kernel(site, KernelCall { kernel,
                    invocation: Invocation::Binary { operation: kernel, left: left.clone(), right: right.clone() },
                    operands: vec![left, right] }, dest)
            }
            Op::Select { callable, condition, yes, no, dest } => {
                let args = vec![self.slot(condition)?, self.slot(yes)?, self.slot(no)?];
                self.kernel(site, KernelCall { kernel: Kernel::Where,
                    invocation: Invocation::Bound { callable: self.slot(callable)?, args: args.clone() }, operands: args }, dest)
            }
            Op::Return(slot) => { let value = self.slot(slot)?; self.phase = Phase::Finished; Ok(Action::Return(value)) }
        }
    }

    fn store(&mut self, dest: usize, value: Object) -> Result<Action, String> {
        *self.slots.get_mut(dest).ok_or("invalid native destination")? = Some(value);
        self.pc = self.pc.checked_add(1).ok_or("native instruction overflow")?;
        self.drive()
    }

    fn mature(&mut self, site: Site, invocation: Invocation, dest: usize, reason: String) -> Result<Action, String> {
        if self.mode == Mode::NativeMath {
            self.phase = Phase::Finished;
            return Ok(Action::Unsupported { site, reason });
        }
        self.mature_operations += 1;
        self.issue(site, Effect::InvokeMature { invocation, reason }, Continuation::Store(dest))
    }

    fn accept(&mut self, pending: Pending, reply: Reply) -> Result<Action, String> {
        if let Reply::Raised(exception) = reply {
            if exception == 0 { return Err("zero host exception handle".into()); }
            self.phase = Phase::Finished;
            return Ok(Action::HostRaised { exception });
        }
        match (pending.continuation, reply) {
            (Continuation::FloatCheck, Reply::Truth(true)) => { self.phase = Phase::Float; self.drive() }
            (Continuation::FloatCheck, Reply::Truth(false)) => {
                self.phase = Phase::StringCheck;
                self.issue(Site { procedure: self.procedure, line: 266, operation: 1 },
                    Effect::IsInstance { object: self.slot(1)?, class: TypeTest::String }, Continuation::StringCheck)
            }
            (Continuation::StringCheck, Reply::Truth(true)) => { self.phase = Phase::String; self.drive() }
            (Continuation::StringCheck, Reply::Truth(false)) => { self.phase = Phase::Finished; Ok(Action::Return(Object::None)) }
            (Continuation::Store(dest), Reply::Object(handle)) if handle != 0 => self.store(dest, Object::Handle(handle)),
            (Continuation::Kernel { call, dest }, Reply::Prepared(values)) => {
                if values.len() != call.operands.len() { return Err("native operand reply count differs".into()); }
                if values.iter().any(|value| {
                    !matches!(value.kind.as_str(), "scalar" | "series" | "array" | "list")
                        || value.kind == "scalar" && value.cells.len() != 1
                        || value.kind == "series" && value.index.len() != value.cells.len()
                        || value.kind != "series" && !value.index.is_empty()
                }) { return Err("malformed typed operand shape".into()); }
                match calculate(call.kernel.name(), &values) {
                    Ok(value) => {
                        self.native_operations += 1;
                        self.issue(pending.site, Effect::Materialize { call, value }, Continuation::Store(dest))
                    }
                    Err(reason) => self.mature(pending.site, call.invocation, dest, reason),
                }
            }
            (Continuation::Kernel { call, dest }, Reply::NotRepresented(reason)) if !reason.is_empty() =>
                self.mature(pending.site, call.invocation, dest, reason),
            (Continuation::Reciprocal { operand, dest }, Reply::BuiltinFloat(bits)) => {
                let number = f64::from_bits(bits);
                if number == 0.0 {
                    self.phase = Phase::Finished;
                    return Ok(Action::Raise { class: "ZeroDivisionError", message: "float division by zero" });
                }
                let result = if number.is_finite() { Some(1.0 / number) } else { None };
                if let Some(result) = result.filter(|value| value.is_finite()) {
                    self.native_operations += 1;
                    self.store(dest, Object::Float(result.to_bits()))
                } else {
                    self.mature(pending.site, Invocation::Binary { operation: Kernel::Divide,
                        left: Object::Integer(1), right: operand }, dest,
                        "nonfinite Python float division retains its exact mature scalar effect".into())
                }
            }
            (Continuation::Reciprocal { operand, dest }, Reply::OtherScalar) => self.mature(pending.site,
                Invocation::Binary { operation: Kernel::Divide, left: Object::Integer(1), right: operand }, dest,
                "float subclass or NumPy scalar retains its exact single scalar effect".into()),
            _ => Err("wrong reply for outstanding native procedure effect".into()),
        }
    }

    pub fn resume(&mut self, token: u64, reply: Reply) -> Result<Action, String> {
        if matches!(self.phase, Phase::Finished | Phase::Poisoned | Phase::Ready) {
            self.phase = Phase::Poisoned;
            return Err("native procedure is not awaiting a reply".into());
        }
        let Some(pending) = self.pending.take() else {
            self.phase = Phase::Poisoned;
            return Err("no outstanding native procedure effect".into());
        };
        if pending.token != token {
            self.phase = Phase::Poisoned;
            return Err("native procedure reply token differs".into());
        }
        let result = self.accept(pending, reply);
        if result.is_err() { self.phase = Phase::Poisoned; self.pending = None; }
        result
    }
}

#[cfg(test)]
#[path = "procedure_tests.rs"]
mod tests;
