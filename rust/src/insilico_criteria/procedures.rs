//! Named PriVA procedure control, pinned to original b8a093a9.
//!
//! 27 October 2026 CT: the user requires complete PriVA-owned functions in
//! named Rust modules. The older column session is a transitional control:
//! Python still orders its branches. Here Rust owns sequencing, constants,
//! masks, branches and publication. The host only services visited mature
//! effects/ownership and one actual unsupported operation, never the original
//! criterion or a completed-prefix replay. No input observation is cached.
use super::{Cell, Value, calculate};

pub const SOURCE_SHA256: &str = "b8a093a97d0ad9c4730540b2a57a60a302794bf84a42ba304a22e782053a81a2";
pub type Object = u64;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Site {
    pub function: &'static str,
    pub line: u32,
}

#[derive(Clone, Debug)]
pub enum Target {
    Binary(&'static str),
    Invert,
    Bound {
        function: Object,
        receiver: Option<Object>,
    },
    Assignment,
}

#[derive(Clone, Debug)]
pub struct Call {
    pub site: Site,
    pub operation: &'static str,
    pub target: Target,
    pub args: Vec<Object>,
    pub keywords: Vec<(&'static str, Object)>,
}

/// A mature object handle retains its real identity, callbacks and exceptions.
/// `prepare` observes ONLY the already-evaluated operands at this operation.
/// It must not read other frame columns, speculate or cache observations across
/// callbacks. `None` means this ONE operation needs its actual mature binding.
/// `commit_native` only reconstructs/updates ownership; it computes no rule.
/// The concrete transport must distinguish scientific exceptions from corrupt
/// or dead-child faults and must not retry either a call or a written prefix.
pub trait Host {
    type Error;
    fn literal(&mut self, site: Site, value: Cell) -> Result<Object, Self::Error>;
    fn column(
        &mut self,
        site: Site,
        frame: Object,
        key: &'static str,
    ) -> Result<Object, Self::Error>;
    fn attribute(
        &mut self,
        site: Site,
        object: Object,
        name: &'static str,
    ) -> Result<Object, Self::Error>;
    fn global(
        &mut self,
        site: Site,
        module: Option<&'static str>,
        name: &'static str,
    ) -> Result<Object, Self::Error>;
    fn length(&mut self, site: Site, object: Object) -> Result<Object, Self::Error>;
    /// Small mature language/ownership effects. Action and formatting pieces
    /// are chosen by Rust; this is never an expression/function-body evaluator.
    fn effect(
        &mut self,
        site: Site,
        action: &'static str,
        args: &[Object],
        text: &[&str],
    ) -> Result<Object, Self::Error>;
    fn truth(&mut self, site: Site, object: Object) -> Result<bool, Self::Error>;
    fn is_none(&mut self, site: Site, object: Object) -> Result<bool, Self::Error>;
    fn release(&mut self, site: Site, objects: &[Object]) -> Result<(), Self::Error>;
    fn prepare(&mut self, call: &Call) -> Result<Option<Vec<Value>>, Self::Error>;
    fn commit_native(&mut self, call: &Call, value: Value) -> Result<Object, Self::Error>;
    fn mature(&mut self, call: &Call) -> Result<Object, Self::Error>;
}

pub struct Engine<'a, H: Host> {
    host: &'a mut H,
}

impl<'a, H: Host> Engine<'a, H> {
    pub fn new(host: &'a mut H) -> Self {
        Self { host }
    }

    fn execute(&mut self, call: Call) -> Result<Object, H::Error> {
        if let Some(values) = self.host.prepare(&call)? {
            if let Ok(value) = calculate(call.operation, &values) {
                return self.host.commit_native(&call, value);
            }
        }
        // No mature effects have happened in the pure typed-kernel attempt.
        // An unrepresented dtype/shape executes its actual operation once,
        // with the same already-evaluated objects and source position.
        self.host.mature(&call)
    }

    fn binary(
        &mut self,
        site: Site,
        op: &'static str,
        left: Object,
        right: Object,
    ) -> Result<Object, H::Error> {
        self.execute(Call {
            site,
            operation: op,
            target: Target::Binary(op),
            args: vec![left, right],
            keywords: vec![],
        })
    }

    fn method0(
        &mut self,
        site: Site,
        receiver: Object,
        name: &'static str,
    ) -> Result<Object, H::Error> {
        let function = self.host.attribute(site, receiver, name)?;
        let result = self.execute(Call {
            site,
            operation: name,
            target: Target::Bound {
                function,
                receiver: Some(receiver),
            },
            args: vec![receiver],
            keywords: vec![],
        })?;
        self.host.release(site, &[function])?;
        Ok(result)
    }

    fn bound(
        &mut self,
        site: Site,
        operation: &'static str,
        function: Object,
        receiver: Option<Object>,
        args: Vec<Object>,
        keywords: Vec<(&'static str, Object)>,
    ) -> Result<Object, H::Error> {
        let result = self.execute(Call {
            site,
            operation,
            target: Target::Bound { function, receiver },
            args,
            keywords,
        })?;
        self.host.release(site, &[function])?;
        Ok(result)
    }

    fn invert(&mut self, site: Site, value: Object) -> Result<Object, H::Error> {
        self.execute(Call {
            site,
            operation: "invert",
            target: Target::Invert,
            args: vec![value],
            keywords: vec![],
        })
    }

    fn fillna(&mut self, site: Site, receiver: Object, value: Cell) -> Result<Object, H::Error> {
        let function = self.host.attribute(site, receiver, "fillna")?;
        let value = self.host.literal(site, value)?;
        self.bound(
            site,
            "fillna",
            function,
            Some(receiver),
            vec![receiver, value],
            vec![],
        )
    }

    fn as_str(&mut self, site: Site, receiver: Object) -> Result<Object, H::Error> {
        let function = self.host.attribute(site, receiver, "astype")?;
        let dtype = self.host.global(site, None, "str")?;
        self.bound(
            site,
            "as_str",
            function,
            Some(receiver),
            vec![receiver, dtype],
            vec![],
        )
    }

    fn contains(
        &mut self,
        site: Site,
        receiver: Object,
        pattern: &'static str,
        explicit_regex: bool,
    ) -> Result<Object, H::Error> {
        let accessor = self.host.attribute(site, receiver, "str")?;
        let function = self.host.attribute(site, accessor, "contains")?;
        let pattern = self.host.literal(site, Cell::Text(pattern.into()))?;
        let keywords = if explicit_regex {
            vec![("regex", self.host.literal(site, Cell::Bool(true))?)]
        } else {
            vec![]
        };
        self.bound(
            site,
            "contains",
            function,
            Some(accessor),
            vec![accessor, pattern],
            keywords,
        )
    }

    fn strip(&mut self, site: Site, receiver: Object) -> Result<Object, H::Error> {
        let accessor = self.host.attribute(site, receiver, "str")?;
        let function = self.host.attribute(site, accessor, "strip")?;
        self.bound(
            site,
            "strip",
            function,
            Some(accessor),
            vec![accessor],
            vec![],
        )
    }

    fn numeric(&mut self, site: Site, receiver: Object) -> Result<Object, H::Error> {
        let function = self.host.global(site, Some("pd"), "to_numeric")?;
        self.numeric_bound(site, function, receiver)
    }

    fn numeric_bound(
        &mut self,
        site: Site,
        function: Object,
        receiver: Object,
    ) -> Result<Object, H::Error> {
        let errors = self.host.literal(site, Cell::Text("coerce".into()))?;
        self.bound(
            site,
            "numeric",
            function,
            None,
            vec![receiver],
            vec![("errors", errors)],
        )
    }

    fn logical_not(&mut self, site: Site, receiver: Object) -> Result<Object, H::Error> {
        let function = self.host.global(site, Some("np"), "logical_not")?;
        self.bound(site, "not", function, None, vec![receiver], vec![])
    }

    fn column_contains(
        &mut self,
        site: Site,
        frame: Object,
        key: &'static str,
        pattern: &'static str,
    ) -> Result<Object, H::Error> {
        let value = self.host.column(site, frame, key)?;
        self.contains(site, value, pattern, false)
    }

    fn variant_mask(
        &mut self,
        site: Site,
        frame: Object,
        pattern: &'static str,
    ) -> Result<Object, H::Error> {
        let first = self.column_contains(site, frame, "Consequence", pattern)?;
        // np.logical_not's binding precedes evaluation of its getter argument.
        let function = self.host.global(site, Some("np"), "logical_not")?;
        let stop = self.column_contains(site, frame, "Consequence", "stop")?;
        let stop = self.bound(site, "not", function, None, vec![stop], vec![])?;
        let value = self.binary(site, "and", first, stop)?;
        let function = self.host.global(site, Some("np"), "logical_not")?;
        let frameshift = self.column_contains(site, frame, "Consequence", "frameshift")?;
        let frameshift = self.bound(site, "not", function, None, vec![frameshift], vec![])?;
        self.binary(site, "and", value, frameshift)
    }

    fn compare_literal(
        &mut self,
        site: Site,
        op: &'static str,
        left: Object,
        value: Cell,
    ) -> Result<Object, H::Error> {
        let right = self.host.literal(site, value)?;
        self.binary(site, op, left, right)
    }

    fn clinvar(
        &mut self,
        site: Site,
        frame: Object,
        confidence: Object,
        pattern: &'static str,
    ) -> Result<Object, H::Error> {
        let significance = self.host.column(site, frame, "CLNSIG")?;
        let significance = self.fillna(site, significance, Cell::Text(String::new()))?;
        let significance = self.contains(site, significance, pattern, false)?;
        let review = self.host.column(site, frame, "CLNREVSTAT")?;
        let function = self.host.attribute(site, review, "map")?;
        let ignore = self.host.literal(site, Cell::Text("ignore".into()))?;
        let review = self.bound(
            site,
            "map",
            function,
            Some(review),
            vec![review, confidence],
            vec![("na_action", ignore)],
        )?;
        let review = self.compare_literal(site, "ge", review, Cell::Int(2))?;
        self.binary(site, "and", significance, review)
    }

    fn emit(&mut self, site: Site, function: Object, message: Object) -> Result<(), H::Error> {
        self.bound(site, "log", function, None, vec![message], vec![])?;
        // The actual handler may change bindings or operands. No later value
        // was prefetched; the host rechecks the mature binding profile here.
        self.host.effect(site, "after_log", &[], &[])?;
        Ok(())
    }

    fn logger(&mut self, site: Site, level: &'static str) -> Result<Object, H::Error> {
        let logger = self.host.global(site, None, "logger")?;
        self.host.attribute(site, logger, level)
    }

    fn zeros(&mut self, site: Site, frame: Object) -> Result<Object, H::Error> {
        // Original CALL order: np.zeros binding, len(df), dtype=int binding.
        let function = self.host.global(site, Some("np"), "zeros")?;
        let count = self.host.length(site, frame)?;
        let dtype = self.host.global(site, None, "int")?;
        let result = self.execute(Call {
            site,
            operation: "zeros",
            target: Target::Bound {
                function,
                receiver: None,
            },
            args: vec![count],
            keywords: vec![("dtype", dtype)],
        })?;
        self.host.release(site, &[function, count, dtype])?;
        Ok(result)
    }

    fn assign(
        &mut self,
        site: Site,
        array: Object,
        mask: Object,
        value: Object,
    ) -> Result<Object, H::Error> {
        self.execute(Call {
            site,
            operation: "assign",
            target: Target::Assignment,
            args: vec![array, mask, value],
            keywords: vec![],
        })
    }
}

fn bp7(line: u32) -> Site {
    Site {
        function: "BP7_criteria",
        line,
    }
}
fn pp3(line: u32) -> Site {
    Site {
        function: "PP3_BP4_criteria",
        line,
    }
}
fn score(line: u32) -> Site {
    Site {
        function: "_coerce_score_column",
        line,
    }
}

/// The entire nested original helper, including coercion diagnostics and its
/// lazy warning/error frontier. Rust owns the branch and message composition;
/// the actual Python formatting/logger callbacks remain mature effects.
fn coerce_score<H: Host>(
    e: &mut Engine<'_, H>,
    frame: Object,
    column: Object,
    missing_default: Object,
) -> Result<Object, H::Error> {
    let raw = e.host.effect(score(45), "getitem", &[frame, column], &[])?;
    let numeric = e.numeric(score(46), raw)?;
    let nonmissing = e.method0(score(47), raw, "notna")?;
    let text = e.as_str(score(47), raw)?;
    let text = e.strip(score(47), text)?;
    let nonempty = e.compare_literal(score(47), "ne", text, Cell::Text(String::new()))?;
    let nonempty = e.binary(score(47), "and", nonmissing, nonempty)?;
    let invalid = e.method0(score(48), numeric, "isna")?;
    let invalid = e.binary(score(48), "and", nonempty, invalid)?;
    let any = e.method0(score(49), invalid, "any")?;
    if e.host.truth(score(49), any)? {
        // LOAD_GLOBAL sorted precedes the actual boolean Series getter.
        let sorted = e.host.global(score(50), None, "sorted")?;
        let selected = e.host.effect(score(50), "getitem", &[raw, invalid], &[])?;
        let selected = e.as_str(score(50), selected)?;
        let selected = e.strip(score(50), selected)?;
        let unique = e.method0(score(50), selected, "unique")?;
        let ordered = e.bound(score(50), "sorted", sorted, None, vec![unique], vec![])?;
        let examples = e.host.effect(score(50), "slice5", &[ordered], &[])?;
        let fill_message = if e.host.is_none(score(53), missing_default)? {
            e.host
                .literal(score(54), Cell::Text("leaving as missing".into()))?
        } else {
            e.host.effect(
                score(52),
                "format",
                &[missing_default],
                &["filling default=", ""],
            )?
        };
        let warning = e.logger(score(56), "warning")?;
        // Each FORMAT_VALUE is visited before the following source operand.
        let column_text = e.host.effect(score(57), "format", &[column], &["", ""])?;
        let count = e.method0(score(57), invalid, "sum")?;
        let count_text = e.host.effect(score(57), "format", &[count], &["", ""])?;
        let examples_text = e.host.effect(score(58), "format", &[examples], &["", ""])?;
        let fill_text = e
            .host
            .effect(score(58), "format", &[fill_message], &["", ""])?;
        let message = e.host.effect(
            score(57),
            "join_text",
            &[column_text, count_text, examples_text, fill_text],
            &[
                "PP3/BP4 score coercion: column '",
                "' has ",
                " non-numeric values (examples: ",
                "); treating as missing and ",
                "",
            ],
        )?;
        e.emit(score(56), warning, message)?;
    }
    if e.host.is_none(score(60), missing_default)? {
        return Ok(numeric);
    }
    let fillna = e.host.attribute(score(62), numeric, "fillna")?;
    e.bound(
        score(62),
        "fillna",
        fillna,
        Some(numeric),
        vec![numeric, missing_default],
        vec![],
    )
}

/// Whole original PP3/BP4 function. Every scientific branch, stage, literal
/// threshold, mask composition and output publication lives in this procedure.
/// There is no Python target-function call and no old completed-prefix replay.
pub fn pp3_bp4_criteria<H: Host>(
    host: &mut H,
    frame: Object,
    pvs1: Object,
    confidence: Object,
) -> Result<(Object, Object), H::Error> {
    let mut e = Engine::new(host);
    let missense = e.variant_mask(pp3(32), frame, "missense_variant")?;
    let splice = e.variant_mask(pp3(33), frame, "splice")?;
    let utr = e.column_contains(pp3(34), frame, "Consequence", "5_prime_UTR_variant")?;
    let utr = e.fillna(pp3(34), utr, Cell::Bool(false))?;
    let missense = e.fillna(pp3(35), missense, Cell::Bool(false))?;
    let splice = e.fillna(pp3(36), splice, Cell::Bool(false))?;

    // df.get is bound before its eagerly evaluated default Series. The
    // original .index access/cache and Series construction are mature here.
    let get = e.host.attribute(pp3(37), frame, "get")?;
    let lof_key = e.host.literal(pp3(38), Cell::Text("LoF".into()))?;
    let series = e.host.global(pp3(38), Some("pd"), "Series")?;
    let empty = e.host.literal(pp3(38), Cell::Text(String::new()))?;
    let index = e.host.attribute(pp3(38), frame, "index")?;
    let default = e.bound(
        pp3(38),
        "series_default",
        series,
        None,
        vec![empty],
        vec![("index", index)],
    )?;
    let lof = e.bound(
        pp3(37),
        "get",
        get,
        Some(frame),
        vec![frame, lof_key, default],
        vec![],
    )?;
    let lof = e.fillna(pp3(39), lof, Cell::Text(String::new()))?;
    let lof = e.as_str(pp3(39), lof)?;
    let lof = e.contains(pp3(39), lof, r"(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)", true)?;
    let loftee_splice = e.binary(pp3(42), "and", splice, lof)?;
    // Evaluate original annotations/defaults at this definition frontier;
    // the resulting Python helper is NEVER called. Its body is above in Rust.
    e.host
        .effect(pp3(44), "define_score_helper", &[frame], &[])?;

    let column = e.host.literal(pp3(64), Cell::Text("PrimateAI".into()))?;
    let zero = e.host.literal(pp3(64), Cell::Int(0))?;
    let primate = coerce_score(&mut e, frame, column, zero)?;
    let column = e
        .host
        .literal(pp3(65), Cell::Text("am_pathogenicity".into()))?;
    let zero = e.host.literal(pp3(65), Cell::Int(0))?;
    let am = coerce_score(&mut e, frame, column, zero)?;
    let column = e.host.literal(pp3(66), Cell::Text("CADD_phred".into()))?;
    let none = e.host.literal(pp3(66), Cell::None)?;
    let cadd = coerce_score(&mut e, frame, column, none)?;
    let column = e
        .host
        .literal(pp3(67), Cell::Text("CADD_reg_phred".into()))?;
    let none = e.host.literal(pp3(67), Cell::None)?;
    let cadd_reg = coerce_score(&mut e, frame, column, none)?;
    let cadd_high = e.compare_literal(pp3(68), "ge", cadd, Cell::Int(20))?;
    let cutoff = e.host.literal(pp3(71), Cell::Float(15.0))?;
    let cadd_present = e.method0(pp3(72), cadd, "notna")?;
    let cadd_low = e.binary(pp3(72), "lt", cadd, cutoff)?;
    let cadd_low = e.binary(pp3(72), "and", cadd_present, cadd_low)?;
    let cadd_reg_high = e.binary(pp3(73), "ge", cadd_reg, cutoff)?;
    let cadd_reg_high = e.fillna(pp3(73), cadd_reg_high, Cell::Bool(false))?;
    let cadd_reg_not_high = e.invert(pp3(73), cadd_reg_high)?;
    let splice_lof = e.host.column(pp3(74), frame, "splicing_lof")?;
    let splice_lof = e.fillna(pp3(74), splice_lof, Cell::Bool(false))?;
    let splice_lof = e.binary(pp3(74), "or", splice_lof, loftee_splice)?;

    // Bind concat BEFORE the list comprehension; to_numeric is bound before
    // each actual df[c] getter, exactly as in the original comprehension.
    let concat = e.host.global(pp3(80), Some("pd"), "concat")?;
    let mut delta_columns = Vec::with_capacity(4);
    for key in [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
    ] {
        let numeric = e.host.global(pp3(81), Some("pd"), "to_numeric")?;
        let raw = e.host.column(pp3(81), frame, key)?;
        delta_columns.push(e.numeric_bound(pp3(81), numeric, raw)?);
    }
    let columns = e.host.effect(pp3(81), "list", &delta_columns, &[])?;
    let axis = e.host.literal(pp3(81), Cell::Int(1))?;
    let matrix = e.bound(
        pp3(80),
        "concat",
        concat,
        None,
        vec![columns],
        vec![("axis", axis)],
    )?;
    let max = e.host.attribute(pp3(82), matrix, "max")?;
    let axis = e.host.literal(pp3(82), Cell::Int(1))?;
    let ds = e.bound(
        pp3(82),
        "max_columns",
        max,
        Some(matrix),
        vec![matrix],
        vec![("axis", axis)],
    )?;
    let ds_present = e.method0(pp3(83), ds, "notna")?;
    let ds_low = e.compare_literal(pp3(83), "lt", ds, Cell::Float(0.1))?;
    let ds_no_impact = e.binary(pp3(83), "and", ds_present, ds_low)?;

    let primate_high = e.compare_literal(pp3(86), "gt", primate, Cell::Float(0.8))?;
    let mut pp3_mask = e.binary(pp3(86), "and", primate_high, missense)?;
    let not_splice = e.logical_not(pp3(87), splice)?;
    let mut other_high = e.binary(pp3(87), "and", cadd_high, not_splice)?;
    let not_missense = e.logical_not(pp3(87), missense)?;
    other_high = e.binary(pp3(87), "and", other_high, not_missense)?;
    let not_utr = e.logical_not(pp3(87), utr)?;
    other_high = e.binary(pp3(87), "and", other_high, not_utr)?;
    pp3_mask = e.binary(pp3(86), "or", pp3_mask, other_high)?;
    let am_class = e.host.column(pp3(88), frame, "am_class")?;
    let am_class = e.fillna(pp3(88), am_class, Cell::Text(String::new()))?;
    let am_class = e.contains(pp3(88), am_class, "pathogenic", false)?;
    let am_class = e.binary(pp3(88), "and", am_class, missense)?;
    pp3_mask = e.binary(pp3(86), "or", pp3_mask, am_class)?;
    let mut vep = e.host.column(pp3(89), frame, "vep_consq_lof")?;
    let not_splice = e.logical_not(pp3(89), splice)?;
    vep = e.binary(pp3(89), "and", vep, not_splice)?;
    let not_missense = e.logical_not(pp3(89), missense)?;
    vep = e.binary(pp3(89), "and", vep, not_missense)?;
    let not_utr = e.logical_not(pp3(89), utr)?;
    vep = e.binary(pp3(89), "and", vep, not_utr)?;
    pp3_mask = e.binary(pp3(86), "or", pp3_mask, vep)?;
    let splice_or_cadd = e.binary(pp3(90), "or", splice_lof, cadd_high)?;
    let splice_evidence = e.binary(pp3(90), "and", splice_or_cadd, splice)?;
    let utr_lof = e.host.column(pp3(90), frame, "5UTR_lof")?;
    let utr_evidence = e.binary(pp3(90), "and", utr_lof, utr)?;
    let class_evidence = e.binary(pp3(90), "or", splice_evidence, utr_evidence)?;
    pp3_mask = e.binary(pp3(86), "or", pp3_mask, class_evidence)?;
    let benign = e.clinvar(pp3(91), frame, confidence, "enign")?;
    let not_benign = e.invert(pp3(92), benign)?;
    pp3_mask = e.binary(pp3(92), "and", pp3_mask, not_benign)?;
    if !e.host.is_none(pp3(98), pvs1)? {
        let weak = e.compare_literal(pp3(103), "lt", pvs1, Cell::Int(3))?;
        pp3_mask = e.binary(pp3(103), "and", pp3_mask, weak)?;
        let logger = e.logger(pp3(104), "info")?;
        let strong = e.compare_literal(pp3(104), "ge", pvs1, Cell::Int(3))?;
        let count = e.method0(pp3(104), strong, "sum")?;
        let message = e.host.effect(
            pp3(104),
            "format",
            &[count],
            &[
                "PP3 double-counting prevention: blocked ",
                " variants with PVS1_Strong or stronger",
            ],
        )?;
        e.emit(pp3(104), logger, message)?;
    }

    let primate_low = e.compare_literal(pp3(106), "lt", primate, Cell::Float(0.8))?;
    let primate_low = e.fillna(pp3(106), primate_low, Cell::Bool(true))?;
    let am_low = e.compare_literal(pp3(106), "lt", am, Cell::Float(0.564))?;
    let am_low = e.fillna(pp3(106), am_low, Cell::Bool(true))?;
    let missense_benign = e.binary(pp3(106), "and", primate_low, am_low)?;
    let missense_benign = e.binary(pp3(106), "and", missense_benign, missense)?;
    let splice_benign = e.logical_not(pp3(107), splice_lof)?;
    let splice_benign = e.binary(pp3(107), "and", splice_benign, splice)?;
    let cadd_or_ds = e.binary(pp3(109), "or", cadd_low, ds_no_impact)?;
    let splice_benign = e.binary(pp3(107), "and", splice_benign, cadd_or_ds)?;
    let logical_not = e.host.global(pp3(110), Some("np"), "logical_not")?;
    let utr_lof = e.host.column(pp3(110), frame, "5UTR_lof")?;
    let utr_lof = e.fillna(pp3(110), utr_lof, Cell::Bool(false))?;
    let utr_benign = e.bound(pp3(110), "not", logical_not, None, vec![utr_lof], vec![])?;
    let utr_benign = e.binary(pp3(110), "and", utr_benign, utr)?;
    let logical_not = e.host.global(pp3(111), Some("np"), "logical_not")?;
    let vep_lof = e.host.column(pp3(111), frame, "vep_consq_lof")?;
    let vep_lof = e.fillna(pp3(111), vep_lof, Cell::Bool(false))?;
    let mut other_benign = e.bound(pp3(111), "not", logical_not, None, vec![vep_lof], vec![])?;
    other_benign = e.binary(pp3(111), "and", other_benign, cadd_low)?;
    other_benign = e.binary(pp3(111), "and", other_benign, cadd_reg_not_high)?;
    let not_splice = e.logical_not(pp3(111), splice)?;
    other_benign = e.binary(pp3(111), "and", other_benign, not_splice)?;
    let not_missense = e.logical_not(pp3(111), missense)?;
    other_benign = e.binary(pp3(111), "and", other_benign, not_missense)?;
    let not_utr = e.logical_not(pp3(111), utr)?;
    other_benign = e.binary(pp3(111), "and", other_benign, not_utr)?;
    let mut bp4_mask = e.binary(pp3(112), "or", missense_benign, splice_benign)?;
    bp4_mask = e.binary(pp3(112), "or", bp4_mask, utr_benign)?;
    bp4_mask = e.binary(pp3(112), "or", bp4_mask, other_benign)?;
    let pathogenic = e.clinvar(pp3(113), frame, confidence, "athogenic")?;
    let not_pathogenic = e.invert(pp3(114), pathogenic)?;
    bp4_mask = e.binary(pp3(114), "and", bp4_mask, not_pathogenic)?;
    let not_pp3 = e.invert(pp3(115), pp3_mask)?;
    bp4_mask = e.binary(pp3(115), "and", bp4_mask, not_pp3)?;
    if !e.host.is_none(pp3(116), pvs1)? {
        let asarray = e.host.global(pp3(119), Some("np"), "asarray")?;
        let pvs_array = e.bound(pp3(119), "asarray", asarray, None, vec![pvs1], vec![])?;
        let strong = e.compare_literal(pp3(119), "ge", pvs_array, Cell::Int(3))?;
        let int = e.host.global(pp3(120), None, "int")?;
        let suppressed = e.binary(pp3(120), "and", bp4_mask, strong)?;
        let suppressed = e.method0(pp3(120), suppressed, "sum")?;
        let suppressed = e.bound(pp3(120), "int", int, None, vec![suppressed], vec![])?;
        let not_strong = e.invert(pp3(121), strong)?;
        bp4_mask = e.binary(pp3(121), "and", bp4_mask, not_strong)?;
        let logger = e.logger(pp3(122), "info")?;
        let message = e.host.effect(
            pp3(122),
            "format",
            &[suppressed],
            &[
                "BP4 suppression: blocked ",
                " variants with PVS1_Strong or stronger",
            ],
        )?;
        e.emit(pp3(122), logger, message)?;
    }
    let pp3_array = e.zeros(pp3(124), frame)?;
    let bp4_array = e.zeros(pp3(125), frame)?;
    let one = e.host.literal(pp3(126), Cell::Int(1))?;
    let pp3_array = e.assign(pp3(126), pp3_array, pp3_mask, one)?;
    let one = e.host.literal(pp3(127), Cell::Int(1))?;
    let bp4_array = e.assign(pp3(127), bp4_array, bp4_mask, one)?;
    Ok((pp3_array, bp4_array))
}

/// Entire original BP7 procedure. Its literal OR and repeated mature getters
/// are intentional, even when an earlier getter returned the same column.
/// The original frame/axes/custom getters are NOT admitted or prefetched here.
/// Native result ownership/failure is serviced at each visited source site.
pub fn bp7_criteria<H: Host>(host: &mut H, frame: Object) -> Result<Object, H::Error> {
    let mut e = Engine::new(host);
    let consequence = e.host.column(bp7(133), frame, "Consequence")?;
    let synonymous_literal = e
        .host
        .literal(bp7(133), Cell::Text("synonymous_variant".into()))?;
    let synonymous = e.binary(bp7(133), "eq", consequence, synonymous_literal)?;
    e.host
        .release(bp7(133), &[consequence, synonymous_literal])?;

    let length = e.host.column(bp7(134), frame, "splicing_len_changing")?;
    let false_literal = e.host.literal(bp7(134), Cell::Bool(false))?;
    let length_false = e.binary(bp7(134), "eq", length, false_literal)?;
    e.host.release(bp7(134), &[length, false_literal])?;
    let lof = e.host.column(bp7(134), frame, "splicing_lof")?;
    let lof_missing = e.method0(bp7(134), lof, "isna")?;
    e.host.release(bp7(134), &[lof])?;
    let first_disjunct = e.binary(bp7(134), "or", length_false, lof_missing)?;
    e.host.release(bp7(134), &[length_false, lof_missing])?;
    let lof_again = e.host.column(bp7(134), frame, "splicing_lof")?;
    let false_again = e.host.literal(bp7(134), Cell::Bool(false))?;
    let lof_false = e.binary(bp7(134), "eq", lof_again, false_again)?;
    e.host.release(bp7(134), &[lof_again, false_again])?;
    let second_disjunct = e.binary(bp7(134), "or", first_disjunct, lof_false)?;
    e.host.release(bp7(134), &[first_disjunct, lof_false])?;
    let length_again = e.host.column(bp7(134), frame, "splicing_len_changing")?;
    let length_missing = e.method0(bp7(134), length_again, "isna")?;
    e.host.release(bp7(134), &[length_again])?;
    let no_splicing = e.binary(bp7(134), "or", second_disjunct, length_missing)?;
    e.host
        .release(bp7(134), &[second_disjunct, length_missing])?;

    let conservation = e.host.column(bp7(135), frame, "Conservation")?;
    let cutoff = e.host.literal(bp7(135), Cell::Int(5))?;
    let not_conserved = e.binary(bp7(135), "le", conservation, cutoff)?;
    e.host.release(bp7(135), &[conservation, cutoff])?;
    let array = e.zeros(bp7(136), frame)?;
    // STORE_SUBSCR evaluates RHS before the target/mask expression.
    let one = e.host.literal(bp7(137), Cell::Int(1))?;
    let first_mask = e.binary(bp7(137), "and", synonymous, no_splicing)?;
    let mask = e.binary(bp7(137), "and", first_mask, not_conserved)?;
    e.host.release(bp7(137), &[first_mask])?;
    let output = e.assign(bp7(137), array, mask, one)?;
    e.host.release(bp7(137), &[one, array, mask])?;
    e.host
        .release(bp7(138), &[synonymous, no_splicing, not_conserved])?;
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    struct FakeHost {
        values: BTreeMap<Object, Value>,
        columns: BTreeMap<&'static str, Object>,
        next: Object,
        events: Vec<String>,
        getter_count: usize,
        fail_getter: Option<usize>,
    }

    impl FakeHost {
        fn put(&mut self, value: Value) -> Object {
            let id = self.next;
            self.next += 1;
            self.values.insert(id, value);
            id
        }
        fn example() -> Self {
            let mut h = Self {
                values: BTreeMap::new(),
                columns: BTreeMap::new(),
                next: 10,
                events: vec![],
                getter_count: 0,
                fail_getter: None,
            };
            let columns = [
                (
                    "Consequence",
                    "object",
                    vec![
                        Cell::Text("synonymous_variant".into()),
                        Cell::Text("synonymous_variant".into()),
                        Cell::Text("synonymous_variant".into()),
                        Cell::Text("synonymous_variant".into()),
                        Cell::Text("missense_variant".into()),
                        Cell::Text("synonymous_variant".into()),
                    ],
                ),
                (
                    "splicing_len_changing",
                    "bool",
                    vec![
                        Cell::Bool(true),
                        Cell::Bool(true),
                        Cell::Bool(true),
                        Cell::Bool(false),
                        Cell::Bool(false),
                        Cell::Bool(false),
                    ],
                ),
                (
                    "splicing_lof",
                    "boolean",
                    vec![
                        Cell::Bool(false),
                        Cell::Bool(true),
                        Cell::Na,
                        Cell::Bool(true),
                        Cell::Bool(false),
                        Cell::Bool(false),
                    ],
                ),
                (
                    "Conservation",
                    "float64",
                    vec![
                        Cell::Float(5.0),
                        Cell::Float(5.0),
                        Cell::Float(5.0),
                        Cell::Float(5.0),
                        Cell::Float(4.0),
                        Cell::Float(5.0001),
                    ],
                ),
            ];
            for (name, dtype, cells) in columns {
                let index = (0..cells.len()).map(|i| format!("i{i:020}")).collect();
                let id = h.put(Value {
                    kind: "series".into(),
                    dtype: dtype.into(),
                    index,
                    cells,
                });
                h.columns.insert(name, id);
            }
            h
        }
    }

    impl Host for FakeHost {
        type Error = String;
        fn literal(&mut self, _site: Site, value: Cell) -> Result<Object, String> {
            Ok(self.put(Value::scalar(value)))
        }
        fn column(
            &mut self,
            _site: Site,
            _frame: Object,
            key: &'static str,
        ) -> Result<Object, String> {
            self.events.push(format!("get:{key}"));
            self.getter_count += 1;
            if self.fail_getter == Some(self.getter_count) {
                return Err("actual getter failure".into());
            }
            let prototype = self.columns.get(key).copied().ok_or("missing column")?;
            Ok(self.put(self.values[&prototype].clone()))
        }
        fn attribute(
            &mut self,
            _site: Site,
            receiver: Object,
            _name: &'static str,
        ) -> Result<Object, String> {
            Ok(self.put(self.values[&receiver].clone()))
        }
        fn global(
            &mut self,
            _site: Site,
            _module: Option<&'static str>,
            _name: &'static str,
        ) -> Result<Object, String> {
            Ok(self.put(Value::scalar(Cell::None)))
        }
        fn length(&mut self, _site: Site, _object: Object) -> Result<Object, String> {
            self.events.push("len".into());
            Ok(self.put(Value::scalar(Cell::Int(6))))
        }
        fn effect(
            &mut self,
            _site: Site,
            _action: &'static str,
            _args: &[Object],
            _text: &[&str],
        ) -> Result<Object, String> {
            Err("unexpected mature effect".into())
        }
        fn truth(&mut self, _site: Site, _object: Object) -> Result<bool, String> {
            Err("unexpected mature truth test".into())
        }
        fn is_none(&mut self, _site: Site, _object: Object) -> Result<bool, String> {
            Err("unexpected mature identity test".into())
        }
        fn release(&mut self, _site: Site, objects: &[Object]) -> Result<(), String> {
            for object in objects {
                self.values
                    .remove(object)
                    .ok_or("released foreign/duplicate object")?;
            }
            Ok(())
        }
        fn prepare(&mut self, call: &Call) -> Result<Option<Vec<Value>>, String> {
            Ok(Some(
                call.args.iter().map(|id| self.values[id].clone()).collect(),
            ))
        }
        fn commit_native(&mut self, call: &Call, value: Value) -> Result<Object, String> {
            self.events.push(format!("native:{}", call.operation));
            Ok(self.put(value))
        }
        fn mature(&mut self, call: &Call) -> Result<Object, String> {
            Err(format!("mature-required:{}", call.operation))
        }
    }

    #[test]
    fn bp7_named_procedure_keeps_literal_or_and_exact_conservation_edge() {
        let mut host = FakeHost::example();
        let result = bp7_criteria(&mut host, 0).unwrap();
        assert_eq!(
            host.values[&result].cells,
            [
                Cell::Int(1),
                Cell::Int(0),
                Cell::Int(1),
                Cell::Int(1),
                Cell::Int(0),
                Cell::Int(0)
            ]
        );
    }

    #[test]
    fn bp7_repeats_actual_getters_in_source_order_without_prefetch() {
        let mut host = FakeHost::example();
        bp7_criteria(&mut host, 0).unwrap();
        let actual: Vec<_> = host
            .events
            .iter()
            .filter(|event| event.starts_with("get:"))
            .map(String::as_str)
            .collect();
        assert_eq!(
            actual,
            [
                "get:Consequence",
                "get:splicing_len_changing",
                "get:splicing_lof",
                "get:splicing_lof",
                "get:splicing_len_changing",
                "get:Conservation"
            ]
        );
    }

    #[test]
    fn bp7_actual_getter_failure_stops_before_later_reads_or_output_allocation() {
        let mut host = FakeHost::example();
        host.fail_getter = Some(4);
        assert_eq!(
            bp7_criteria(&mut host, 0),
            Err("actual getter failure".into())
        );
        assert_eq!(host.getter_count, 4);
        assert!(
            !host
                .events
                .iter()
                .any(|event| event == "get:Conservation" || event == "native:zeros")
        );
    }

    #[test]
    fn bp7_allocates_before_final_mask_and_publishes_only_after_assignment() {
        let mut host = FakeHost::example();
        let result = bp7_criteria(&mut host, 0).unwrap();
        let zeros = host
            .events
            .iter()
            .position(|event| event == "native:zeros")
            .unwrap();
        let first_and = host
            .events
            .iter()
            .position(|event| event == "native:and")
            .unwrap();
        assert!(zeros < first_and);
        assert_eq!(host.events.last().unwrap(), "native:assign");
        assert_eq!(host.values[&result].kind, "array");
    }

    #[test]
    fn bp7_releases_consumed_temporaries_without_reusing_a_released_handle() {
        let mut host = FakeHost::example();
        let result = bp7_criteria(&mut host, 0).unwrap();
        assert_eq!(host.values.len(), host.columns.len() + 1);
        assert!(host.values.contains_key(&result));
        for prototype in host.columns.values() {
            assert!(host.values.contains_key(prototype));
        }
    }
}
