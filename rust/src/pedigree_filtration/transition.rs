//! Explicit stateless v2 transport, delegating the original v1 decision bodies.
//!
//! 2026-10-27: root-approved pure-transition reuse. The current record's parser,
//! getter, logging and write events stay on the host. No hidden Filter survives
//! a transition. Raw GT tags retain list/tuple and bool/int distinctions even
//! though the original mathematical Genotype representation coalesces some.
use super::*;

pub(super) const VERSION: &str = "priva-pedigree-transition-v2";

fn prefixed(fields: &Fields, prefix: &str) -> Fields {
    fields
        .iter()
        .filter_map(|(key, value)| {
            key.strip_prefix(prefix)
                .map(|key| (key.to_owned(), value.clone()))
        })
        .collect()
}

fn add_prefixed(output: &mut Fields, prefix: &str, fields: &Fields) {
    for (key, value) in fields {
        output.insert(format!("{prefix}{key}"), value.clone());
    }
}

fn phase_name(phase: Phase) -> (&'static str, usize) {
    match phase {
        Phase::Idle => ("Idle", 0),
        Phase::Patient(index) => ("Patient", index),
        Phase::Control(index) => ("Control", index),
        Phase::Chrom => ("Chrom", 0),
        Phase::Proband => ("Proband", 0),
        Phase::Father => ("Father", 0),
        Phase::Mother => ("Mother", 0),
        Phase::DebugAuto => ("DebugAuto", 0),
        Phase::DebugDrop => ("DebugDrop", 0),
        Phase::Done => ("Done", 0),
    }
}

fn parse_phase(name: &str, index: usize) -> Result<Phase> {
    if !["Patient", "Control"].contains(&name) && index != 0 {
        return Err(invalid("nonindexed decision phase has an index"));
    }
    Ok(match name {
        "Idle" => Phase::Idle,
        "Patient" => Phase::Patient(index),
        "Control" => Phase::Control(index),
        "Chrom" => Phase::Chrom,
        "Proband" => Phase::Proband,
        "Father" => Phase::Father,
        "Mother" => Phase::Mother,
        "DebugAuto" => Phase::DebugAuto,
        "DebugDrop" => Phase::DebugDrop,
        "Done" => Phase::Done,
        _ => return Err(invalid("unknown decision phase")),
    })
}

/// Reconstruct only a private calculation value; never acquire a PED/VCF input.
fn context_filter(context: &Fields) -> Result<Filter> {
    let mut fields = context.clone();
    for role in ["father", "mother"] {
        let kind = take(&mut fields, &format!("{role}.original_kind"))?;
        if !["none", "str"].contains(&kind.as_str())
            || (kind == "none" && fields.get(role).map(String::as_str) != Some(""))
        {
            return Err(invalid("invalid original parent kind"));
        }
    }
    if take(&mut fields, "sex.original_kind")? != "str" {
        return Err(invalid("unadmitted parsed sex kind"));
    }
    let state = Filter::new(&mut fields)?;
    finished(fields)?;
    Ok(state)
}

fn raw_gt(fields: &Fields, prefix: &str) -> Fields {
    fields
        .iter()
        .filter(|(name, _)| name.starts_with(&format!("{prefix}.")))
        .map(|(name, value)| (name.clone(), value.clone()))
        .collect()
}

fn initial_raw() -> Fields {
    ["proband", "father", "mother"]
        .into_iter()
        .map(|role| (format!("{role}.kind"), "none".into()))
        .collect()
}

fn snapshot(state: &Filter, raw: &Fields) -> Fields {
    let (name, index) = phase_name(state.phase);
    let mut result = raw.clone();
    result.insert("phase".into(), name.into());
    result.insert("index".into(), index.to_string());
    result.insert("chrom_type".into(), state.chrom_type.clone());
    result
}

fn restore(context: &Fields, state: &Fields) -> Result<Filter> {
    let mut result = context_filter(context)?;
    let mut fields = state.clone();
    let name = take(&mut fields, "phase")?;
    let index = number(&mut fields, "index")?;
    result.phase = parse_phase(&name, index)?;
    if (name == "Patient" && index >= result.patients.len())
        || (name == "Control" && index >= result.controls.len())
    {
        return Err(invalid("decision participant index exceeds context"));
    }
    result.chrom_type = take(&mut fields, "chrom_type")?;
    if !["", "AUTOSOMAL", "X", "Y"].contains(&result.chrom_type.as_str()) {
        return Err(invalid("invalid decision chromosome classification"));
    }
    result.proband_gt = genotype(&mut fields, "proband")?;
    result.father_gt = genotype(&mut fields, "father")?;
    result.mother_gt = genotype(&mut fields, "mother")?;
    finished(fields)?;
    // Actual record ordinals are checked by the host on EVERY visited begin.
    // The only old ordinal read is begin's predecessor+1 check. A private
    // 0 -> 1 delegate exercises the identical guard without making the pure
    // semantic state depend on a record identifier or retaining a record.
    result.ordinal = 0;
    Ok(result)
}

fn initialize(context: &Fields) -> Result<Fields> {
    let state = context_filter(context)?;
    Ok(snapshot(&state, &initial_raw()))
}

/// Pure: the inputs remain immutable, including on a source-error outcome.
fn execute(context: &Fields, state: &Fields, operation: &str, message: &Fields) -> Result<Fields> {
    if ![
        "genotype",
        "chromosome",
        "route",
        "source",
        "filter-begin",
        "filter-reply",
    ]
    .contains(&operation)
    {
        return Err(invalid("unqualified pure transition operation"));
    }
    if message.contains_key("op") || message.contains_key("sequence") {
        return Err(invalid("reserved pure transition field"));
    }
    let mut filter = if context.is_empty() {
        if !state.is_empty() {
            return Err(invalid("decision state without initialized context"));
        }
        None
    } else {
        Some(restore(context, state)?)
    };
    let mut raw: Fields = state
        .iter()
        .filter(|(name, _)| {
            ["proband.", "father.", "mother."]
                .iter()
                .any(|prefix| name.starts_with(prefix))
        })
        .map(|(name, value)| (name.clone(), value.clone()))
        .collect();
    // Preserve the exact typed message when an existing rule stores a GT,
    // including its value/tag before a subsequent has_variant error.
    let stored_role = if operation == "filter-reply" {
        filter.as_ref().and_then(|state| match state.phase {
            Phase::Proband => Some("proband"),
            Phase::Father => Some("father"),
            Phase::Mother => Some("mother"),
            _ => None,
        })
    } else {
        None
    };
    let mut fields = message.clone();
    if operation == "filter-begin" {
        if fields.contains_key("ordinal") {
            return Err(invalid(
                "actual ordinal belongs to the checked outer envelope",
            ));
        }
        fields.insert("ordinal".into(), "1".into());
    }
    fields.insert("op".into(), operation.into());
    fields.insert("sequence".into(), "1".into());
    // This is the preserved source-derived decision implementation, not a
    // second Python implementation, guessed action, or cached read observation.
    let (mut outcome, stop) = dispatch_request(fields, &mut filter)?;
    if stop || outcome.remove("sequence").as_deref() != Some("1") {
        return Err(invalid("invalid pure transition envelope"));
    }
    if let Some(role) = stored_role {
        raw.retain(|name, _| !name.starts_with(&format!("{role}.")));
        for (name, value) in raw_gt(message, "gt") {
            raw.insert(format!("{role}.{}", &name[3..]), value);
        }
    }
    let mut result = Fields::new();
    if let Some(state) = filter {
        add_prefixed(&mut result, "state.", &snapshot(&state, &raw));
    }
    add_prefixed(&mut result, "outcome.", &outcome);
    Ok(result)
}

fn transport(mut fields: Fields) -> Result<(Fields, bool)> {
    let op = fields
        .get("op")
        .cloned()
        .ok_or_else(|| invalid("missing operation"))?;
    if ["family", "parents"].contains(&op.as_str()) {
        // One-time Polars family/parent selection remains exactly the old code.
        return dispatch_request(fields, &mut None);
    }
    take(&mut fields, "op")?;
    let sequence = number(&mut fields, "sequence")?;
    let mut stop = false;
    let mut output = match op.as_str() {
        "hello" => {
            if take(&mut fields, "protocol")? != VERSION {
                return Err(invalid("unqualified stateless pedigree protocol"));
            }
            outcome(Ok(VERSION))
        }
        "quit" => {
            stop = true;
            outcome(Ok("closed"))
        }
        "kernel-init" => {
            let context = prefixed(&fields, "context.");
            fields.retain(|name, _| !name.starts_with("context."));
            let state = initialize(&context)?;
            let mut result = outcome(Ok("ready"));
            add_prefixed(&mut result, "context.", &context);
            add_prefixed(&mut result, "state.", &state);
            result
        }
        "transition" => {
            let operation = take(&mut fields, "procedure")?;
            let context = prefixed(&fields, "context.");
            let state = prefixed(&fields, "state.");
            let message = prefixed(&fields, "message.");
            fields.retain(|name, _| {
                !["context.", "state.", "message."]
                    .iter()
                    .any(|prefix| name.starts_with(prefix))
            });
            let mut result = execute(&context, &state, &operation, &message)?;
            result.insert("status".into(), "ok".into());
            result
        }
        _ => return Err(invalid("unknown stateless pedigree transport operation")),
    };
    finished(fields)?;
    output.insert("sequence".into(), sequence.to_string());
    Ok((output, stop))
}

fn serve_stateless(input: &mut impl Read, output: &mut impl Write) -> Result<()> {
    let mut sequence = 1;
    while let Some(fields) = read_frame(input)? {
        if fields.get("sequence") != Some(&sequence.to_string()) {
            return Err(invalid("stateless pedigree wire sequence changed"));
        }
        let hello = fields.get("op").map(String::as_str) == Some("hello");
        if hello != (sequence == 1) {
            return Err(invalid("stateless pedigree handshake order changed"));
        }
        let (response, stop) = transport(fields)?;
        write_frame(output, response)?;
        if stop {
            return Ok(());
        }
        sequence += 1;
    }
    Err(invalid(
        "stateless pedigree session ended before explicit close",
    ))
}

pub(super) fn dispatch() -> Result<()> {
    serve_stateless(&mut io::stdin().lock(), &mut io::stdout().lock())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(values: &[(&str, &str)]) -> Fields {
        values
            .iter()
            .map(|(key, value)| ((*key).into(), (*value).into()))
            .collect()
    }

    fn context() -> Fields {
        fields(&[
            ("patients.count", "1"),
            ("patients.0", "p"),
            ("controls.count", "0"),
            ("proband", "p"),
            ("father", "f"),
            ("mother", "m"),
            ("sex", "1"),
            ("father.original_kind", "str"),
            ("mother.original_kind", "str"),
            ("sex.original_kind", "str"),
        ])
    }

    fn state_after(response: &Fields) -> Fields {
        prefixed(response, "state.")
    }

    #[test]
    fn inputs_are_immutable_and_complete_outcomes_are_repeatable() {
        let context = context();
        let state = initialize(&context).unwrap();
        let before = state.clone();
        let left = execute(&context, &state, "filter-begin", &Fields::new()).unwrap();
        let right = execute(&context, &state, "filter-begin", &Fields::new()).unwrap();
        assert_eq!(state, before);
        assert_eq!(left, right);
        assert_eq!(left["outcome.action"], "patient-gt");
        assert_eq!(state_after(&left)["phase"], "Patient");
        assert!(!state_after(&left).contains_key("ordinal"));
    }

    #[test]
    fn stored_source_error_keeps_exact_typed_successor() {
        let context = context();
        let mut state = initialize(&context).unwrap();
        state.insert("phase".into(), "Proband".into());
        let response = execute(
            &context,
            &state,
            "filter-reply",
            &fields(&[("gt.kind", "int"), ("gt.value", "19")]),
        )
        .unwrap();
        assert_eq!(response["outcome.status"], "source-error");
        assert_eq!(response["outcome.class"], "TypeError");
        assert_eq!(response["state.proband.kind"], "int");
        assert_eq!(response["state.proband.value"], "19");
        assert_eq!(response["state.phase"], "Proband");
        assert_eq!(state["proband.kind"], "none");
    }

    #[test]
    fn tuple_list_and_boolean_integer_tags_are_not_erased() {
        for kind in ["tuple", "list"] {
            for allele in ["bool", "int"] {
                let context = context();
                let mut state = initialize(&context).unwrap();
                state.insert("phase".into(), "Father".into());
                let response = execute(
                    &context,
                    &state,
                    "filter-reply",
                    &fields(&[
                        ("gt.kind", kind),
                        ("gt.count", "1"),
                        ("gt.0.kind", allele),
                        ("gt.0.value", "1"),
                    ]),
                )
                .unwrap();
                assert_eq!(response["state.father.kind"], kind);
                assert_eq!(response["state.father.0.kind"], allele);
                assert_eq!(response["state.phase"], "Mother");
            }
        }
    }

    #[test]
    fn constructor_side_effect_and_duplicate_guard_are_preserved() {
        let context = context();
        let mut state = initialize(&context).unwrap();
        state.insert("phase".into(), "Chrom".into());
        let response = execute(
            &context,
            &state,
            "chromosome",
            &fields(&[("chrom", "chrX")]),
        )
        .unwrap();
        assert_eq!(response["state.chrom_type"], "X");
        let next = state_after(&response);
        assert!(execute(&context, &next, "chromosome", &fields(&[("chrom", "chrX")])).is_err());
        assert!(execute(&context, &state, "filter-reply", &fields(&[("ack", "1")])).is_err());
    }

    #[test]
    fn unused_proband_failure_is_after_debug_ack_and_parent_checks() {
        let context = context();
        let mut state = initialize(&context).unwrap();
        state.insert("phase".into(), "DebugAuto".into());
        state.insert("chrom_type".into(), "AUTOSOMAL".into());
        state.insert("proband.kind".into(), "int".into());
        state.insert("proband.value".into(), "7".into());
        for role in ["father", "mother"] {
            state.insert(format!("{role}.kind"), "tuple".into());
            state.insert(format!("{role}.count"), "1".into());
            state.insert(format!("{role}.0.kind"), "int".into());
            state.insert(format!("{role}.0.value"), "0".into());
        }
        let response = execute(&context, &state, "filter-reply", &fields(&[("ack", "1")])).unwrap();
        assert_eq!(response["outcome.status"], "source-error");
        assert_eq!(response["state.phase"], "DebugAuto");
        for role in ["father", "mother"] {
            state.retain(|name, _| !name.starts_with(&format!("{role}.")));
            state.insert(format!("{role}.kind"), "none".into());
        }
        let response = execute(&context, &state, "filter-reply", &fields(&[("ack", "1")])).unwrap();
        assert_eq!(response["outcome.value"], "UNKNOWN");
    }

    #[test]
    fn direct_autosomal_does_not_add_dispatcher_proband_gate() {
        let response = execute(
            &Fields::new(),
            &Fields::new(),
            "source",
            &fields(&[
                ("method", "auto"),
                ("sex", ""),
                ("proband.kind", "none"),
                ("father.kind", "tuple"),
                ("father.count", "1"),
                ("father.0.kind", "int"),
                ("father.0.value", "0"),
                ("mother.kind", "tuple"),
                ("mother.count", "1"),
                ("mother.0.kind", "int"),
                ("mother.0.value", "0"),
            ]),
        )
        .unwrap();
        assert_eq!(response["outcome.value"], "DENOVO");
        assert!(prefixed(&response, "state.").is_empty());
    }

    #[test]
    fn unknown_state_message_context_and_actual_ordinals_fail_closed() {
        let context = context();
        let state = initialize(&context).unwrap();
        let mut wrong = state.clone();
        wrong.insert("foreign".into(), "state".into());
        assert!(execute(&context, &wrong, "filter-begin", &Fields::new()).is_err());
        assert!(
            execute(
                &context,
                &state,
                "filter-begin",
                &fields(&[("ordinal", "1")])
            )
            .is_err()
        );
        assert!(execute(&Fields::new(), &state, "genotype", &Fields::new()).is_err());
        assert!(execute(&context, &state, "filter-reply", &fields(&[("ack", "1")])).is_err());
        let mut wrong = context.clone();
        wrong.insert("father.original_kind".into(), "other".into());
        assert!(initialize(&wrong).is_err());
    }
}
