//! Pure control-program checks. Scripted effect replies are deliberately not
//! claims about actual Python packages or equality to original procedures.

use super::*;
use crate::frequency_criteria::Cell;

fn begin(mode: Mode) -> (Program, Action) {
    let mut program = Program::new(Procedure::GnomadRareAf, mode, 101, 102).unwrap();
    let action = program.start().unwrap();
    (program, action)
}

fn respond(program: &mut Program, action: &Action, reply: Reply) -> Action {
    let Action::Effect { token, .. } = action else { panic!("effect required"); };
    program.resume(*token, reply).unwrap()
}

fn first_kernel(mode: Mode) -> (Program, Action) {
    let (mut program, action) = begin(mode);
    let mut action = respond(&mut program, &action, Reply::Truth(true));
    for handle in 103..106 { action = respond(&mut program, &action, Reply::Object(handle)); }
    assert!(matches!(action, Action::Effect { effect: Effect::PrepareKernel(_), .. }));
    (program, action)
}

fn trace_route(is_float: bool, is_string: bool, cutoff: f64, fail_getter: Option<usize>) -> (Program, Vec<Effect>, Action) {
    let (mut program, mut action) = begin(Mode::Mixed);
    let mut trace = Vec::new();
    let mut next_handle = 200;
    let mut getters = 0;
    for _ in 0..200 {
        let Action::Effect { effect, .. } = &action else { return (program, trace, action); };
        trace.push(effect.clone());
        let reply = match effect {
            Effect::IsInstance { class: TypeTest::Float, .. } => Reply::Truth(is_float),
            Effect::IsInstance { class: TypeTest::String, .. } => Reply::Truth(is_string),
            Effect::InspectBuiltinFloat(_) => Reply::BuiltinFloat(cutoff.to_bits()),
            Effect::PrepareKernel(_) => Reply::NotRepresented("scripted effect frontier".into()),
            Effect::GetItem { .. } => {
                getters += 1;
                if fail_getter == Some(getters) { Reply::Raised(999) }
                else { next_handle += 1; Reply::Object(next_handle) }
            }
            Effect::BindNumpyWhere | Effect::GetAttribute { .. } | Effect::InvokeMature { .. }
                | Effect::Materialize { .. } => { next_handle += 1; Reply::Object(next_handle) }
        };
        action = respond(&mut program, &action, reply);
    }
    panic!("native procedure did not terminate within its finite tape");
}

fn getter_keys(trace: &[Effect]) -> Vec<Object> {
    trace.iter().filter_map(|effect| match effect {
        Effect::GetItem { key, .. } => Some(key.clone()), _ => None,
    }).collect()
}

#[test]
fn unimplemented_procedures_are_not_advertised() {
    assert_eq!(Procedure::named("gnomAD_rare_AF").unwrap().name(), "gnomAD_rare_AF");
    assert!(Procedure::named("BA1_criteria").is_err());
    assert!(Procedure::named("ACMG_criteria_assign").is_err());
}

#[test]
fn other_cutoff_type_returns_none_without_frame_reads() {
    let (_program, trace, result) = trace_route(false, false, 0.0, None);
    assert_eq!(result, Action::Return(Object::None));
    assert_eq!(trace.len(), 2);
    assert!(getter_keys(&trace).is_empty());
}

#[test]
fn float_branch_captures_where_before_eager_operands() {
    let (program, trace, result) = trace_route(true, false, 0.001, None);
    assert!(matches!(result, Action::Return(Object::Handle(_))));
    assert!(matches!(trace[1], Effect::BindNumpyWhere));
    assert_eq!(getter_keys(&trace), vec![Object::Text("gnomAD_joint_AN_max"),
        Object::Text("gnomAD_joint_AF_max"), Object::Text("gnomAD_joint_AF")]);
    let selections: Vec<_> = trace.iter().filter_map(|effect| match effect {
        Effect::PrepareKernel(call) if call.kernel == Kernel::Where => Some(call), _ => None,
    }).collect();
    assert_eq!(selections.len(), 1);
    assert!(matches!(&selections[0].invocation,
        Invocation::Bound { callable: Object::Handle(201), args } if args.len() == 3));
    assert_eq!(program.native_operations, 1); // exact scalar reciprocal
}

#[test]
fn string_branch_repeats_every_cutoff_getter() {
    let (_program, trace, result) = trace_route(false, true, 0.0, None);
    assert!(matches!(result, Action::Return(Object::Handle(_))));
    assert_eq!(getter_keys(&trace), vec![Object::Text("gnomAD_joint_AN_max"), Object::Handle(102),
        Object::Text("gnomAD_joint_AF_max"), Object::Handle(102),
        Object::Text("gnomAD_joint_AF"), Object::Handle(102)]);
    let fills: Vec<_> = trace.iter().filter_map(|effect| match effect {
        Effect::PrepareKernel(call) if call.kernel == Kernel::FillNa => Some(call.operands[1].clone()), _ => None,
    }).collect();
    assert_eq!(fills, vec![Object::Integer(1_000_000), Object::Float(0x3ee4_f8b5_88e3_68f1),
        Object::Integer(0), Object::Integer(0), Object::Integer(0), Object::Integer(0)]);
}

#[test]
fn zero_cutoff_raises_after_an_fill_before_af_getters() {
    for zero in [0.0_f64, -0.0_f64] {
        let (_program, trace, result) = trace_route(true, false, zero, None);
        assert_eq!(result, Action::Raise { class: "ZeroDivisionError", message: "float division by zero" });
        assert_eq!(getter_keys(&trace), vec![Object::Text("gnomAD_joint_AN_max")]);
        assert!(trace.iter().any(|effect| matches!(effect, Effect::InvokeMature {
            invocation: Invocation::Bound { args, .. }, .. } if args == &vec![Object::Integer(1_000_000)])));
    }
}

#[test]
fn getter_failure_retains_exact_prefix_and_exception_handle() {
    let (_program, trace, result) = trace_route(false, true, 0.0, Some(4));
    assert_eq!(result, Action::HostRaised { exception: 999 });
    assert_eq!(getter_keys(&trace), vec![Object::Text("gnomAD_joint_AN_max"), Object::Handle(102),
        Object::Text("gnomAD_joint_AF_max"), Object::Handle(102)]);
    assert!(!trace.iter().any(|effect| matches!(effect,
        Effect::PrepareKernel(call) if call.kernel == Kernel::Where)));
}

#[test]
fn float_reciprocal_is_native_and_retains_negative_value() {
    let (program, trace, _result) = trace_route(true, false, -0.125, None);
    let comparisons: Vec<_> = trace.iter().filter_map(|effect| match effect {
        Effect::PrepareKernel(call) if call.kernel == Kernel::GreaterEqual => Some(call), _ => None,
    }).collect();
    assert_eq!(comparisons[0].operands[1], Object::Float((-8.0_f64).to_bits()));
    assert_eq!(program.native_operations, 1);
}

#[test]
fn nonfinite_scalar_is_an_explicit_mature_effect() {
    let (_program, trace, _result) = trace_route(true, false, f64::NAN, None);
    assert_eq!(trace.iter().filter(|effect| matches!(effect, Effect::InvokeMature {
        invocation: Invocation::Binary { operation: Kernel::Divide,
            left: Object::Integer(1), right: Object::Handle(102) }, .. })).count(), 1);
}

#[test]
fn native_kernel_result_is_consumed_by_boxing() {
    let (mut program, action) = first_kernel(Mode::NativeMath);
    let column = Value { kind: "series".into(), dtype: "float64".into(), index: vec!["srow".into()],
        cells: vec![Cell::Float(f64::NAN)] };
    let action = respond(&mut program, &action, Reply::Prepared(vec![column,
        Value::scalar(Cell::Int(1_000_000))]));
    let Action::Effect { effect: Effect::Materialize { value, .. }, .. } = &action else { panic!("native boxing required"); };
    assert_eq!(value.cells, vec![Cell::Float(1_000_000.0)]);
    assert_eq!(value.index, vec!["srow"]);
    assert_eq!(program.native_operations, 1);
    let next = respond(&mut program, &action, Reply::Object(300));
    assert!(matches!(next, Action::Effect { effect: Effect::InspectBuiltinFloat(Object::Handle(102)), .. }));
}

#[test]
fn native_math_unsupported_stops_without_mature_retry() {
    let (mut program, action) = first_kernel(Mode::NativeMath);
    let result = respond(&mut program, &action, Reply::NotRepresented("unknown extension".into()));
    assert!(matches!(result, Action::Unsupported { reason, .. } if reason == "unknown extension"));
    assert_eq!(program.mature_operations, 0);
    assert!(program.resume(4, Reply::Object(500)).is_err());
}

#[test]
fn mixed_mode_invokes_only_the_outstanding_operation_once() {
    let (mut program, action) = first_kernel(Mode::Mixed);
    let action = respond(&mut program, &action, Reply::NotRepresented("explicit package frontier".into()));
    assert!(matches!(&action, Action::Effect { effect: Effect::InvokeMature {
        invocation: Invocation::Bound { callable: Object::Handle(105), args }, .. }, .. }
        if args == &vec![Object::Integer(1_000_000)]));
    let result = respond(&mut program, &action, Reply::Raised(800));
    assert_eq!(result, Action::HostRaised { exception: 800 });
    assert_eq!(program.mature_operations, 1);
}

#[test]
fn malformed_operand_count_poisoned_before_math() {
    let (mut program, action) = first_kernel(Mode::Mixed);
    let Action::Effect { token, .. } = action else { unreachable!() };
    assert!(program.resume(token, Reply::Prepared(vec![])).is_err());
    assert!(program.resume(token, Reply::NotRepresented("retry".into())).is_err());
    assert_eq!(program.native_operations, 0);
    assert_eq!(program.mature_operations, 0);
}

#[test]
fn malformed_typed_shape_never_reaches_a_kernel_or_mature_call() {
    let (mut program, action) = first_kernel(Mode::Mixed);
    let Action::Effect { token, .. } = action else { unreachable!() };
    let broken = Value { kind: "scalar".into(), dtype: "scalar".into(), index: vec![], cells: vec![] };
    assert!(program.resume(token, Reply::Prepared(vec![broken, Value::scalar(Cell::Int(0))])).is_err());
    assert_eq!(program.native_operations, 0);
    assert_eq!(program.mature_operations, 0);
}

#[test]
fn wrong_reply_type_and_token_never_recover() {
    let (mut program, action) = begin(Mode::Mixed);
    let Action::Effect { token, .. } = action else { unreachable!() };
    assert!(program.resume(token, Reply::Object(200)).is_err());
    assert!(program.resume(token, Reply::Truth(true)).is_err());
    let (mut program, action) = begin(Mode::Mixed);
    let Action::Effect { token, .. } = action else { unreachable!() };
    assert!(program.resume(token + 1, Reply::Truth(true)).is_err());
    assert!(program.resume(token, Reply::Truth(true)).is_err());
}

#[test]
fn duplicate_token_and_second_start_are_rejected() {
    let (mut program, action) = begin(Mode::Mixed);
    let Action::Effect { token, .. } = action else { unreachable!() };
    let _next = program.resume(token, Reply::Truth(true)).unwrap();
    assert!(program.resume(token, Reply::Object(200)).is_err());
    assert!(program.start().is_err());
}

#[test]
fn late_materialization_failure_does_not_replay_native_math() {
    let (mut program, action) = first_kernel(Mode::Mixed);
    let values = vec![Value { kind: "array".into(), dtype: "float64".into(), index: vec![],
        cells: vec![Cell::Float(4.0)] }, Value::scalar(Cell::Int(1_000_000))];
    let action = respond(&mut program, &action, Reply::Prepared(values));
    let result = respond(&mut program, &action, Reply::Raised(803));
    assert_eq!(result, Action::HostRaised { exception: 803 });
    assert_eq!(program.native_operations, 1);
    assert_eq!(program.mature_operations, 0);
}

#[test]
fn completed_procedure_cannot_accept_another_reply() {
    let (mut program, _trace, result) = trace_route(true, false, 0.001, None);
    assert!(matches!(result, Action::Return(_)));
    assert!(program.resume(1000, Reply::Object(1001)).is_err());
    assert!(program.start().is_err());
}
