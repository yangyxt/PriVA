//! Column preparation and decisions around the unchanged mature binomial calls.
//! No probability distribution is evaluated here. Unsupported dtype/shape or
//! warning behavior returns an explicit mixed-operation request to the host.
use super::columns::{aligned, boolean_cells, boolean_column, numeric_cells, numeric_column};
use super::values::{Cell, Result, Value};
use polars_core::chunked_array::ops::arity;
use polars_core::prelude::*;

fn prototype<'a>(left: &'a Value, right: &'a Value) -> &'a Value {
    if left.kind == "scalar" || right.kind == "series" && left.kind != "series" {
        right
    } else {
        left
    }
}

fn integral(value: &Value) -> bool {
    matches!(value.dtype.as_str(), "int64" | "Int64")
        || value.kind == "scalar" && matches!(value.cells.as_slice(), [Cell::Int(_)])
}

fn integers(cells: &[Cell]) -> Result<Int64Chunked> {
    let cells = cells.iter().map(|cell| match cell {
        Cell::Int(value) => Ok(Some(*value)),
        Cell::Na | Cell::None => Ok(None),
        _ => Err("noninteger native column".into()),
    }).collect::<Result<Vec<_>>>()?;
    Ok(Int64Chunked::from_iter_options("integer".into(), cells.into_iter()))
}

fn arithmetic(op: &str, left: &Value, right: &Value) -> Result<Value> {
    let (index, a, b) = aligned(left, right)?;
    let proto = prototype(left, right);
    if integral(left) && integral(right) && op != "div" {
        let a = integers(&a)?;
        let b = integers(&b)?;
        let result: Int64Chunked = arity::binary_elementwise(&a, &b,
            |a: Option<i64>, b: Option<i64>| a.zip(b).map(|(a,b)| match op {
                "add" => a.wrapping_add(b), "sub" => a.wrapping_sub(b),
                "mul" => a.wrapping_mul(b), "minimum" => a.min(b),
                "maximum" => a.max(b), _ => unreachable!(),
            }));
        let nullable = left.dtype == "Int64" || right.dtype == "Int64";
        let mut result = proto.column(result.into_iter().map(|cell|
            cell.map(Cell::Int).unwrap_or(Cell::Na)).collect(),
            if nullable { "Int64" } else { "int64" })?;
        result.index = index;
        return Ok(result);
    }
    if a.iter().chain(&b).any(|cell| !cell.missing() && cell.number().is_none()) {
        return Err("object/string arithmetic remains mature at this source site".into());
    }
    // Scalar Python zero division differs from vector NumPy/Pandas division.
    if op == "div" && left.kind == "scalar" && right.kind == "scalar"
        && b[0].number() == Some(0.0) {
        return Err("scalar zero division retains original exception".into());
    }
    let a = numeric_column(&a)?;
    let b = numeric_column(&b)?;
    let result: Float64Chunked = match op {
        "add" => &a + &b, "sub" => &a - &b, "mul" => &a * &b, "div" => &a / &b,
        "minimum" | "maximum" => arity::binary_elementwise(&a, &b,
            |a: Option<f64>, b: Option<f64>| a.zip(b).map(|(a,b)| {
                if a.is_nan() { a } else if b.is_nan() { b }
                else if op == "minimum" { if a <= b { a } else { b } }
                else if a >= b { a } else { b }
            })),
        _ => return Err("unknown arithmetic operation".into()),
    };
    let nullable = matches!(left.dtype.as_str(), "Float64" | "Int64")
        || matches!(right.dtype.as_str(), "Float64" | "Int64");
    let mut value = proto.column(numeric_cells(&result, nullable),
        if nullable { "Float64" } else { "float64" })?;
    value.index = index;
    Ok(value)
}

fn cast(value: &Value, dtype: &str) -> Result<Value> {
    match dtype {
        "int64" => {
            if value.cells.iter().any(|cell| cell.missing()
                || cell.number().is_none_or(|n| !n.is_finite()
                    || !(-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&n))) {
                return Err("invalid integer cast/warning remains mature".into());
            }
            let column = numeric_column(&value.cells)?;
            let result: Int64Chunked = arity::unary_elementwise(&column,
                |n: Option<f64>| n.map(|n| n as i64));
            value.column(result.into_no_null_iter().map(Cell::Int).collect(), "int64")
        }
        "bool" => {
            if value.cells.iter().any(|cell| matches!(cell, Cell::Na)) {
                return Err("pandas.NA truth conversion remains mature".into());
            }
            let cells = value.cells.iter().map(|cell| cell.truth().map(|v| Cell::Bool(v.unwrap())))
                .collect::<Result<Vec<_>>>()?;
            let column = boolean_column(&cells)?;
            value.column(boolean_cells(&column, false, false), "bool")
        }
        "boolean" => {
            if value.cells.iter().any(|cell| !cell.missing() && !matches!(cell, Cell::Bool(_))) {
                return Err("nonboolean extension conversion remains mature".into());
            }
            value.column(boolean_cells(&boolean_column(&value.cells)?, true, false), "boolean")
        }
        _ => Err("unrepresented source cast".into()),
    }
}

fn conditional(args: &[Value]) -> Result<Value> {
    if args.len() != 3 { return Err("where arity".into()); }
    let size = args.iter().filter(|a| a.kind != "scalar").map(|a| a.cells.len())
        .next().unwrap_or(1);
    if args.iter().any(|a| a.kind != "scalar" && a.cells.len() != size) {
        return Err("where broadcast shape requires mature NumPy".into());
    }
    let expand = |value: &Value| if value.kind == "scalar" {
        vec![value.cells[0].clone(); size]
    } else { value.cells.clone() };
    let mask = expand(&args[0]);
    if mask.iter().any(|cell| !matches!(cell, Cell::Bool(_))) {
        return Err("where nonboolean/missing mask remains mature".into());
    }
    let a = expand(&args[1]);
    let b = expand(&args[2]);
    let mask = boolean_column(&mask)?;
    let cells: Vec<_> = mask.into_iter().zip(a.into_iter().zip(b))
        .map(|(selected,(a,b))| if selected == Some(true) { a } else { b }).collect();
    let extension = args[1..].iter().any(|v| matches!(v.dtype.as_str(), "boolean" | "Int64" | "Float64" | "object"));
    let dtype = if extension { "object" }
        else if cells.iter().all(|c| matches!(c, Cell::Bool(_))) { "bool" }
        else if args[1..].iter().any(|v| v.dtype == "float64" || v.cells.iter().any(|c| matches!(c, Cell::Float(_)))) { "float64" }
        else if cells.iter().all(|c| matches!(c, Cell::Int(_))) { "int64" }
        else { return Err("unrepresented where dtype promotion".into()); };
    Ok(Value {kind:"array".into(), dtype:dtype.into(), index:Vec::new(), cells})
}

fn clip(args: &[Value]) -> Result<Value> {
    if args.len()!=3 { return Err("clip arity".into()); }
    let first=&args[0];
    let bound=|v:&Value| -> Result<Option<f64>> {
        if v.kind != "scalar" { return Err("vector clip remains mature".into()); }
        if matches!(v.cells[0],Cell::None) { Ok(None) }
        else { v.cells[0].number().map(Some).ok_or("nonnumeric clip bound".into()) }
    };
    let lower=bound(&args[1])?; let upper=bound(&args[2])?;
    if lower.is_some_and(f64::is_nan) || upper.is_some_and(f64::is_nan) {
        return Err("NaN clip bounds remain mature".into());
    }
    if integral(first) {
        if [lower,upper].into_iter().flatten().any(|v| v.fract()!=0.0 || v.abs()>9_007_199_254_740_992.0) {
            return Err("integer clip promotion remains mature".into());
        }
        let result:Int64Chunked=arity::unary_elementwise(&integers(&first.cells)?,|v:Option<i64>|
            v.map(|mut n| {if let Some(l)=lower {n=n.max(l as i64)} if let Some(u)=upper {n=n.min(u as i64)} n}));
        return first.column(result.into_iter().map(|n|n.map(Cell::Int).unwrap_or(Cell::Na)).collect(),&first.dtype);
    }
    let result:Float64Chunked=arity::unary_elementwise(&numeric_column(&first.cells)?, |v:Option<f64>|
        v.map(|mut n| {if let Some(l)=lower {if n<l {n=l}} if let Some(u)=upper {if n>u {n=u}} n}));
    first.column(numeric_cells(&result,first.dtype=="Float64"),&first.dtype)
}

pub fn calculate(op:&str,args:&[Value])->Option<Result<Value>> {
    let known=["add","sub","mul","div","minimum","maximum","round","cast_int64","cast_bool","cast_boolean","where","clip","isnan","replace_zero_nan","mask_missing","mask_value"];
    if !known.contains(&op) {return None;}
    Some((|| {
        let first=args.first().ok_or("missing binomial-preparation operand")?;
        match op {
            "add"|"sub"|"mul"|"div"|"minimum"|"maximum" => {
                if args.len()!=2 {return Err("arithmetic arity".into());}
                arithmetic(op,first,&args[1])
            }
            "round" => {
                if integral(first) {return Ok(first.clone());}
                let column=numeric_column(&first.cells)?;
                let result:Float64Chunked=arity::unary_elementwise(&column,|n:Option<f64>|n.map(f64::round_ties_even));
                first.column(numeric_cells(&result,first.dtype=="Float64"),&first.dtype)
            }
            "cast_int64"=>cast(first,"int64"), "cast_bool"=>cast(first,"bool"), "cast_boolean"=>cast(first,"boolean"),
            "where"=>conditional(args), "clip"=>clip(args),
            "isnan"=> {
                if first.cells.iter().any(|c| !c.missing() && c.number().is_none()) {return Err("object isnan remains mature".into());}
                let cells=first.cells.iter().map(|c|match c {Cell::Na=>Cell::Na, _=>Cell::Bool(c.missing())}).collect();
                first.column(cells,if matches!(first.dtype.as_str(),"Float64"|"Int64"|"boolean") {"boolean"} else {"bool"})
            }
            "replace_zero_nan"=> {
                let result:Float64Chunked=arity::unary_elementwise(&numeric_column(&first.cells)?,|n:Option<f64>|n.filter(|n|*n!=0.0));
                first.column(numeric_cells(&result,first.dtype=="Float64"),if first.dtype=="Float64" {"Float64"} else {"float64"})
            }
            "mask_missing"|"mask_value"=> {
                if args.len()!=(if op=="mask_missing" {2}else{3}) {return Err("mask arity".into());}
                let (_,_,mask)=aligned(first,&args[1])?;
                if mask.iter().any(|c|!matches!(c,Cell::Bool(_))) {return Err("missing mask requires mature semantics".into());}
                let other=if op=="mask_value" {
                    let (_,_,v)=aligned(first,&args[2])?;v
                } else {vec![if first.dtype=="boolean" {Cell::Na}else{Cell::Float(f64::NAN)};first.cells.len()]};
                let cells=boolean_column(&mask)?.into_iter().zip(first.cells.iter().zip(other))
                    .map(|(mask,(a,b))|if mask==Some(true){b}else{a.clone()}).collect();
                if !matches!(first.dtype.as_str(),"float64"|"Float64"|"boolean") {return Err("mask dtype promotion remains mature".into());}
                first.column(cells,&first.dtype)
            }
            _=>unreachable!(),
        }
    })())
}
