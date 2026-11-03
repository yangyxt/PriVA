//! Frequency-owned string, grouping, cache-map and criterion publication phases.
//! Their ordering comes from the pinned original source lowering, not a new
//! model or a fused shortcut through the original eager calculations.
use super::columns::{boolean_column, numeric_column};
use super::values::{Cell, Result, Value};
use polars_core::chunked_array::ops::arity;
use polars_core::prelude::*;
use std::collections::BTreeMap;

fn key(cell:&Cell)->Result<String> {
    match cell {
        Cell::Text(s)=>Ok(format!("s{s}")),
        Cell::Int(n)=>Ok(format!("i{:020}",i128::from(*n)+(1i128<<63))),
        _=>Err("dictionary/grouping key equality requires mature handling".into()),
    }
}

fn mapped(args:&[Value], default_zero:bool)->Result<Value> {
    if args.len()!=3 || args[1].cells.len()!=args[2].cells.len() {return Err("cache map arity".into());}
    let mut map=BTreeMap::new();
    for (k,v) in args[1].cells.iter().zip(&args[2].cells) {
        if map.insert(key(k)?,v.clone()).is_some(){return Err("mapping key collision".into());}
    }
    let cells=args[0].cells.iter().map(|c| {
        if c.missing(){return Ok(if default_zero {Cell::Int(0)}else{Cell::Float(f64::NAN)});}
        Ok(map.get(&key(c)?).cloned().unwrap_or(if default_zero {Cell::Int(0)}else{Cell::Float(f64::NAN)}))
    }).collect::<Result<Vec<_>>>()?;
    let dtype=if cells.iter().all(|c|matches!(c,Cell::Int(_))){"int64"}
        else if cells.iter().all(|c| c.number().is_some() || c.missing()) {"float64"}
        else {"object"};
    args[0].column(cells,dtype)
}

fn group_max(args:&[Value])->Result<Value> {
    if args.len()!=2 || args[0].cells.len()!=args[1].cells.len() || args[0].index!=args[1].index {
        return Err("groupby alignment remains mature".into());
    }
    let mut groups:BTreeMap<String,(Cell,Vec<Option<f64>>)>=BTreeMap::new();
    let mut key_type=None;
    for (v,k) in args[0].cells.iter().zip(&args[1].cells) {
        if k.missing(){continue;}
        let kind=matches!(k,Cell::Text(_));
        if key_type.is_some_and(|old|old!=kind){return Err("mixed group sort remains mature".into());}
        key_type=Some(kind);
        groups.entry(key(k)?).or_insert_with(||(k.clone(),Vec::new())).1
            .push(if v.missing(){None}else{Some(v.number().ok_or("nonnumeric group maximum")?)});
    }
    let mut cells=Vec::new();
    for (_,(k,v)) in groups {
        let column=Float64Chunked::from_iter_options("group".into(),v.into_iter());
        cells.push(k);cells.push(Cell::Float(column.max().unwrap_or(f64::NAN)));
    }
    Ok(Value{kind:"pairs".into(),dtype:"object".into(),index:Vec::new(),cells})
}

fn exon_map(args:&[Value])->Result<Value> {
    if args.len()!=5 || args[0].cells.len()!=args[1].cells.len()
        || args[0].index!=args[1].index || args[2].cells.len()!=args[3].cells.len()
        || args[2].cells.len()!=args[4].cells.len() {
        return Err("exon cache dimensions differ".into());
    }
    let mut mapping=BTreeMap::new();
    for ((feature,exon),value) in args[2].cells.iter().zip(&args[3].cells).zip(&args[4].cells) {
        if !matches!((feature,exon),(Cell::Text(_),Cell::Text(_))) {
            return Err("nontext exon cache key remains mature".into());
        }
        if mapping.insert((key(feature)?,key(exon)?),value.clone()).is_some() {
            return Err("duplicate exon cache key".into());
        }
    }
    let mut cells=Vec::new();
    for (feature,exon) in args[0].cells.iter().zip(&args[1].cells) {
        if !matches!((feature,exon),(Cell::Text(_),Cell::Text(_))) {
            return Err("nontext exon row key remains mature".into());
        }
        cells.push(mapping.get(&(key(feature)?,key(exon)?)).cloned().unwrap_or(Cell::Float(f64::NAN)));
    }
    let dtype=if cells.is_empty(){"float64"}
        else if cells.iter().all(|c|matches!(c,Cell::Int(_))){"int64"}
        else if cells.iter().all(|c|c.number().is_some()||c.missing()) {"float64"}
        else {return Err("exon result dtype inference remains mature".into());};
    args[0].column(cells,dtype)
}

pub fn calculate(op:&str,args:&[Value])->Option<Result<Value>> {
    let known=["lower","upper","removeprefix","isin","contains_literal","condition_token","cache_map_zero","cache_map_nan","exon_map","group_max_pairs","full","zeros_bool","assign_array","filter_none"];
    if !known.contains(&op){return None;}
    Some((|| {
        let first=args.first().ok_or("missing frequency-rule operand")?;
        match op {
            "lower"|"upper"|"removeprefix"|"contains_literal"|"condition_token"=> {
                let parameter=if args.len()>1 {args[1].text()?}else{""};
                let text=first.cells.iter().map(|c|match c {Cell::Text(s) if s.is_ascii()=>Ok(Some(s.as_str())),Cell::Na|Cell::None=>Ok(None),_=>Err("non-ASCII/custom string remains mature".into())}).collect::<Result<Vec<_>>>()?;
                if text.iter().any(Option::is_none) && first.dtype!="string" {return Err("object string missingness remains mature".into());}
                let column=StringChunked::from_iter_options("text".into(),text.into_iter());
                if op=="contains_literal" || op=="condition_token" {
                    let result:BooleanChunked=arity::unary_elementwise(&column,|s:Option<&str>|s.map(|s| {
                        if op=="contains_literal" {s.contains(parameter)} else {
                            let mut parts=s.split(';').peekable();
                            while let Some(p)=parts.next(){if p==parameter || parts.peek().is_none() && p.strip_suffix('\n')==Some(parameter){return true;}}
                            false
                        }
                    }));
                    return first.column(result.into_iter().map(|v|v.map(Cell::Bool).unwrap_or(Cell::Na)).collect(),if first.dtype=="string"{"boolean"}else{"bool"});
                }
                let result:StringChunked=arity::unary_elementwise(&column,|s:Option<&str>|s.map(|s|match op {
                    "lower"=>s.to_ascii_lowercase(),"upper"=>s.to_ascii_uppercase(),_=>s.strip_prefix(parameter).unwrap_or(s).to_owned(),
                }));
                first.column(result.into_iter().map(|v|v.map(|s|Cell::Text(s.into())).unwrap_or(Cell::Na)).collect(),&first.dtype)
            }
            "isin"=> {
                if args.len()!=2{return Err("isin arity".into());}
                // The host admits only string choices. None, pd.NA and NaN
                // operands cannot match such a choice and produce false.
                let choices=args[1].cells.iter().map(|c|match c {
                    Cell::Text(s)=>Ok(s.as_str()),
                    _=>Err("nontext isin choices remain mature".into()),
                }).collect::<Result<std::collections::BTreeSet<_>>>()?;
                let text=first.cells.iter().map(|c|match c {
                    Cell::Text(s)=>Ok(Some(s.as_str())),
                    c if c.missing()=>Ok(None),
                    _=>Err("nontext isin remains mature".into()),
                }).collect::<Result<Vec<_>>>()?;
                let column=StringChunked::from_iter_options("membership".into(),text.into_iter());
                let result:BooleanChunked=arity::unary_elementwise(&column,|s:Option<&str>| {
                    Some(s.is_some_and(|s|choices.contains(s)))
                });
                first.column(result.into_iter().map(|v|Cell::Bool(v.expect("membership is non-null"))).collect(),"bool")
            }
            "cache_map_zero"=>mapped(args,true),"cache_map_nan"=>mapped(args,false),
            "exon_map"=>exon_map(args),
            "group_max_pairs"=>group_max(args),
            "filter_none"=>first.column(first.cells.iter().map(|c|Cell::Bool(!matches!(c,Cell::None))).collect(),"bool"),
            "full"|"zeros_bool"=> {
                let [Cell::Int(n)]=first.cells.as_slice() else{return Err("array size integer required".into());};
                let n=usize::try_from(*n).map_err(|_|"negative array size")?;
                if n>16_777_216{return Err("bounded array size exceeded".into());}
                let (value,dtype)=if op=="zeros_bool"{(Cell::Bool(false),"bool")}else{
                    if args.len()!=2 || args[1].kind!="scalar"{return Err("full scalar required".into());}
                    let v=args[1].cells[0].clone();let dtype=if matches!(v,Cell::Float(_)){"float64"}else if matches!(v,Cell::Int(_)){"int64"}else{return Err("full dtype remains mature".into());};(v,dtype)
                };
                Ok(Value{kind:"array".into(),dtype:dtype.into(),index:Vec::new(),cells:vec![value;n]})
            }
            "assign_array"=> {
                if args.len()!=3 || first.kind!="array" || args[1].cells.len()!=first.cells.len(){return Err("array assignment shape remains mature".into());}
                let mask=boolean_column(&args[1].cells)?;
                if mask.null_count()!=0{return Err("missing assignment mask remains mature".into());}
                let selected=mask.sum().unwrap_or(0) as usize;
                if args[2].kind!="scalar" && args[2].cells.len()!=selected{return Err("assignment RHS broadcasting remains mature".into());}
                let mut rhs=0;let mut cells=first.cells.clone();
                for (i,m) in mask.into_iter().enumerate(){if m==Some(true){cells[i]=args[2].cells[if args[2].kind=="scalar"{0}else{rhs}].clone();rhs+=1;}}
                if first.dtype=="float64"{let _=numeric_column(&cells)?;}
                else if first.dtype=="int64" && cells.iter().any(|c|!matches!(c,Cell::Int(_))) {return Err("assignment dtype cast remains mature".into());}
                else if first.dtype=="bool" && cells.iter().any(|c|!matches!(c,Cell::Bool(_))) {return Err("Boolean assignment cast remains mature".into());}
                first.column(cells,&first.dtype)
            }
            _=>unreachable!(),
        }
    })())
}
