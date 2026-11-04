use super::*;
use crate::python_args::parse_long_option_occurrences;
use std::collections::HashMap;
const HELP: &str = "usage: priva installer-config-paths {resolve,managed-roots,portable-updates,resolve-cache-paths,defaults,resource-groups,adoption-plan,adopt-sources,render,bundled-resource,reverse-chain-source,source-format,source-needed,cache-identity,sha256,gnomad-members,gnomad-validate,hash-order}\n\
--config YAML --source-root CHECKOUT [--helper-root CHECKOUT] [--check] [--output PATH]\n\
resolve [--managed-updates-json JSON]; library plans [--updates-json JSON]; adopt-sources --plan-json JSON --cache-root PATH\n\
scalar procedures --key NAME --value-json JSON; file procedures --path FILE\n\
Use an explicit PYTHONHASHSEED profile to preserve the original replacement-table iteration order.";
fn required<'a>(options: &'a HashMap<String, String>, key: &str) -> Result<&'a str> {
    options
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| invalid(format!("missing --{key}")))
}
fn paths(text: &str) -> Result<Paths> {
    let v = json_text::parse(text, false).map_err(|e| error("JSONDecodeError", e))?;
    v.as_object()
        .ok_or_else(|| error("TypeError", "path updates must be a mapping"))?
        .iter()
        .map(|(k, v)| {
            v.as_str()
                .map(|v| {
                    (
                        crate::python_unicode::encode(k),
                        crate::python_unicode::encode(v),
                    )
                })
                .ok_or_else(|| error("TypeError", format!("{k} must be a path string")))
        })
        .collect()
}
fn print(value: &Json, sort: bool) -> Result<()> {
    let mut value = value.clone();
    if sort {
        if let Json::Object(values) = &mut value {
            values.sort_by(|a, b| a.0.cmp(&b.0));
        }
    }
    println!(
        "{}",
        receipt_json::RecordJson::from_json(&value).dumps(true)?
    );
    Ok(())
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let commands = [
        "resolve",
        "managed-roots",
        "portable-updates",
        "resolve-cache-paths",
        "defaults",
        "resource-groups",
        "adoption-plan",
        "adopt-sources",
        "render",
        "bundled-resource",
        "reverse-chain-source",
        "source-format",
        "source-needed",
        "cache-identity",
        "sha256",
        "gnomad-members",
        "gnomad-validate",
        "hash-order",
    ];
    if !commands.contains(&command) {
        return Err(invalid(HELP));
    }
    let Some(pairs) = parse_long_option_occurrences(
        &args[1..],
        &[
            "config",
            "source-root",
            "helper-root",
            "managed-updates-json",
            "output",
            "updates-json",
            "plan-json",
            "cache-root",
            "key",
            "value-json",
            "path",
            "assembly",
            "reference",
        ],
        &["check", "no-source-root"],
    )
    .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let options: HashMap<_, _> = pairs
        .into_iter()
        .map(|(k, v)| {
            let value = if [
                "config",
                "source-root",
                "helper-root",
                "output",
                "cache-root",
                "path",
                "reference",
                "key",
            ]
            .contains(&k.as_str())
            {
                crate::python_unicode::encode(&v)
            } else {
                v
            };
            (k, value)
        })
        .collect();
    for key in match command {
        "hash-order" => vec![],
        "adopt-sources" => vec!["plan-json", "cache-root"],
        "bundled-resource" => vec!["key", "value-json"],
        "reverse-chain-source" => vec!["value-json"],
        "source-format" => vec!["key", "path"],
        "source-needed" => vec!["key", "config"],
        "cache-identity" | "sha256" | "gnomad-members" => vec!["path"],
        "gnomad-validate" => vec!["path", "assembly"],
        "render" => vec!["config", "updates-json"],
        _ => vec!["config", "source-root"],
    } {
        required(&options, key)?;
    }
    crate::builder_error::finish((|| {
        let root = options.get("source-root").map(|v| path_value(v));
        let helper = options
            .get("helper-root")
            .map(|v| path_value(v))
            .or_else(|| root.clone())
            .unwrap_or_else(|| PathBuf::from("."));
        let ctx = Context::new(&helper)?;
        if command == "hash-order" {
            let order =
                crate::python_set_order::StringSet::literal(ctx.seed, REPLACEMENT).into_order();
            let difference = crate::python_set_order::StringSet::literal(ctx.seed, REPLACEMENT)
                .difference_small(&["hgnc_table"]);
            return print(
                &Json::Object(vec![
                    (
                        "literal".into(),
                        Json::Array(order.into_iter().map(Json::Str).collect()),
                    ),
                    (
                        "difference".into(),
                        Json::Array(difference.into_iter().map(Json::Str).collect()),
                    ),
                ]),
                false,
            );
        }
        let root = if options.contains_key("no-source-root") {
            None
        } else {
            root.as_deref()
        };
        let file = options.get("path").map(|v| path_value(v));
        match command {
            "cache-identity" => {
                return print(&files::cache_identity(file.as_deref().unwrap())?, false);
            }
            "sha256" => {
                println!("{}", files::sha256_file(file.as_deref().unwrap())?);
                return Ok(());
            }
            "gnomad-members" => {
                return print(
                    &Json::Object(
                        gnomad::family_paths(file.as_deref().unwrap())?
                            .into_iter()
                            .map(|(k, v)| (k, Json::Str(v.display().to_string())))
                            .collect(),
                    ),
                    false,
                );
            }
            "gnomad-validate" => {
                let r = options
                    .get("reference")
                    .map_or(Value::Missing, |v| Value::Text(v));
                return print(
                    &gnomad::validate_family(
                        file.as_deref().unwrap(),
                        required(&options, "assembly")?,
                        r,
                    )?,
                    false,
                );
            }
            "source-format" => {
                return print(
                    &plan::source_format(
                        &ctx,
                        required(&options, "key")?,
                        required(&options, "path")?,
                    )?
                    .map_or(Json::Null, |v| Json::Str(v.into())),
                    false,
                );
            }
            "bundled-resource" | "reverse-chain-source" => {
                let c = Config::parse(&format!("value: {}", required(&options, "value-json")?))?;
                let v = c.get("value");
                let p = if command == "reverse-chain-source" {
                    Some(plan::reverse_chain_source(&ctx, v, root)?)
                } else {
                    Some(plan::bundled_resource(
                        &ctx,
                        required(&options, "key")?,
                        v.string()
                            .ok_or_else(|| error("TypeError", "value must be a path string"))?,
                        root,
                    )?)
                    .flatten()
                };
                return print(
                    &p.map_or(Json::Null, |p| Json::Str(p.display().to_string())),
                    false,
                );
            }
            "adopt-sources" => {
                let v = json_text::parse(required(&options, "plan-json")?, false)
                    .map_err(|e| error("JSONDecodeError", e))?;
                let Json::Array(items) = v else {
                    return Err(error("TypeError", "adoption plan must be a list"));
                };
                let mut plan = vec![];
                for item in items {
                    let getstr = |k| {
                        item.get(k)
                            .and_then(Json::as_str)
                            .ok_or_else(|| error("KeyError", repr(k)))
                    };
                    plan.push(plan::Adoption {
                        key: getstr("key")?.into(),
                        source: path_value(&crate::python_unicode::encode(getstr("source")?)),
                        destination: path_value(&crate::python_unicode::encode(getstr(
                            "destination",
                        )?)),
                        kind: getstr("kind")?.into(),
                    });
                }
                adopt::adopt_sources(&ctx, &plan, Path::new(required(&options, "cache-root")?))?;
                return Ok(());
            }
            _ => {}
        }
        let path = path_value(required(&options, "config")?);
        if command == "resolve" {
            if let Some(v) = adopt::resolve_config(
                &ctx,
                &path,
                root.ok_or_else(|| invalid("missing --source-root"))?,
                options.contains_key("check"),
                options.get("managed-updates-json").map(String::as_str),
                options.get("output").map(Path::new),
            )? {
                print(&pathmap(&v), true)?;
            }
            return Ok(());
        }
        let original = files::read_text(&path)?;
        if command == "render" {
            let updated =
                adopt::render_config(&original, &paths(required(&options, "updates-json")?)?)?;
            print!("{updated}");
            return Ok(());
        }
        let config = Config::parse(&original)?;
        if command == "managed-roots" {
            if let Some(v) = adopt::configure_cache_paths(
                &ctx,
                &path,
                &config,
                if options.contains_key("check") {
                    "--check"
                } else {
                    "write"
                },
                root,
            )? {
                print(&pathmap(&v), true)?;
            }
            return Ok(());
        }
        if command == "source-needed" {
            return print(
                &Json::Bool(plan::source_is_needed(required(&options, "key")?, &config)),
                false,
            );
        }
        let value = match command {
            "resolve-cache-paths" => {
                let (d, u) = plan::resolve_cache_paths(&ctx, &config, root)?;
                Json::Array(vec![pathmap(&d), pathmap(&u)])
            }
            "portable-updates" => {
                pathmap(&defaults::portable_updates(&ctx, &config, root.unwrap())?)
            }
            "defaults" => pathmap(&defaults::portable_defaults(&config, root.unwrap())?),
            "resource-groups" => {
                let dirs = paths(required(&options, "updates-json")?)?;
                pathmap(&plan::resource_groups(&ctx, &config, &dirs)?)
            }
            "adoption-plan" => {
                let updates = if let Some(v) = options.get("updates-json") {
                    paths(v)?
                } else {
                    defaults::portable_updates(&ctx, &config, root.unwrap())?
                };
                Json::Array(
                    plan::adoption_plan(&ctx, &config, &updates, root.unwrap())?
                        .iter()
                        .map(plan::Adoption::json)
                        .collect(),
                )
            }
            _ => unreachable!(),
        };
        print(&value, false)
    })())
}
