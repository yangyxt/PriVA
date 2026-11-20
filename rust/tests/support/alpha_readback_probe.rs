//! Explicit, ignored readback probe for an already prepared Parquet file.
//! No input is modified. The caller owns source/input pins and timing receipts.
//! Select with the fully qualified Cargo test name and PRIVA_ALPHA_PROBE_*:
//! INPUT is a Parquet path, GROUPS is comma-separated zero-based row groups,
//! MODE is frames, rows, digest, or scan; STRATEGY is none or rowgroups.
//! No production path uses these variables.
use super::*;
use std::path::PathBuf;
use std::time::Instant;

#[test]
#[ignore = "requires an explicitly pinned prepared cache and selected row groups"]
fn bounded_readback() -> Result<()> {
    let input = PathBuf::from(std::env::var("PRIVA_ALPHA_PROBE_INPUT")?);
    let groups = std::env::var("PRIVA_ALPHA_PROBE_GROUPS")?
        .split(',')
        .map(str::parse::<usize>)
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let mode = std::env::var("PRIVA_ALPHA_PROBE_MODE")?;
    let strategy = std::env::var("PRIVA_ALPHA_PROBE_STRATEGY").unwrap_or_else(|_| "none".into());
    assert!(matches!(
        mode.as_str(),
        "frames" | "rows" | "digest" | "scan"
    ));
    assert!(matches!(strategy.as_str(), "none" | "rowgroups"));
    assert!(!groups.is_empty() && groups.windows(2).all(|pair| pair[0] < pair[1]));
    let started = Instant::now();
    let header = metadata(&input)?;
    assert_eq!(header.get("priva.schema").map(String::as_str), Some(SCHEMA));
    let mut reader = ParquetReader::new(File::open(&input)?);
    let meta = reader.get_metadata()?.clone();
    assert!(groups.iter().all(|&group| group < meta.row_groups.len()));
    let metadata_s = started.elapsed().as_secs_f64();
    if mode == "scan" {
        // Exercise the production scanner, not the instrumented phase loop.
        assert_eq!(groups.len(), meta.row_groups.len());
        assert_eq!(strategy, "rowgroups");
        let mut sha = Sha256::new();
        let mut rows = 0;
        scan(&input, |value| {
            digest(&value, &mut sha);
            rows += 1;
            Ok(())
        })?;
        assert_eq!(rows, meta.num_rows);
        println!(
            "ALPHA_READBACK_PROBE {{\"mode\":\"scan\",\"strategy\":\"production-rowgroups\",\"groups\":{},\"rows\":{rows},\"total_s\":{},\"digest\":\"{}\"}}",
            groups.len(),
            started.elapsed().as_secs_f64(),
            sha.finish_hex()
        );
        return Ok(());
    }
    let mut offset = 0;
    let mut rows = 0;
    let mut read_s = 0.0;
    let mut visit_s = 0.0;
    let mut sha = Sha256::new();
    for (index, group) in meta.row_groups.iter().enumerate() {
        let count = group.num_rows();
        if groups.binary_search(&index).is_ok() {
            let before = Instant::now();
            let mut reader = ParquetReader::new(File::open(&input)?)
                .with_slice(Some((offset, count)))
                .read_parallel(if strategy == "none" {
                    ParallelStrategy::None
                } else {
                    ParallelStrategy::RowGroups
                });
            reader.set_metadata(meta.clone());
            let frame = reader.finish()?;
            assert_eq!(frame.schema().as_ref(), frame_schema()?.as_ref());
            assert_eq!(frame.height(), count);
            read_s += before.elapsed().as_secs_f64();
            let before = Instant::now();
            if mode != "frames" {
                super::rows(&frame, |value| {
                    if mode == "digest" {
                        digest(&value, &mut sha);
                    } else {
                        std::hint::black_box(value);
                    }
                    Ok(())
                })?;
            }
            visit_s += before.elapsed().as_secs_f64();
            rows += count;
        }
        offset += count;
    }
    println!(
        "ALPHA_READBACK_PROBE {{\"mode\":\"{mode}\",\"strategy\":\"{strategy}\",\"groups\":{},\"rows\":{rows},\"metadata_s\":{metadata_s},\"read_s\":{read_s},\"visit_s\":{visit_s},\"total_s\":{},\"digest\":\"{}\"}}",
        groups.len(),
        started.elapsed().as_secs_f64(),
        sha.finish_hex()
    );
    Ok(())
}

#[test]
fn ordered_groups_preserve_all_typed_digest_fields() -> Result<()> {
    use crate::alpha_domain::key::Key;
    use crate::vcf_variant::Scalar;
    let folder = std::env::temp_dir().join(format!(
        "priva-alpha-readback-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)?
            .as_nanos()
    ));
    std::fs::create_dir(&folder)?;
    let path = folder.join("groups.parquet");
    let rows: Vec<_> = (0..9)
        .map(|i| Row {
            gene: format!("G{i}"),
            gene_order: i,
            path: vec!["domain".into(), format!("层{i}")],
            member_order: 8 - i,
            item_order: i * 3,
            key: format!("key{i}"),
            kind: "fixture".into(),
            keys: if i % 2 == 0 { vec![] } else { vec!["".into()] },
            typed_keys: vec![
                Key::Text("text".into()),
                Key::Tuple(
                    vec![
                        Scalar::Missing,
                        Scalar::Int(-(i as i64)),
                        Scalar::Float(-0.0),
                        Scalar::Bool(i % 2 == 0),
                        Scalar::Text("é".into()),
                    ]
                    .into(),
                ),
            ],
            values: vec![-0.0, 1.25, f64::from_bits(0x7ff8_0000_0000_00ab)],
            counts: vec![1, 2, 3],
            minimum_int: (i % 2 == 0).then_some(-(i as i64)),
            minimum_float: (i % 2 == 1).then_some(-0.0),
            maximum: Some(i as i64),
            auxiliary_int: None,
        })
        .collect();
    let mut expected = Sha256::new();
    for value in &rows {
        digest(value, &mut expected);
    }
    let mut frame = frame(&rows)?;
    ParquetWriter::new(File::create(&path)?)
        .with_row_group_size(Some(2))
        .set_parallel(false)
        .with_key_value_metadata(Some(KeyValueMetadata::from_static(vec![(
            "priva.schema".into(),
            SCHEMA.into(),
        )])))
        .finish(&mut frame)?;
    let groups = ParquetReader::new(File::open(&path)?)
        .get_metadata()?
        .row_groups
        .len();
    assert!(groups >= 3);
    let mut actual = Sha256::new();
    let mut seen = 0;
    scan(&path, |value| {
        assert_eq!(value.gene_order, seen);
        digest(&value, &mut actual);
        seen += 1;
        Ok(())
    })?;
    assert_eq!(seen, rows.len() as u64);
    assert_eq!(expected.finish_hex(), actual.finish_hex());
    std::fs::remove_dir_all(folder)?;
    Ok(())
}
