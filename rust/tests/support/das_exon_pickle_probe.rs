//! Standalone metadata projection; never imports or executes pickle callables.
#[path = "../../src/das_exon_pickle.rs"]
mod das_exon_pickle;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

fn hex(value: &str) -> String {
    if value.is_empty() {
        "-".to_owned()
    } else {
        value
            .as_bytes()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let options: HashMap<&str, &str> = args[usize::from(!args.is_empty())..]
        .chunks_exact(2)
        .map(|pair| (pair[0].as_str(), pair[1].as_str()))
        .collect();
    if args.first().map(String::as_str) != Some("extract")
        || args.len() != 5
        || options.len() != 2
        || !["--input", "--output"]
            .iter()
            .all(|key| options.contains_key(key))
    {
        return Err("usage: das-exon-pickle-probe extract --input PICKLE --output LINES".into());
    }
    let projection = das_exon_pickle::read(File::open(options["--input"])?)?;
    let mut output = BufWriter::new(File::create_new(options["--output"])?);
    for row in &projection.sets {
        write!(output, "{}", row.path.len())?;
        for key in &row.path {
            write!(output, " {}", hex(key))?;
        }
        write!(output, " {}", row.items.len())?;
        for item in &row.items {
            write!(output, " {}", hex(item))?;
        }
        writeln!(output)?;
    }
    output.flush()?;
    let s = projection.statistics;
    println!(
        "{{\"bytes\":{},\"skipped_buffer_bytes\":{},\"opcodes\":{},\"memo_entries\":{},\"retained_memo_entries\":{},\"dictionaries\":{},\"sets\":{},\"set_items\":{},\"maximum_stack\":{},\"numpy_reductions\":{},\"projected_sets\":{}}}",
        s.bytes,
        s.skipped_buffer_bytes,
        s.opcodes,
        s.memo_entries,
        s.retained_memo_entries,
        s.dictionaries,
        s.sets,
        s.set_items,
        s.maximum_stack,
        s.numpy_reductions,
        projection.sets.len()
    );
    Ok(())
}
