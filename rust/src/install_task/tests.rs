use super::content::{Input, Places, content, listing};
use super::digest::{
    called_commands, called_words, data_names, dynamic_command, script_names, toplevel_lines,
};
use super::plan::{add_months, glob_match, header_release, sorted_lines_sha256};
use super::{delete_key, set_key, yaml_value};
use crate::sha256::sha256_hex;

fn scratch(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("install_task_test_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn message_lines_count_only_command_positions() {
    let text = "f () \n{ \n    log \"run whole_installer then $(helper_one x) | helper_two\";\n    other_call arg_1 2nd\n}\n";
    let words: Vec<String> = called_words(text).into_iter().collect();
    assert!(words.contains(&"helper_one".to_owned()));
    assert!(words.contains(&"helper_two".to_owned()));
    assert!(words.contains(&"other_call".to_owned()));
    assert!(words.contains(&"arg_1".to_owned()));
    assert!(!words.contains(&"whole_installer".to_owned()));
    assert!(!words.contains(&"2nd".to_owned()));
}

#[test]
fn script_and_data_names_follow_grep() {
    let text = r#"python "${SCRIPT_DIR}/tool_a.py" x.sh.bak; $SCRIPT_DIR/sub/b.R run.bash "${DATA_DIR}/clingen/map.${assembly}.tsv" DATA_DIR/x"#;
    assert_eq!(
        script_names(text),
        vec!["tool_a.py", "x.sh", "sub/b.R", "run.bash"]
    );
    assert_eq!(data_names(text), vec!["clingen/map.${assembly}.tsv", "x"]);
}

#[test]
fn toplevel_keeps_lines_outside_functions() {
    let text =
        "A=1\n# comment\nfunction f() {\n    body\n}\nB=2\n\nfunction g() {\n    x\n}\ntrailing\n";
    assert_eq!(toplevel_lines(text), "A=1\nB=2\ntrailing\n");
}

#[test]
fn months_are_calendar_months() {
    assert_eq!(add_months("2026-01-31", 1).unwrap(), "2026-02-28");
    assert_eq!(add_months("2024-01-31", 1).unwrap(), "2024-02-29");
    assert_eq!(add_months("2026-11-17", 3).unwrap(), "2027-02-17");
}

#[test]
fn settings_patterns_match_like_bash() {
    assert!(glob_match("cadd_*md5*", "cadd_GRCh38_snv_anno_md5"));
    assert!(glob_match("clinvar_*", "clinvar_vcf"));
    assert!(!glob_match("clinvar_*", "clinvar"));
    assert!(glob_match("assembly", "assembly"));
}

#[test]
fn release_headers_come_from_the_final_response() {
    let headers = b"HTTP/1.1 302 Found\r\nETag: \"old\"\r\n\r\nHTTP/2 200\r\nLast-Modified: Tue, 01 Jan 2026\r\nETag: \"New\"\r\n\r\n";
    assert_eq!(
        header_release(headers),
        "etag: \"new\"|last-modified: tue, 01 jan 2026|"
    );
    assert_eq!(header_release(b""), "||");
}

#[test]
fn sorted_lines_hash_like_sort_sha256sum() {
    assert_eq!(sorted_lines_sha256(b"b\na\n"), sha256_hex(b"a\nb\n"));
    assert_eq!(sorted_lines_sha256(b"b\na"), sha256_hex(b"a\nb\n"));
    assert_eq!(sorted_lines_sha256(b""), sha256_hex(b""));
}

#[test]
fn configuration_keys_are_set_and_read_like_the_shell_runner() {
    let dir = scratch("config");
    let config = dir.join("c.yaml");
    std::fs::write(
        &config,
        "a: 1  # one\nflag: false\non: true\nlist: [x, 2]\nb : \"text\"\nlast: z\n",
    )
    .unwrap();
    assert_eq!(yaml_value(&config, "a").unwrap(), "1");
    assert_eq!(yaml_value(&config, "flag").unwrap(), "");
    assert_eq!(yaml_value(&config, "on").unwrap(), "true");
    assert_eq!(yaml_value(&config, "list").unwrap(), "[\"x\",2]");
    assert_eq!(yaml_value(&config, "missing").unwrap(), "");
    set_key(&config, "b", "/new/path").unwrap();
    set_key(&config, "added", "v").unwrap();
    assert_eq!(yaml_value(&config, "b").unwrap(), "/new/path");
    delete_key(&config, "flag").unwrap();
    let text = std::fs::read_to_string(&config).unwrap();
    assert_eq!(
        text,
        "a: 1  # one\non: true\nlist: [x, 2]\nb : /new/path\nlast: z\nadded: v\n"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn priva_commands_are_read_from_the_lines_that_run_them() {
    assert_eq!(
        called_commands(r#"  "${PRIVA_BIN}" mondo-scope --out x"#),
        vec!["mondo-scope"]
    );
    assert_eq!(
        called_commands(r#""$PRIVA_CACHE_INSTALL_BINARY" builder-identity --step A"#),
        vec!["builder-identity"]
    );
    assert!(dynamic_command(r#""${PRIVA_BIN}" "$command" --in x"#));
    assert!(!dynamic_command(r#""${PRIVA_BIN}" mondo-scope"#));
}

#[cfg(unix)]
#[test]
fn content_names_results_not_the_folder_they_sit_in() {
    use std::os::unix::fs::symlink;
    let dir = scratch("content");
    let upstream = dir.join("UP/lineage1");
    std::fs::create_dir_all(upstream.join("data")).unwrap();
    std::fs::write(upstream.join("data/ref.fa"), "ACGT\n").unwrap();
    let tree = dir.join("tree");
    std::fs::create_dir_all(tree.join("data")).unwrap();
    std::fs::write(tree.join("data/shipped.tsv"), "x\n").unwrap();
    let inputs = vec![Input {
        id: "UP".into(),
        folder: upstream.clone(),
        content: "abc".into(),
    }];
    let make = |lineage: &str, body: &str| {
        let own = dir.join("STEP").join(lineage);
        let data = own.join("data");
        std::fs::create_dir_all(data.join("sub")).unwrap();
        std::fs::write(data.join("sub/out.tsv"), body).unwrap();
        symlink(upstream.join("data/ref.fa"), data.join("ref.fa")).unwrap();
        symlink("sub/out.tsv", data.join("alias.tsv")).unwrap();
        let keys = vec![
            (
                "out".to_owned(),
                data.join("sub/out.tsv").to_string_lossy().into_owned(),
            ),
            (
                "shipped".to_owned(),
                tree.join("data/shipped.tsv").to_string_lossy().into_owned(),
            ),
        ];
        let places = Places {
            own: own.clone(),
            tree: tree.clone(),
            inputs: &inputs,
        };
        listing(&data, &keys, &places, 2).unwrap()
    };
    let first = make("aaaa", "1\t2\n");
    let moved = make("bbbb", "1\t2\n");
    let changed = make("cccc", "1\t3\n");
    assert_eq!(first, moved);
    assert_eq!(content(&first), content(&moved));
    assert_ne!(content(&first), content(&changed));
    let lines: Vec<&str> = first.lines().collect();
    assert_eq!(lines[0], "l alias.tsv own data/sub/out.tsv");
    assert_eq!(lines[1], "l ref.fa input UP data/ref.fa abc");
    assert_eq!(lines[2], "d sub");
    assert!(lines[3].starts_with("f sub/out.tsv "));
    assert_eq!(lines[4], "k out own data/sub/out.tsv");
    assert!(lines[5].starts_with("k shipped tree data/shipped.tsv "));
    let _ = std::fs::remove_dir_all(&dir);
}
