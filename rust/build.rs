fn main() {
    // The FTP projection uses the same pinned mature SQLite engine as its
    // original Python environment. Vendor source avoids a runtime host-library
    // dependency and makes normal offline bootstrap builds self-contained.
    let sqlite_root = "third_party/sqlite3";
    let mut sqlite = cc::Build::new();
    sqlite
        .file(format!("{sqlite_root}/sqlite3.c"))
        .opt_level(2)
        .debug(false)
        .warnings(false);
    for line in std::fs::read_to_string(format!("{sqlite_root}/compile_options.txt"))
        .expect("pinned SQLite compile options")
        .lines()
    {
        let (name, value) = line.split_once('=').expect("SQLite definition=value");
        sqlite.define(name, value);
    }
    sqlite.compile("priva_sqlite3");
    for library in ["m", "dl", "pthread"] {
        println!("cargo:rustc-link-lib={library}");
    }
    println!("cargo:rerun-if-changed={sqlite_root}");

    // Match PyYAML's mature event parser/emitter. Package export bytes are
    // qualified against SafeDumper; no independent YAML formatting policy.
    let yaml = pkg_config::Config::new()
        .atleast_version("0.2.5")
        .statik(true)
        .cargo_metadata(false)
        .probe("yaml-0.1")
        .expect("building priva requires libyaml >=0.2.5 and its static archive");
    let yaml_dir = std::path::PathBuf::from(
        pkg_config::get_variable("yaml-0.1", "libdir").expect("libyaml libdir is required"),
    );
    assert!(
        yaml_dir.join("libyaml.a").is_file(),
        "libyaml.a is required"
    );
    cc::Build::new()
        .file("src/installed_environment/yaml_bridge.c")
        .includes(&yaml.include_paths)
        .warnings_into_errors(true)
        .compile("priva_yaml_bridge");
    for path in yaml.link_paths.iter().chain(std::iter::once(&yaml_dir)) {
        println!("cargo:rustc-link-search=native={}", path.display());
    }
    println!("cargo:rustc-link-lib=static=yaml");
    println!("cargo:rustc-env=PRIVA_LIBYAML_VERSION={}", yaml.version);
    println!("cargo:rerun-if-changed=src/installed_environment/yaml_bridge.c");
    // ElementTree uses Expat. Keep the same mature XML backend instead of
    // independently changing namespace, encoding, DTD or text semantics.
    // The static dependency follows the binary across execution environments.
    // PKG_CONFIG_PATH can select the build environment's lib/pkgconfig.
    let library = pkg_config::Config::new()
        .atleast_version("2.6.0")
        .statik(true)
        .cargo_metadata(false)
        .probe("expat")
        .expect(
            "building priva requires Expat >=2.6.0, its static library and pkg-config metadata",
        );
    // pkg-config's Rust helper deliberately links /usr libraries dynamically,
    // even with statik(true). Request the static Expat archive explicitly, in
    // both system and Conda layouts. It travels inside the qualified binary.
    let libdir = std::path::PathBuf::from(
        pkg_config::get_variable("expat", "libdir").expect("Expat libdir metadata is required"),
    );
    assert!(
        libdir.join("libexpat.a").is_file(),
        "Expat's static archive libexpat.a is required"
    );
    for path in library.link_paths.iter().chain(std::iter::once(&libdir)) {
        println!("cargo:rustc-link-search=native={}", path.display());
    }
    println!("cargo:rustc-link-lib=static=expat");
    // Keep platform support libraries such as libm at their ordinary linkage.
    for dependency in library.libs.iter().filter(|name| name.as_str() != "expat") {
        println!("cargo:rustc-link-lib={dependency}");
    }
    for arguments in library.ld_args {
        println!("cargo:rustc-link-arg=-Wl,{}", arguments.join(","));
    }
    println!("cargo:rustc-env=PRIVA_EXPAT_VERSION={}", library.version);
    // Python's gzip module compresses through zlib. Linking the same static
    // zlib lets a native cache writer produce Python's compressed bytes.
    let zlib = pkg_config::Config::new()
        .statik(true)
        .cargo_metadata(false)
        .probe("zlib")
        .expect("building priva requires zlib, its static library and pkg-config metadata");
    let zlib_dir = std::path::PathBuf::from(
        pkg_config::get_variable("zlib", "libdir").expect("zlib libdir metadata is required"),
    );
    assert!(
        zlib_dir.join("libz.a").is_file(),
        "zlib's static archive libz.a is required"
    );
    for path in zlib.link_paths.iter().chain(std::iter::once(&zlib_dir)) {
        println!("cargo:rustc-link-search=native={}", path.display());
    }
    println!("cargo:rustc-link-lib=static=z");
    println!("cargo:rustc-env=PRIVA_ZLIB_VERSION={}", zlib.version);
}
