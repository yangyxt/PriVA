fn main() {
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
}
