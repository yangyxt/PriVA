# SQLite source dependency

The unmodified SQLite 3.46.0 amalgamation and public headers are in the public domain. Publisher identities are in `provenance.json`; `compile_options.txt` records the selected engine features and limits.

`rust/build.rs` compiles this source into the optional BioMart FTP projection backend using the pinned `cc` build dependency. The original SQL and transaction order remain in the BioMart module. No runtime SQLite installation or Python adapter is used. The normal source manifest includes this directory for offline builds.

SQLite compile-option labels omit the values of some macros. The explicit zero values for auto-vacuum and recursive triggers preserve the original engine defaults; `provenance.json` also records the effective connection defaults.
