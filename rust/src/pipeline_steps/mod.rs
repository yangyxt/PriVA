//! The steps of the analysis pipeline without PriVA's shell scripts: each
//! runs the tools' own commands (bcftools, tabix, VEP, CADD.sh) and priva's
//! own, with the arguments the shell functions gave them. `priva task run`
//! calls them inside a Nextflow task; `priva task step` runs one by hand.
//!
//! One file per step, and the files they share (`tools`), so that a step's
//! hub lineage, which counts the Rust files the step's own file reaches,
//! moves only when that step's code does. Nothing here names this file's
//! children for the same reason: the dispatch lives in the task code.
pub(crate) mod annotation;
pub(crate) mod cadd;
pub(crate) mod clinvar;
pub(crate) mod control;
pub(crate) mod family;
pub(crate) mod gnomad;
pub(crate) mod preprocess;
pub(crate) mod tools;
pub(crate) mod vep;
