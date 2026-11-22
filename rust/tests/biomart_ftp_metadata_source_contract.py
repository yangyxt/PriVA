#!/usr/bin/env python3
"""Check additive gzip diagnostics and unchanged original BioMart identities.

check --source-root SOURCE --baseline-root QUALIFIED_PREPARED_SOURCE
No producer, acquisition or scientific caller is run.
"""
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
def digest(value):return hashlib.sha256(value).hexdigest()
def check(a):
 name='rust/src/python_io.rs';before=(a.baseline_root/name).read_text().split('#[cfg(test)]')[0];after=(a.source_root/name).read_text().split('#[cfg(test)]')[0]
 constructor='''    /// Keep the mature decoder's reason for helpers that compare zlib's
    /// diagnostic text. Existing callers retain the established default.
    pub(crate) fn with_zlib_diagnostics(source: R) -> Self {
        let mut reader = Self::new(source);
        reader.detailed_errors = true;
        reader
    }

'''
 diagnostic='''                            if self.detailed_errors {
                                if let Some(reason) = cause.message() {
                                    return GzipFailure::Zlib(format!(
                                        "Error -3 while decompressing data: {reason}"
                                    ))
                                    .into_error();
                                }
                            }
'''
 projected=after
 for added in ['    detailed_errors: bool,\n','            detailed_errors: false,\n',constructor,diagnostic]:
  assert projected.count(added)==1;projected=projected.replace(added,'')
 assert projected.count('.map_err(|cause| {')==1;projected=projected.replace('.map_err(|cause| {','.map_err(|_| {')
 assert projected==before
 allowed={'rust/src/python_io.rs','rust/src/biomart/cli.rs','rust/src/biomart/ftp/mod.rs'};same={}
 for path in (a.baseline_root/'rust/src').rglob('*.rs'):
  relative=str(path.relative_to(a.baseline_root))
  if relative in allowed:continue
  assert path.read_bytes()==(a.source_root/relative).read_bytes(),relative;same[relative]=digest(path.read_bytes())
 for name in ['scripts/biomart_nm_to_enst.py','scripts/biomart_ftp_nm_to_enst.py','rust/Cargo.toml','rust/Cargo.lock','rust/build.rs']:
  assert (a.baseline_root/name).read_bytes()==(a.source_root/name).read_bytes(),name
 spec=importlib.util.spec_from_file_location('pinned_ftp_source',a.source_root/'scripts/biomart_ftp_nm_to_enst.py');source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
 facts=json.loads((a.source_root/'rust/src/biomart/facts.json').read_text());assert source.specification()==facts['ftp_specification']
 print(json.dumps(dict(passed=True,unchanged_gzip_default_body_sha256=digest(before.encode()),unchanged_existing_rust_files=same,original_assets_and_build_dependencies_unchanged=True,source_pins_exact=True),sort_keys=True))
def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(required=True);q=s.add_parser('check');q.add_argument('--source-root',type=Path,required=True);q.add_argument('--baseline-root',type=Path,required=True);q.set_defaults(run=check);a=p.parse_args();a.run(a)
if __name__=='__main__':main()
