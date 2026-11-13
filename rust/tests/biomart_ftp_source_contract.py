#!/usr/bin/env python3
"""Check native FTP constants and mature dependency source pins.

check --source-root CHECKOUT [--baseline-root QUALIFIED_CHECKOUT]
This is a source contract check; it performs no scientific computation.
"""
import argparse,ast,hashlib,json,re,sys
from pathlib import Path
sys.dont_write_bytecode=True

def function(tree,name):return next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
def rust_method(text,name):
 match=re.search(r'(?:pub\(crate\) )?fn '+re.escape(name)+r'\(',text);assert match,name
 start=text.index('{',match.start());depth=1;at=start+1
 while depth:
  if text[at]=='{':depth+=1
  if text[at]=='}':depth-=1
  at+=1
 return text[match.start():at]
def check(a):
 root=a.source_root;source=(root/'scripts/biomart_ftp_nm_to_enst.py').read_text();tree=ast.parse(source)
 # All lookup SQL constants are copied verbatim, including whitespace. The
 # dynamic orphan query has the original alias/key interpolation in Rust.
 expected={};audit=function(tree,'relationship_audit')
 for node in audit.body:
  if isinstance(node,ast.Assign) and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='scalar' and isinstance(node.value.args[1],ast.Constant):expected[ast.literal_eval(node.targets[0].slice)]=node.value.args[1].value
 expected['invalid_transcript_rows']=next(n.value.args[0].left.args[1].value for n in audit.body if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.args[0],ast.Compare) and isinstance(n.value.args[0].left,ast.Call))
 assert json.loads((root/'rust/src/biomart/ftp/queries.json').read_text())==expected
 constants={}
 for node in tree.body:
  if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ['DATASET','PROJECTIONS','QUERY']:
   exec(compile(ast.Module(body=[node],type_ignores=[]),'<source-constants>','exec'),{'__builtins__':{}},constants)
 projections=json.loads((root/'rust/src/biomart/ftp/projections.json').read_text());assert list(projections)==list(constants['PROJECTIONS'])
 assert projections==json.loads(json.dumps(constants['PROJECTIONS']))
 facts=json.loads((root/'rust/src/biomart/facts.json').read_text());assert facts['ftp_specification']['query_sql']==constants['QUERY']
 vendor=root/'rust/third_party/sqlite3';record=json.loads((vendor/'provenance.json').read_text())
 for name,item in record['files'].items():
  payload=(vendor/name).read_bytes();assert len(payload)==item['bytes'] and hashlib.sha256(payload).hexdigest()==item['sha256']
 assert hashlib.sha3_256((vendor/'sqlite3.c').read_bytes()).hexdigest()==record['sqlite3_c_sha3_256']
 methods={}
 if a.baseline_root:
  before=(a.baseline_root/'rust/src/python_lines.rs').read_text();after=(root/'rust/src/python_lines.rs').read_text()
  for name in ['new','fill','next','readline','fill_for_readline','next_csv']:
   old=rust_method(before,name);new=rust_method(after,name)
   if name=='new':new=new.replace('            latin1: false,\n','')
   if name=='fill':
    extra='''        if self.latin1 {
            self.text = self.pending.iter().map(|byte| char::from(*byte)).collect();
            self.pending.clear();
            self.at = 0;
            return Ok(!self.eof || !self.text.is_empty());
        }
'''
    assert extra in new;new=new.replace(extra,'',1)
   assert old==new,name
   methods[name]=hashlib.sha256(old.encode()).hexdigest()
 print(json.dumps(dict(sql_constants_exact=True,projection_constants_order_exact=True,official_sqlite_source_pin_exact=True,existing_utf8_method_bodies=methods)))
def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(required=True);q=s.add_parser('check');q.add_argument('--source-root',type=Path,required=True);q.add_argument('--baseline-root',type=Path);q.set_defaults(run=check);a=p.parse_args();a.run(a)
if __name__=='__main__':main()
