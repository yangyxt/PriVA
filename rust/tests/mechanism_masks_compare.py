#!/usr/bin/env python3
"""Compare opt-in Rust typed transport to the retained actual Python oracle."""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path

def compare(args):
    args.output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads(args.oracle_manifest.read_text())
    def verify():
        for name,sha in manifest.items():
            assert hashlib.sha256((args.oracle/name).read_bytes()).hexdigest()==sha,name
    verify()
    oracle=json.loads((args.oracle/'oracle.json').read_text())
    results=[]
    for case in oracle['cases']:
        name=case['name'];source=args.oracle/(name+'.input.csv');out=args.output/(name+'.actual.csv')
        cmd=[str(args.binary),'mechanism-masks','--input',str(source),'--output',str(out)]
        for key,option in [('var_plausible_patho_mechs','--profile-dtype'),('variant_effect','--effect-dtype'),('variant_lof_score','--lof-dtype'),('variant_gof_score','--gof-dtype'),('variant_dn_score','--dn-dtype')]:
            dtype=case['dtypes'].get(key,'object')
            if dtype in ('object','string','Float64','Int64','boolean'):cmd.extend([option,dtype])
        before=hashlib.sha256(source.read_bytes()).hexdigest()
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
        (args.output/(name+'.stderr.log')).write_text(p.stderr)
        if 'exception' in case:
            error=case['exception'];expected=f"priva: {error['type']}: {error['message']}"
            ok=p.returncode==2 and p.stderr.strip()==expected and not out.exists()
        else:
            expected=list(csv.reader((args.oracle/(name+'.expected.csv')).open(newline='')))
            actual=list(csv.reader(out.open(newline=''))) if out.exists() else None
            ok=p.returncode==0 and actual==expected
        assert hashlib.sha256(source.read_bytes()).hexdigest()==before
        results.append(dict(name=name,rows=case['rows'],passed=ok,returncode=p.returncode))
    verify()
    result=dict(oracle_manifest_sha256=hashlib.sha256(args.oracle_manifest.read_bytes()).hexdigest(),completed=all(x['passed'] for x in results),cases=results,oracle_sha256=hashlib.sha256((args.oracle/'oracle.json').read_bytes()).hexdigest(),binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),scientific_acceptance=False,scope='Typed primitive-cell differential boundary; not full ranking or family parity')
    (args.output/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if not result['completed']:raise SystemExit(1)

def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True);c=s.add_parser('compare')
    for name in ('binary','oracle','oracle-manifest','output'):c.add_argument('--'+name,type=Path,required=True)
    compare(p.parse_args())
if __name__=='__main__':main()
