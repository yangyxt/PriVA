#!/usr/bin/env python3
"""Exact InterPro mapping producer comparison against the unchanged Python.

python interpro_builder_compare.py compare --binary PRIVA --repo FROZEN_SOURCE \
    --work-root EMPTY_DIRECTORY [--xml PINNED_PUBLISHER_XML_GZ]

The optional XML is an explicitly supplied, already acquired input. This test
never downloads resources or changes an installed cache. All successful cases
compare every native value with CPython's independently written trusted pickle.
Failures retain exit status, exception category and output preservation. XML
parser diagnostics are recorded, not required to mimic a Python traceback.
"""
from __future__ import annotations
import argparse
import gzip
import interpro_encoding_tables
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import cache_values_compare as values

BASE='''<interprodb><release/><interpro id="IPR1" type="Domain" short_name="short">
<name> Full name </name><class_list><classification class_type="GO" id="GO:1">
<category> process </category><description> binding </description></classification></class_list>
<abstract>first<b>bold</b> tail <![CDATA[<literal>]]> end</abstract>
<member_list><db_xref db="Pfam" dbkey="PF1"/><db_xref db="Pfam" dbkey="PF1"/>
<db_xref db="SMART" dbkey="SM2"/></member_list></interpro></interprodb>'''
SENTINEL=b'preserve prior output\n'


def save(path, value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')


def measured(command, folder, label):
    started=time.monotonic()
    with (folder/f'{label}.stdout').open('x') as stdout,(folder/f'{label}.stderr').open('x') as stderr:
        proc=subprocess.Popen(command,cwd=folder,stdout=stdout,stderr=stderr,
                              env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
        pid,status,usage=os.wait4(proc.pid,0)
        proc.returncode=os.waitstatus_to_exitcode(status)
    return dict(command=command,returncode=proc.returncode,wall_s=time.monotonic()-started,peak_rss_kib=usage.ru_maxrss)


def fixtures():
    cases=[]
    def add(name, xml=BASE, **kw):cases.append(dict(name=name,xml=xml,**kw))
    add('full-fields-duplicate-members')
    add('missing-id',BASE.replace(' id="IPR1"',''))
    add('empty-id',BASE.replace('id="IPR1"','id=""'))
    add('empty-defaults',BASE.replace('type="Domain"','type=""').replace('short_name="short"','short_name=""'))
    add('no-default-attributes',BASE.replace(' type="Domain" short_name="short"',''))
    add('mixed-name-text',BASE.replace(' Full name ','head<b>child</b>tail'))
    add('empty-name',BASE.replace(' Full name ',''))
    add('name-whitespace',BASE.replace(' Full name ',' \u00a0\u2003name\u3000 '))
    add('duplicate-name',BASE.replace('</name>','</name><name>ignored</name>',1))
    add('duplicate-class-list',BASE.replace('</class_list>','</class_list><class_list><classification class_type="GO" id="GO:ignored"/></class_list>',1))
    add('duplicate-member-list',BASE.replace('</member_list>','</member_list><member_list><db_xref db="OTHER" dbkey="ignored"/></member_list>',1))
    add('nested-member-not-direct',BASE.replace('<db_xref db="SMART" dbkey="SM2"/>','<wrapper><db_xref db="SMART" dbkey="SM2"/></wrapper>'))
    add('multiple-interpro',BASE.replace('</interprodb>',BASE[BASE.index('<interpro id='):BASE.index('</interprodb>')]+'</interprodb>'))
    add('nested-interpro-cleared',BASE.replace('<abstract>','<abstract><interpro id="IPRinside"><member_list><db_xref db="A" dbkey="B"/></member_list><name>inner</name></interpro>after'))
    add('nested-first-interpro-empty',BASE.replace('<release/>','<interpro/>'))
    for pad in (0,16000,16200,16300,16370,32760):
        nested='<interpro id="OUT"><abstract>before<interpro id="IN"><member_list><db_xref db="D" dbkey="IN"/></member_list></interpro>after</abstract><member_list><db_xref db="D" dbkey="OUT"/></member_list></interpro>'
        add(f'nested-feed-boundary-{pad}','<root>'+' '*pad+nested+'</root>')
    add('interpro-root',BASE.removeprefix('<interprodb><release/>').removesuffix('</interprodb>'))
    add('empty-interpro-first-end','<interpro/>')
    add('all-no-members','<root><interpro id="I"><name>X</name></interpro></root>')
    add('default-namespace',BASE.replace('<interprodb>','<interprodb xmlns="urn:test">'))
    add('prefixed-namespace',BASE.replace('<interprodb>','<interprodb xmlns:p="urn:test">').replace('<interpro id','<p:interpro id').replace('</interpro>','</p:interpro>'))
    add('default-namespace-reset',BASE.replace('<interprodb>','<interprodb xmlns="urn:test">').replace('<interpro id','<interpro xmlns="" id'))
    add('namespace-attribute-not-id',BASE.replace('<interprodb>','<interprodb xmlns:p="urn:test">').replace(' id="IPR1"',' p:id="IPR1"'))
    add('empty-go-id',BASE.replace('id="GO:1"','id=""'))
    add('lowercase-go-ignored',BASE.replace('class_type="GO"','class_type="go"'))
    add('missing-go-fields',BASE.replace('<category> process </category><description> binding </description>',''))
    add('go-first-text',BASE.replace(' binding ','before<i>child</i>after'))
    add('duplicate-go',BASE.replace('</classification>','</classification><classification class_type="GO" id="GO:1"><category> process </category><description> binding </description></classification>'))
    add('empty-member-db',BASE.replace('db="Pfam"','db=""'))
    add('unicode-db-upper',BASE.replace('Pfam','stra\u00dfe').replace('PF1','\u03b1\u65e5'))
    add('member-whitespace-preserved',BASE.replace('db="Pfam"','db=" Pfam "').replace('dbkey="PF1"','dbkey=" PF1 "'))
    add('xml-escaped',BASE.replace(' Full name ','&lt;&gt;&quot;&apos;&amp; &#65; &#x3b1;'))
    add('comments-pi-concatenate',BASE.replace(' Full name ','a<!--x-->b<?pi value?>c<![CDATA[d]]>e'))
    for sep,name in [('\r\n','crlf'),('\r','cr')]:add(name,BASE.replace('\n',sep))
    add('attribute-linebreak-normalization',BASE.replace('dbkey="PF1"','dbkey="A\tB\r\nC"'))
    add('attribute-escaped-linebreak',BASE.replace('dbkey="PF1"','dbkey="A&#9;B&#10;C"'))
    add('external-dtd-not-loaded','<!DOCTYPE interprodb SYSTEM "file:///not-loaded.dtd">'+BASE)
    add('internal-text-entity','<!DOCTYPE interprodb [<!ENTITY word "word &amp; more">]>'+BASE.replace(' Full name ','&word;'))
    add('internal-markup-entity','<!DOCTYPE interprodb [<!ENTITY word "one<b>two</b>three">]>'+BASE.replace('<abstract>','<abstract>&word;'))
    add('dtd-default-attributes','<!DOCTYPE interprodb [<!ATTLIST interpro type CDATA "DefaultType">]>'+BASE.replace(' type="Domain"',''))
    add('external-entity-refused','<!DOCTYPE interprodb [<!ENTITY secret SYSTEM "file:///not-loaded.txt">]>'+BASE.replace(' Full name ','&secret;'))
    for enc,label in [('utf-8','utf8-declaration'),('utf-16','utf16-bom'),('utf-16-le','utf16le'),('utf-16-be','utf16be'),('iso-8859-1','latin1'),('windows-1252','cp1252')]:
        text='<?xml version="1.0" encoding="'+enc+'"?>'+BASE.replace(' Full name ','caf\u00e9')
        add(label,text.encode(enc))
    add('utf8-bom',b'\xef\xbb\xbf'+BASE.encode())
    # The Expat fallback is CPython's whole 0..255 byte decode, not WHATWG
    # label equivalence. Exercise every pinned standard-library codec label
    # through the actual XML procedure, including parse and codec failures.
    for label in interpro_encoding_tables.catalogue():
        body=BASE.encode().replace(b' Full name ',b'probe \xa3\xe9\x81\xfe')
        add('codec_'+label,b'<?xml version="1.0" encoding="'+label.encode()+b'"?>'+body)
    for label in ('latin-1','LATIN_1','UTF-16-LE','Utf_16_Be','cp-1252','windows_1252'):
        body=BASE.encode().replace(b' Full name ',b'probe \xa3\xe9\x81\xfe')
        add('alias_'+label,b'<?xml version="1.0" encoding="'+label.encode()+b'"?>'+body)
    add('empty-document','')
    add('whitespace-document',' \n\t')
    add('truncated-xml',BASE[:-13])
    add('mismatched-tag',BASE.replace('</abstract>','</bad>'))
    add('duplicate-attribute',BASE.replace('id="IPR1"','id="IPR1" id="IPR2"'))
    add('unbound-prefix',BASE.replace('<name>','<p:name>').replace('</name>','</p:name>'))
    add('unknown-entity',BASE.replace(' Full name ','&missing;'))
    add('invalid-character',BASE.replace(' Full name ','&#0;'))
    add('multiple-roots',BASE+BASE)
    add('bad-utf8',BASE.encode().replace(b'Full',b'\xffFull'))
    add('missing-file',missing=True)
    add('directory-input',directory=True)
    add('input-not-gzip',raw=BASE.encode())
    add('empty-gzip-file',raw=b'')
    add('truncated-gzip',raw=gzip.compress(BASE.encode(),mtime=0)[:-4])
    damaged=bytearray(gzip.compress(BASE.encode(),mtime=0));damaged[-8]^=1
    add('bad-gzip-crc',raw=bytes(damaged))
    add('gzip-zero-padding',raw=gzip.compress(BASE.encode(),mtime=0)+b'\0'*17)
    add('gzip-trailing-junk',raw=gzip.compress(BASE.encode(),mtime=0)+b'bad')
    split=BASE.encode();add('multiple-gzip-members',raw=gzip.compress(split[:100],mtime=0)+gzip.compress(split[100:],mtime=0))
    add('prior-output-mode',mode=0o640)
    add('missing-parent',missing_parent=True)
    add('directory-output',output_directory=True)
    rng=random.Random(831301)
    for n in range(12):
        parts=['<root><release/>']
        for i in range(rng.randrange(1,100)):
            members=''.join(f'<db_xref db="db{rng.randrange(3)}" dbkey="KEY{rng.randrange(15)}"/>' for _ in range(rng.randrange(5)))
            parts.append(f'<interpro id="I{i}"><name>N{i}</name><member_list>{members}</member_list><abstract>A<b>B</b>C</abstract></interpro>')
        add(f'seeded-{n:02d}',''.join(parts)+'</root>')
    return cases


def reference(args):
    spec=importlib.util.spec_from_file_location('interpro_reference',args.source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    try:
        result=module.parse_interpro_xml_for_mapping(str(args.xml))
        module._write_interpro_mapping_atomic(result,str(args.output))
    except Exception as error:
        save(args.result,dict(ok=False,category=type(error).__name__,message=str(error)))
        raise SystemExit(1)
    save(args.result,dict(ok=True,keys=len(result),records=sum(map(len,result.values()))))


def compare(args):
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    sources={name:values.sha256(args.repo/name) for name in ('scripts/protein_domain_mapping.py','rust/tests/interpro_builder_compare.py','rust/tests/cache_values_compare.py')}
    initial_binary=values.sha256(args.binary)
    if (args.repo/'rust/src/interpro_mapping_builder/encodings.rs').read_text()!=interpro_encoding_tables.render():
        raise RuntimeError('pinned CPython codec tables differ from the reference environment')
    cases=fixtures()
    if args.xml:cases.append(dict(name='pinned-publisher',supplied=args.xml))
    results=[];failed=[]
    for case in cases:
        directory=root/case['name'];directory.mkdir()
        xml=directory/'input.xml.gz'
        if case.get('supplied'):shutil.copy2(case['supplied'],xml)
        elif case.get('directory'):xml.mkdir()
        elif not case.get('missing'):
            data=case['xml'].encode() if isinstance(case['xml'],str) else case['xml']
            xml.write_bytes(case.get('raw',gzip.compress(data,mtime=0)))
        output_parent=directory/'missing' if case.get('missing_parent') else directory
        pyout=output_parent/'python.pkl';rsout=output_parent/'native.rkyv'
        for path in (pyout,rsout):
            if case.get('missing_parent'):continue
            if case.get('output_directory'):path.mkdir()
            else:path.write_bytes(SENTINEL);path.chmod(case.get('mode',0o600))
        identity=values.sha256(xml) if xml.is_file() else None
        py=measured([sys.executable,str(Path(__file__).resolve()),'reference','--source',str(args.repo/'scripts/protein_domain_mapping.py'),'--xml',str(xml),'--output',str(pyout),'--result',str(directory/'reference.json')],directory,'python')
        native=measured([str(args.binary),'interpro-mapping','build','--xml',str(xml),'--output',str(rsout)],directory,'native')
        oracle=json.loads((directory/'reference.json').read_text())
        record=dict(name=case['name'],python=py,native=native,reference=oracle,match=False)
        if py['returncode']==native['returncode']==0:
            result=values.compare(SimpleNamespace(binary=args.binary,cache='interpro_mapping_pickle',pickle=pyout,archive=rsout,report=directory/'values.json',mismatches=directory/'mismatches.jsonl'))
            record.update(comparison=result,match=result['passed'] and (pyout.stat().st_mode&0o777)==(rsout.stat().st_mode&0o777))
        elif py['returncode']==native['returncode']==1:
            errors=(directory/'native.stderr').read_text()
            category=oracle['category']
            preserved=all((path.is_file() and path.read_bytes()==SENTINEL) if not case.get('missing_parent') and not case.get('output_directory') else not path.exists() if case.get('missing_parent') else path.is_dir() for path in (pyout,rsout))
            record.update(native_category_match=bool(re.search(r'\b'+re.escape(category)+r':',errors)),preserved=preserved)
            record['match']=record['native_category_match'] and preserved
        record['input_unchanged']=identity is None or identity==values.sha256(xml)
        record['match'] &= record['input_unchanged']
        save(directory/'case.json',record);results.append(record)
        if not record['match']:failed.append(case['name'])
    stable=sources=={name:values.sha256(args.repo/name) for name in sources} and initial_binary==values.sha256(args.binary)
    result=dict(passed=not failed and stable,cases=len(results),failed=failed,source_stable=stable,sources=sources,binary_sha256=initial_binary,results=results,
                scope='Full XML-mapping procedure; original VCF-format report workflow is separate')
    save(root/'comparison.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('sources','results')}),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('reference')
    for key in ('source','xml','output','result'):p.add_argument('--'+key,type=Path,required=True)
    p.set_defaults(function=reference)
    p=sub.add_parser('compare')
    for key in ('binary','repo','work-root'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--xml',type=Path);p.set_defaults(function=compare)
    args=parser.parse_args();args.function(args)

if __name__=='__main__':main()
