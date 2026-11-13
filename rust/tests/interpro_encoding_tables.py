#!/usr/bin/env python3
"""Freeze CPython's Expat fallback byte maps; no production Python dependency.

generate --output ../src/interpro_mapping_builder/encodings.rs
check --output ../src/interpro_mapping_builder/encodings.rs

CPython's PyUnknownEncodingHandler decodes bytes 0..255 together with replacement,
then requires exactly 256 Unicode characters. This includes its failure behavior;
WHATWG label aliases and per-byte decoding are not equivalent. Standard-library
codec tables are pinned here as small portable runtime data. Custom codecs
registered by an application are outside the original standalone script profile.
Reference: CPython 3.12 Modules/pyexpat.c, PyUnknownEncodingHandler (PSF license).
"""
from __future__ import annotations
import argparse
import encodings
import encodings.aliases
import json
import pkgutil
from pathlib import Path
import sys
import warnings


def catalogue():
    names=set(encodings.aliases.aliases)|set(encodings.aliases.aliases.values())
    names.update(item.name for item in pkgutil.iter_modules(encodings.__path__)
                 if not item.name.startswith('_') and item.name != 'aliases')
    result={}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for name in sorted(names):
            try:
                decoded=bytes(range(256)).decode(name,errors='replace')
                if len(decoded)!=256:
                    value=('error','ValueError','multi-byte encodings are not supported')
                else:
                    value=('map',tuple(-1 if char=='\ufffd' else ord(char) for char in decoded))
            except Exception as error:
                value=('error',type(error).__name__,str(error))
            result[name]=value
    return result


def render():
    if sys.version_info[:2]!=(3,12):raise RuntimeError('regenerate with the pinned CPython 3.12 profile')
    values=catalogue();maps={};lines=[
        '//! CPython 3.12 standard-library Expat fallback maps, including codec errors.',
        '//! Reproduce with rust/tests/interpro_encoding_tables.py; CPython (PSF license).',
        '//! No Python interpreter is used by the native XML builder.',
        '',
        'pub(super) enum Codec {',
        "    Map(&'static [i32; 256]),",
        "    Error(&'static str, &'static str),",
        '}',
        '',
        '#[rustfmt::skip]',
        'pub(super) fn lookup(name: &str) -> Option<Codec> {',
        '    match name {',
    ]
    for name,value in values.items():
        if value[0]=='map':
            index=maps.setdefault(value[1],len(maps))
            answer=f'Codec::Map(&MAP{index})'
        else:
            answer='Codec::Error('+json.dumps(value[1])+', '+json.dumps(value[2],ensure_ascii=False)+')'
        lines.append('        '+json.dumps(name)+' => Some('+answer+'),')
    lines += ['        _ => None,','    }','}','']
    for values,index in maps.items():
        lines += ['#[rustfmt::skip]',f'const MAP{index}: [i32; 256] = [']
        lines += ['    '+', '.join(map(str,values[start:start+16]))+',' for start in range(0,256,16)]
        lines += ['];','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    for name in ('generate','check'):
        command=sub.add_parser(name);command.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();expected=render()
    if args.command=='generate':args.output.write_text(expected)
    elif args.output.read_text()!=expected:raise SystemExit('CPython codec profile differs from the committed maps')
    print(json.dumps(dict(command=args.command,codec_labels=len(catalogue()),bytes=len(expected.encode()),passed=True)))


if __name__=='__main__':main()
