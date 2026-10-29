#!/usr/bin/env python
# -*- coding: ASCII -*-
# Adapted from CADD-scripts v1.7.2 src/scripts/extract_scored.py.
# Copyright (c) University of Washington, Hudson-Alpha Institute for
# Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved.
# Permission is hereby granted, to all non-commercial users and licensees of
# CADD, to use, copy, modify, merge, and distribute this software subject to
# retaining the copyright and permission notices. See CADD-LICENSE for the
# complete upstream notice and disclaimer.
#
# PriVA change: upstream asks the tabix index once per variant, and with the
# published .tbi index every query unpacks the prescored table from the start
# of its 16 kb index window. This version reads the input in batches, merges
# positions on one chromosome that lie within SPAN_GAP of each other, and reads
# each merged region once. When PriVA's position index (TABLE.posidx, built by
# scripts/cadd_prescore_index.py during main_install) is present, a region is
# read from the last indexed row before it, a few BGZF blocks away. It writes
# the same rows in the same order: a row is kept for a variant only when its
# Pos, Ref and Alt strings equal the variant's, exactly as upstream compares
# them. Variants that upstream cannot query (unknown chromosome, position
# outside the index) still take the upstream code path.

import sys
import os
import pysam
import struct
import zlib
from array import array
from bisect import bisect_left
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="Path to scored variants.")
parser.add_option("-i", "--input", dest="input", help="Read variants from vcf file (default stdin)", default=None)
parser.add_option("--found_out", dest="found_out", help="Write found variants to file (default: stdout)", default=None)
parser.add_option("--header", dest="header", help="Write full header to output (default none)",
                  default=False, action="store_true")
(options, args) = parser.parse_args()

if options.input:
    stdin = open(options.input, 'r')
else:
    stdin = sys.stdin

if options.found_out:
    found_out = open(options.found_out, 'w')
else:
    found_out = sys.stdout

fpos, fref, falt = 1, 2, 3
if os.path.exists(options.path) and os.path.exists(options.path+".tbi"):
    filename = options.path
    sys.stderr.write("Opening %s...\n" % (filename))
    regionTabix = pysam.Tabixfile(filename, 'r')
    header = list(regionTabix.header)
    for line in header:
        if options.header:
            found_out.write(line+"\n")
        try:
            fref = line.split('\t').index('Ref')
            falt = line.split('\t').index('Alt')
        except ValueError:
            pass
else:
    raise IOError("No valid file with pre-scored variants.\n")

POSIDX_MAGIC = b"PRIVA-POSIDX\x01\x00\x00\x00"


class PositionIndex(object):
    """Checkpoints (position, virtual offset) per chromosome from TABLE.posidx."""

    def __init__(self, table, path):
        offset_type = 'L' if array('L').itemsize == 8 else 'Q'
        if sys.byteorder != 'little' or array('I').itemsize != 4:
            raise ValueError('unsupported platform for the position index')
        self.chromosomes = {}
        with open(path, 'rb') as stream:
            if stream.read(len(POSIDX_MAGIC)) != POSIDX_MAGIC:
                raise ValueError('not a PriVA position index')
            size, every, count = struct.unpack('<QII', stream.read(16))
            if size != os.path.getsize(table):
                raise ValueError('the index belongs to a different table')
            for _ in range(count):
                name = stream.read(struct.unpack('<I', stream.read(4))[0])
                if not isinstance(name, str):
                    name = name.decode('ascii')
                rows = struct.unpack('<Q', stream.read(8))[0]
                positions, offsets = array('I'), array(offset_type)
                positions.fromfile(stream, rows)
                offsets.fromfile(stream, rows)
                self.chromosomes[name] = (positions, offsets)
        self.data = open(table, 'rb')

    def start(self, chrom, first):
        """Virtual offset of the last indexed row before first, or None."""
        entry = self.chromosomes.get(chrom)
        if entry is None:
            return None
        index = bisect_left(entry[0], first) - 1
        return entry[1][index] if index >= 0 else None

    def rows(self, start, chrom, first, last):
        """Rows of chrom with first <= Pos <= last, in file order, from start on."""
        self.data.seek(start >> 16)
        skip = start & 0xFFFF
        pending = ''
        while True:
            header = self.data.read(18)
            if len(header) < 18:
                # The end of the file: a last row without a newline still counts.
                lines, pending = [pending] if pending else [], ''
            else:
                size = struct.unpack('<H', header[16:18])[0] + 1
                text = zlib.decompress(self.data.read(size - 18)[:-8], -15)
                if not isinstance(text, str):
                    text = text.decode('ascii')
                if skip:
                    text, skip = text[skip:], 0
                lines = (pending + text).split('\n')
                pending = lines.pop()
            for line in lines:
                tab = line.find('\t')
                if line[:tab] != chrom:
                    return
                position = int(line[tab+1:line.find('\t', tab+1)])
                if position > last:
                    return
                if position >= first:
                    yield line
            if len(header) < 18:
                return


position_index = None
if os.path.exists(filename + ".posidx"):
    try:
        position_index = PositionIndex(filename, filename + ".posidx")
    except (IOError, OSError, ValueError, struct.error) as error:
        sys.stderr.write("Ignoring %s.posidx: %s\n" % (filename, error))

# Positions closer than this many bases share one streamed region. A new
# query unpacks up to 16 kb of rows with the .tbi index but only a few BGZF
# blocks with the position index, so the position index merges only close
# positions.
SPAN_GAP = 256 if position_index else 4096
# Input lines held in memory before their results are written.
BATCH_SIZE = 20000
# Tabix coordinates are below 2^29; larger positions take the upstream path.
MAX_POSITION = 1 << 29
contigs = set(regionTabix.contigs)


def upstream_lookup(line):
    """The unchanged upstream handling of one input variant line."""
    try:
        fields = line.split('\t')
        found = False
        chrom = fields[0]
        pos = int(fields[1])
        lref, allele = fields[-2], fields[-1]
        for regionHit in regionTabix.fetch(chrom, pos-1, pos):
            vfields = regionHit.rstrip().split('\t')
            if (vfields[fref] == lref) and (vfields[falt] == allele) and (vfields[fpos] == fields[1]):
                found_out.write(regionHit+"\n")
                found = True

        if not found:
            sys.stdout.write(line + '\n')

    except ValueError:
        sys.stderr.write('Encountered uncovered chromosome\n')
        sys.stdout.write(line + '\n')


def spans(positions):
    """Merge sorted positions into (first, last) regions."""
    first = last = positions[0]
    for pos in positions[1:]:
        if pos - last > SPAN_GAP:
            yield first, last
            first = pos
        last = pos
    yield first, last


def region_rows(chrom, first, last):
    """Rows with first <= Pos <= last, as regionTabix.fetch(chrom, first-1, last) returns them."""
    if position_index is not None:
        start = position_index.start(chrom, first)
        if start is not None:
            return position_index.rows(start, chrom, first, last)
    return regionTabix.fetch(chrom, first-1, last)


def merged_lookup(batch):
    """Return {batch index: matching rows in file order} for merged variants."""
    by_chrom = {}
    for index, entry in enumerate(batch):
        if entry[0] == 'merged':
            by_chrom.setdefault(entry[2], []).append(index)
    found = {}
    for chrom, indices in by_chrom.items():
        # Key by the variant's position string: upstream keeps a row only when
        # the row's Pos string equals it, which also fixes the row's position.
        wanted = {}
        for index in indices:
            wanted.setdefault(batch[index][3], []).append(index)
        positions = sorted(set(batch[index][4] for index in indices))
        for first, last in spans(positions):
            for regionHit in region_rows(chrom, first, last):
                start = regionHit.find('\t')
                end = regionHit.find('\t', start+1)
                key = regionHit[start+1:end] if end >= 0 else regionHit[start+1:].rstrip()
                candidates = wanted.get(key)
                if not candidates:
                    continue
                vfields = regionHit.rstrip().split('\t')
                for index in candidates:
                    entry = batch[index]
                    if (vfields[fref] == entry[5]) and (vfields[falt] == entry[6]) and (vfields[fpos] == entry[3]):
                        found.setdefault(index, []).append(regionHit)
    return found


def flush(batch):
    found = merged_lookup(batch)
    for index, entry in enumerate(batch):
        kind, line = entry[0], entry[1]
        if kind == 'header':
            sys.stdout.write(line + '\n')
        elif kind == 'upstream':
            upstream_lookup(line)
        else:
            rows = found.get(index)
            if rows:
                for regionHit in rows:
                    found_out.write(regionHit+"\n")
            else:
                sys.stdout.write(line + '\n')


batch = []
for line in stdin:
    line = line.rstrip('\n\r')
    if line.startswith('#'):
        batch.append(('header', line))
    else:
        try:
            fields = line.split('\t')
            pos = int(fields[1])
        except ValueError:
            batch.append(('upstream', line))
        except Exception:
            # Upstream stops here after writing every earlier result.
            flush(batch)
            raise
        else:
            if fields[0] in contigs and 1 <= pos < MAX_POSITION:
                batch.append(('merged', line, fields[0], fields[1], pos, fields[-2], fields[-1]))
            else:
                batch.append(('upstream', line))
    if len(batch) >= BATCH_SIZE:
        flush(batch)
        batch = []
flush(batch)

if options.input:
    stdin.close()

if options.found_out:
    found_out.close()
