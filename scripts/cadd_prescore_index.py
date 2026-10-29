#!/usr/bin/env python3
"""Build PriVA's position index for a CADD prescored table.

cadd_prescore_index.py build TABLE [--threads N] [--every K]
cadd_prescore_index.py check TABLE

The published .tbi index has 16 kb windows, so CADD's lookup unpacks up to
16 kb of rows (tens of MB of text in the genome-wide SNV table) to reach one
position. TABLE.posidx lists, for every K-th BGZF block of the table, the
chromosome and position of the first row that starts in that block and that
row's virtual offset. PriVA's extract_scored.py (scripts/cadd_scripts) jumps
to the last listed row before a position and reads forward, decompressing
about K blocks. The index is built in one parallel pass over the table and
published atomically; it records the table's size, so an index never pairs
with a different table.

Layout (little-endian): MAGIC, table size (uint64), K (uint32), number of
chromosomes (uint32); per chromosome in file order: name length (uint32),
name, count (uint64), positions (count uint32), virtual offsets (count uint64).
"""
import argparse
from array import array
from multiprocessing import Pool
import os
from pathlib import Path
import struct
import sys
import tempfile
import zlib

MAGIC = b"PRIVA-POSIDX\x01\x00\x00\x00"
# Every BGZF block written by htslib starts with these 16 bytes, then BSIZE.
BLOCK_START = bytes.fromhex("1f8b08040000000000ff060042430200")


def block_size(header):
    if len(header) < 18 or header[:16] != BLOCK_START:
        raise ValueError("Not a BGZF block")
    return struct.unpack("<H", header[16:18])[0] + 1


def first_block(stream, start, end):
    """Return the offset of the first BGZF block starting at or after start."""
    if start == 0:
        return 0
    stream.seek(start)
    window = stream.read(1 << 20)
    found = window.find(BLOCK_START)
    while found >= 0 and start + found < end:
        offset = start + found
        stream.seek(offset)
        try:
            size = block_size(stream.read(18))
            stream.seek(offset + size)
            following = stream.read(18)
            # A real block is followed by another block or ends the file.
            if not following or following[:16] == BLOCK_START:
                return offset
        except ValueError:
            pass
        found = window.find(BLOCK_START, found + 1)
    return end


def scan(task):
    """Checkpoints for blocks in [start, end): [(chromosome, positions, virtual offsets)]."""
    path, start, end, every = task
    checkpoints = []
    with open(path, "rb", buffering=1 << 20) as stream:
        offset = first_block(stream, start, end)
        count = 0
        while offset < end:
            stream.seek(offset)
            header = stream.read(18)
            if not header:
                break
            size = block_size(header)
            if count % every == 0:
                text = zlib.decompressobj(-15).decompress(stream.read(size - 18))
                newline = text.find(b"\n")
                # The first row starting inside this block, if complete enough
                # to show its chromosome and position.
                if 0 <= newline < len(text) - 1:
                    fields = text[newline + 1:].split(b"\t", 2)
                    if len(fields) == 3 and not fields[0].startswith(b"#"):
                        chromosome = fields[0].decode("ascii")
                        if not checkpoints or checkpoints[-1][0] != chromosome:
                            checkpoints.append((chromosome, array("I"), array("Q")))
                        checkpoints[-1][1].append(int(fields[1]))
                        checkpoints[-1][2].append((offset << 16) | (newline + 1))
            count += 1
            offset += size
    return checkpoints


def build(table, threads=1, every=4):
    table = Path(table)
    size = table.stat().st_size
    parts = max(1, threads) * 4
    bounds = [size * part // parts for part in range(parts + 1)]
    tasks = [(str(table), bounds[part], bounds[part + 1], every) for part in range(parts)]
    with Pool(max(1, threads)) as pool:
        results = pool.map(scan, tasks)
    order, chromosomes = [], {}
    for checkpoints in results:
        for chromosome, part_positions, part_offsets in checkpoints:
            if chromosome not in chromosomes:
                order.append(chromosome)
                chromosomes[chromosome] = (array("I"), array("Q"))
            elif order[-1] != chromosome:
                raise ValueError(f"{table} is not grouped by chromosome at {chromosome}")
            positions, offsets = chromosomes[chromosome]
            previous = positions[-1] if positions else 0
            for position in part_positions:
                if position < previous:
                    raise ValueError(f"{table} is not sorted by position at {chromosome}:{position}")
                previous = position
            positions.extend(part_positions)
            offsets.extend(part_offsets)
    if sys.byteorder != "little":
        raise ValueError("The position index layout is little-endian")
    target = Path(str(table) + ".posidx")
    handle, temporary = tempfile.mkstemp(prefix="." + target.name + ".", dir=target.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(MAGIC + struct.pack("<QII", size, every, len(order)))
            for chromosome in order:
                positions, offsets = chromosomes[chromosome]
                name = chromosome.encode("ascii")
                stream.write(struct.pack("<I", len(name)) + name + struct.pack("<Q", len(positions)))
                positions.tofile(stream)
                offsets.tofile(stream)
        if table.stat().st_size != size:
            raise ValueError(f"{table} changed while it was indexed")
        os.chmod(temporary, 0o644)
        os.replace(temporary, target)
    finally:
        if os.path.lexists(temporary):
            os.unlink(temporary)
    return target, sum(len(chromosomes[name][0]) for name in order)


def check(table):
    """Return True when TABLE.posidx exists and belongs to TABLE."""
    index = Path(str(table) + ".posidx")
    try:
        with index.open("rb") as stream:
            header = stream.read(len(MAGIC) + 16)
    except OSError:
        return False
    return (len(header) == len(MAGIC) + 16 and header[:len(MAGIC)] == MAGIC
            and struct.unpack("<Q", header[len(MAGIC):len(MAGIC) + 8])[0] == Path(table).stat().st_size)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("procedure", choices=("build", "check"))
    parser.add_argument("table", type=Path, help="BGZF-compressed, tabix-sorted prescored table")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--every", type=int, default=4, help="Index every K-th BGZF block (default 4)")
    args = parser.parse_args()
    try:
        if args.procedure == "check":
            sys.exit(0 if check(args.table) else 1)
        target, count = build(args.table, args.threads, max(1, args.every))
        print(f"Wrote {target} with {count} positions")
    except (OSError, ValueError) as error:
        sys.exit("ERROR: " + str(error))


if __name__ == "__main__":
    main()
