#!/usr/bin/env python3
"""Reorder top-level CLI descriptor elements into Slicer canonical order.

Operates on raw text. ElementTree round-tripping destroys comments and
CDATA sections, and reformats the whole file; neither is acceptable here.
Only the relative order of whole element blocks changes -- every block's
bytes, and everything outside the header region, are preserved verbatim.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cli_xml_lib as lib


def _block_pattern(tag):
    return re.compile(
        rf"(?P<lead>[ \t]*)<{re.escape(tag)}(?:\s[^>]*)?>.*?</{re.escape(tag)}>[ \t]*\n?",
        re.DOTALL,
    )


def reorder_text(source):
    """Return source with header elements in canonical order."""
    open_match = re.search(r"<executable[^>]*>[ \t]*\n?", source)
    if not open_match:
        return source
    header_start = open_match.end()

    params_match = re.search(r"[ \t]*<parameters", source[header_start:])
    header_end = header_start + (
        params_match.start() if params_match else len(source) - header_start
    )

    header = source[header_start:header_end]

    found = []
    for tag in lib.CANONICAL_ORDER:
        match = _block_pattern(tag).search(header)
        if match:
            found.append((tag, match.group(0), match.start(), match.end()))

    if not found:
        return source

    found.sort(key=lambda item: item[2])
    source_order = [tag for tag, _, _, _ in found]
    canonical_order = sorted(source_order, key=lib.CANONICAL_ORDER.index)
    if source_order == canonical_order:
        return source

    # Each gap between two source-adjacent blocks travels with the block
    # that precedes it, so interleaved comments stay pinned to their
    # original neighbour even after the blocks are reordered.
    leading_gap = header[: found[0][2]]
    trailing_gap = header[found[-1][3] :]
    block_of = {tag: block for tag, block, _, _ in found}
    gap_after = {
        found[i][0]: header[found[i][3] : found[i + 1][2]]
        for i in range(len(found) - 1)
    }

    canonical_tags = sorted(source_order, key=lib.CANONICAL_ORDER.index)
    ordered = leading_gap + "".join(
        block_of[tag] + gap_after.get(tag, "") for tag in canonical_tags
    )
    return source[:header_start] + ordered + trailing_gap + source[header_end:]


def main(argv):
    paths = [Path(p) for p in argv[1:]] if len(argv) > 1 else lib.descriptor_paths()
    changed = 0
    for path in paths:
        original = path.read_text()
        updated = reorder_text(original)
        if updated != original:
            path.write_text(updated)
            changed += 1
            print(f"reordered {path}")
    print(f"{changed} file(s) reordered.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
