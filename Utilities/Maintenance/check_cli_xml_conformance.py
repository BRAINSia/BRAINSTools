#!/usr/bin/env python3
"""Fail if any CLI descriptor violates Slicer Execution Model conventions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cli_xml_lib as lib


def main():
    violations = lib.check_all()
    for violation in violations:
        print(f"ERROR: {violation}", file=sys.stderr)
    if violations:
        print(f"\n{len(violations)} conformance violation(s).", file=sys.stderr)
        return 1
    print(f"OK: {len(lib.descriptor_paths())} descriptors conformant.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
