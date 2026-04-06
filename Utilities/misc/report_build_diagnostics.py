#!/usr/bin/env python3
"""Extract and print build warnings/errors from CTest Build.xml.

ctest_build() writes compiler diagnostics to Build.xml for CDash but
does not print them to stdout.  This script reads Build.xml from the
most recent test run, filters to the project source tree only (excluding
ExternalProject dependency build trees), and prints the <Text> content
of every <Warning> and <Error> element so diagnostics appear directly
in CI logs.

Exit codes:
  0  — no BRAINSTools diagnostics found
  1  — any BRAINSTools warnings or errors found (intended to fail CI)
  2  — usage error

Usage:
    python report_build_diagnostics.py <build-directory> [<source-directory>]

Arguments:
    build-directory   Root of the CMake binary tree (contains Testing/TAG).
    source-directory  Optional project source root.  Only warnings/errors
                      whose SourceFile starts with this path are reported.
                      Paths under <source-directory>/build/ are excluded so
                      ExternalProject dependency diagnostics are suppressed.
                      Defaults to the parent of <build-directory>.
"""

import os
import sys
import xml.etree.ElementTree as ET


def _is_project_file(source_file: str, source_dir: str, build_dir: str) -> bool:
    """Return True if source_file belongs to the project, not to an EP."""
    if not source_file:
        return False
    sf = os.path.realpath(source_file)
    sd = os.path.realpath(source_dir)
    bd = os.path.realpath(build_dir)
    # Must be under the source directory ...
    if not sf.startswith(sd + os.sep) and sf != sd:
        return False
    # ... but NOT under the build directory (ExternalProject trees live there).
    if sf.startswith(bd + os.sep) or sf == bd:
        return False
    return True


def main() -> int:
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(
            f"Usage: {sys.argv[0]} <build-directory> [<source-directory>]",
            file=sys.stderr,
        )
        return 2

    build_dir = os.path.abspath(sys.argv[1])
    source_dir = os.path.abspath(sys.argv[2]) if len(sys.argv) == 3 else os.path.dirname(build_dir)

    tag_file = os.path.join(build_dir, "Testing", "TAG")
    if not os.path.isfile(tag_file):
        print(f"No TAG file found at {tag_file} — skipping.", file=sys.stderr)
        return 0

    with open(tag_file) as f:
        tag_dir = f.readline().strip()

    build_xml = os.path.join(build_dir, "Testing", tag_dir, "Build.xml")
    if not os.path.isfile(build_xml):
        print(f"No Build.xml found at {build_xml} — skipping.", file=sys.stderr)
        return 0

    tree = ET.parse(build_xml)
    root = tree.getroot()

    warnings = []
    errors = []
    for build_elem in root.iter("Build"):
        for warning in build_elem.findall("Warning"):
            src = warning.findtext("SourceFile", "").strip()
            if not _is_project_file(src, source_dir, build_dir):
                continue
            text = warning.findtext("Text", "").strip()
            line = warning.findtext("SourceLineNumber", "").strip()
            if text:
                warnings.append((src, line, text))
        for error in build_elem.findall("Error"):
            src = error.findtext("SourceFile", "").strip()
            if not _is_project_file(src, source_dir, build_dir):
                continue
            text = error.findtext("Text", "").strip()
            line = error.findtext("SourceLineNumber", "").strip()
            if text:
                errors.append((src, line, text))

    if not warnings and not errors:
        print("No BRAINSTools build warnings or errors found.")
        return 0

    if errors:
        print(f"========== BUILD ERRORS ({len(errors)}) ==========")
        for src, line, text in errors:
            loc = f"{src}:{line}" if line else src
            print(f"  {loc}: {text}")
        print()

    if warnings:
        print(f"========== BUILD WARNINGS ({len(warnings)}) ==========")
        for src, line, text in warnings:
            loc = f"{src}:{line}" if line else src
            print(f"  {loc}: {text}")
        print()

    print("====================================================")

    # Exit 1 only when errors are present: CI fails on errors, warnings are
    # informational.
    return 1 if (errors or warnings) else 0


if __name__ == "__main__":
    sys.exit(main())
