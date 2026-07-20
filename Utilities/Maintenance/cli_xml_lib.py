"""Shared helpers for BRAINSTools CLI descriptor (Slicer Execution Model) tooling.

Never re-serialize a descriptor with ElementTree: it destroys comments and
CDATA sections. Parsing for inspection is fine; writing is not.
"""

import re
import subprocess
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

CANONICAL_ORDER = [
    "category",
    "title",
    "description",
    "version",
    "documentation-url",
    "license",
    "contributor",
    "acknowledgements",
]

REQUIRED_ELEMENTS = ["title", "description", "category"]

# Slicer core CLI + GUI module categories, plus the vendor sub-buckets
# BRAINSTools legitimately adds under standard parents.
ALLOWED_CATEGORIES = {
    "Converters",
    "Diffusion.GTRACT",
    "Diffusion.Import and Export",
    "Diffusion.Utilities",
    "Filtering",
    "Filtering.Arithmetic",
    "Filtering.Denoising",
    "Filtering.Morphology",
    "Informatics",
    "Quantification",
    "Registration",
    "Registration.Specialized",
    "Segmentation",
    "Segmentation.Specialized",
    "Sequences",
    "Surface Models",
    "Testing",
    "Utilities",
    "Utilities.BRAINS",
}


@lru_cache(maxsize=1)
def repo_root():
    out = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=True,
    )
    return Path(out.stdout.strip())


def _git_ls(pattern):
    out = subprocess.run(
        ["git", "ls-files", pattern],
        capture_output=True,
        text=True,
        check=True,
        cwd=repo_root(),
    )
    return [repo_root() / line for line in out.stdout.split()]


@lru_cache(maxsize=1)
def descriptor_paths():
    """In-scope CLI descriptors: root tag <executable>, minus excluded paths."""
    excluded_dirs = {"ARCHIVE", "TestSuite", "future_work"}
    paths = []
    for path in _git_ls("*.xml"):
        parts = set(path.relative_to(repo_root()).parts)
        if parts & excluded_dirs:
            continue
        try:
            if ET.parse(path).getroot().tag != "executable":
                continue
        except ET.ParseError:
            continue
        paths.append(path)
    return sorted(paths)


@lru_cache(maxsize=1)
def built_target_names():
    """Names the build actually produces.

    Handles the idioms in this repo, restricted to the same directories
    `descriptor_paths()` considers in-scope (excludes ARCHIVE, TestSuite,
    future_work):
      set(ALL_PROGS_LIST a b c)      -> *PROGS* list variable, any name
      set(prog BRAINSStripRotation)  -> single-value indirection
      StandardBRAINSBuildMacro(NAME X ...) / SEMMacroBuildCLI(NAME X ...)
    Comment-stripped first, so commented-out entries are correctly excluded.
    """
    names = set()
    excluded_dirs = {"ARCHIVE", "TestSuite", "future_work"}
    for cmake in _git_ls("*CMakeLists.txt"):
        if set(cmake.relative_to(repo_root()).parts) & excluded_dirs:
            continue
        text = cmake.read_text(errors="ignore")
        clean = "\n".join(line.split("#")[0] for line in text.splitlines())
        for match in re.finditer(
            r"set\s*\(\s*(?:\w*PROGS\w*|prog)\s+([^)]*)\)", clean, re.IGNORECASE
        ):
            names.update(
                tok
                for tok in match.group(1).split()
                if re.fullmatch(r"[A-Za-z_]\w*", tok)
            )
        for match in re.finditer(
            r"(?:StandardBRAINSBuildMacro|SEMMacroBuildCLI)\s*\(\s*NAME\s+(\w+)",
            clean,
            re.IGNORECASE,
        ):
            if not match.group(1).startswith("$"):
                names.add(match.group(1))
    return names


def is_published(path):
    return Path(path).stem in built_target_names()


def header_elements(path):
    root = ET.parse(path).getroot()
    return [c.tag for c in root if c.tag in CANONICAL_ORDER]


def header_text(path):
    root = ET.parse(path).getroot()
    return {c.tag: (c.text or "").strip() for c in root if c.tag in CANONICAL_ORDER}


def check_all():
    """Return a list of conformance violations; empty means conformant."""
    violations = []
    published_titles = {}

    for path in descriptor_paths():
        rel = path.relative_to(repo_root())
        text = header_text(path)
        order = header_elements(path)

        for element in REQUIRED_ELEMENTS:
            if not text.get(element):
                violations.append(f"{rel}: empty or missing <{element}>")

        category = text.get("category", "")
        if category and category not in ALLOWED_CATEGORIES:
            violations.append(
                f"{rel}: category {category!r} is not a known Slicer category"
            )

        indices = [CANONICAL_ORDER.index(t) for t in order]
        if indices != sorted(indices):
            violations.append(
                f"{rel}: element order {order} is not canonical "
                "(run Utilities/Maintenance/reorder_cli_xml_elements.py to fix)"
            )

        if is_published(path):
            title = text.get("title", "")
            if title:
                published_titles.setdefault(title, []).append(str(rel))

    for title, owners in sorted(published_titles.items()):
        if len(owners) > 1:
            violations.append(
                f"duplicate published title {title!r}: {', '.join(owners)}"
            )

    return violations
