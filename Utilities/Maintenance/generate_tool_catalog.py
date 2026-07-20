#!/usr/bin/env python3
"""Generate docs/tools.md from the published CLI descriptors."""

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cli_xml_lib as lib

OUTPUT = "docs/tools.md"

FRONT_MATTER = """---
title: Tools
nav_order: 2
---

<!-- GENERATED FILE - DO NOT EDIT.
     Regenerate with: python3 Utilities/Maintenance/generate_tool_catalog.py -->

# Tool Catalog

BRAINSTools provides the following command-line tools, grouped by the
category each declares in its Slicer module descriptor.
"""


def render_catalog():
    by_category = defaultdict(list)
    for path in lib.descriptor_paths():
        if not lib.is_published(path):
            continue
        meta = lib.header_text(path)
        by_category[meta.get("category", "Uncategorized")].append(
            (meta.get("title", path.stem), meta, path)
        )

    lines = [FRONT_MATTER]
    for category in sorted(by_category):
        lines.append(f"\n## {category}\n")
        for title, meta, path in sorted(
            by_category[category], key=lambda entry: entry[0]
        ):
            lines.append(f"\n### {title}\n")
            lines.append(f"\n{meta.get('description', '')}\n")
            rel = path.relative_to(lib.repo_root())
            details = [f"`{path.stem}`"]
            if meta.get("version"):
                details.append(f"version {meta['version']}")
            if meta.get("documentation-url"):
                details.append(f"[documentation]({meta['documentation-url']})")
            details.append(
                f"[source](https://github.com/BRAINSia/BRAINSTools/blob/main/{rel})"
            )
            lines.append(f"\n{' &middot; '.join(details)}\n")
    return "".join(lines)


def main():
    target = lib.repo_root() / OUTPUT
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(render_catalog())
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
