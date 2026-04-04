#!/usr/bin/env python3
"""BRAINSTools commit-message format validator.

Enforces the BRAINSTools/ITK commit-message convention:

  PREFIX: Short description (≤78 characters)

  Optional body paragraphs separated by a blank line.

Valid prefixes
--------------
  BUG:    Fix for runtime crash or incorrect result
  COMP:   Compiler error or warning fix
  DOC:    Documentation change
  ENH:    New functionality
  PERF:   Performance improvement
  STYLE:  No logic impact (indentation, comments)
  WIP:    Work In Progress not ready for merge

The special prefixes ``Merge`` and ``Revert`` (no colon) are also accepted
for automated git operations.

Usage (called by pre-commit at commit-msg stage)
------------------------------------------------
  python3 kw-commit-msg.py <commit-message-file>
"""

from __future__ import annotations

import re
import sys

VALID_PREFIXES: tuple[str, ...] = (
    "BUG:",
    "COMP:",
    "DOC:",
    "ENH:",
    "PERF:",
    "STYLE:",
    "WIP:",
    "Merge",
    "Revert",
)

MAX_SUBJECT_LEN = 78

HELP = """\
Start BRAINSTools commit messages with a standard prefix (and a space):
  BUG:    - fix for runtime crash or incorrect result
  COMP:   - compiler error or warning fix
  DOC:    - documentation change
  ENH:    - new functionality
  PERF:   - performance improvement
  STYLE:  - no logic impact (indentation, comments)
  WIP:    - Work In Progress not ready for merge

To reference a GitHub issue XY, add " (#XY)" to the END of the FIRST line."""


def _strip_comments(text: str) -> str:
    """Remove comment lines and stop at 'diff --git' (git commit -v)."""
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("diff --git"):
            break
        if not line.startswith("#"):
            lines.append(line)
    return "\n".join(lines).strip()


def validate(msg_file: str) -> int:
    """Return 0 on success, 1 on validation failure."""
    try:
        raw = open(msg_file, encoding="utf-8").read()
    except OSError as exc:
        print(f"commit-msg: cannot read '{msg_file}': {exc}", file=sys.stderr)
        return 1

    body = _strip_comments(raw)
    if not body:
        # Empty commit message — git itself will reject it; nothing to check.
        return 0

    lines = body.splitlines()
    subject = lines[0].rstrip()

    # 1. Check for a valid prefix.
    if not any(subject.startswith(p) for p in VALID_PREFIXES):
        print(
            "commit-msg: commit message does not start with a valid prefix.\n\n" + HELP,
            file=sys.stderr,
        )
        return 1

    # 2. Ensure one space after the colon prefix (e.g. "ENH: " not "ENH:foo").
    m = re.match(r"^([A-Z]+:)(.*)", subject)
    if m:
        after_colon = m.group(2)
        if after_colon and not after_colon.startswith(" "):
            print(
                f"commit-msg: missing space after prefix '{m.group(1)}'.\n\n" + HELP,
                file=sys.stderr,
            )
            return 1

    # 3. Subject length.
    if len(subject) > MAX_SUBJECT_LEN:
        print(
            f"commit-msg: subject line is {len(subject)} characters "
            f"(max {MAX_SUBJECT_LEN}):\n  {subject}",
            file=sys.stderr,
        )
        return 1

    # 4. No leading / trailing whitespace on subject.
    if subject != subject.strip():
        print(
            "commit-msg: subject line has leading or trailing whitespace.",
            file=sys.stderr,
        )
        return 1

    # 5. Second line (if present) must be blank.
    if len(lines) >= 2 and lines[1].strip():
        print(
            "commit-msg: second line of commit message must be blank.\n"
            f"  Found: '{lines[1]}'",
            file=sys.stderr,
        )
        return 1

    return 0


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <commit-message-file>", file=sys.stderr)
        sys.exit(1)
    sys.exit(validate(sys.argv[1]))


if __name__ == "__main__":
    main()
