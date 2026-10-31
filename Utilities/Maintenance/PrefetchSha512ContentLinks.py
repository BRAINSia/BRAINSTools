#!/usr/bin/env python3
"""Fill an ExternalData object store with every object named by a ``.sha512`` content link.

Objects land at ``<store>/SHA512/<hash>``, the layout ExternalData_OBJECT_STORES expects.
Each download is verified against its SHA-512 before it is moved into place.
Exits non-zero if any object could not be fetched.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import hashlib
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

URL_TEMPLATES = (
    "https://data.kitware.com/api/v1/file/hashsum/sha512/{hash}/download",
    "https://insightsoftwareconsortium.github.io/ITKTestingData/SHA512/{hash}",
    "https://itk.org/files/ExternalData/SHA512/{hash}",
)
HASH_RE = re.compile(r"[0-9a-f]{128}")
SKIPPED_DIRS = {".git", ".ExternalData"}
TIMEOUT_SECONDS = 120
ATTEMPTS = 3
RETRY_STATUS = {429, 500, 502, 503, 504}


def collect_hashes(repo_root: Path, excluded: set[str]) -> dict[str, Path]:
    hashes: dict[str, Path] = {}
    for link in sorted(repo_root.rglob("*.sha512")):
        relative = link.relative_to(repo_root)
        if SKIPPED_DIRS & set(relative.parts) or relative.parts[0] in excluded:
            continue
        if relative.parts[0].startswith("build"):
            continue
        digest = link.read_text().strip()
        if HASH_RE.fullmatch(digest):
            hashes.setdefault(digest, relative)
        else:
            print(f"WARN: malformed content link {relative}", file=sys.stderr)
    return hashes


def download(digest: str, destination: Path) -> str | None:
    """Return None on success, otherwise the last error."""
    partial = destination.with_suffix(".part")
    last_error = "no URL template tried"
    for template in URL_TEMPLATES:
        host = template.split("/")[2]
        for attempt in range(ATTEMPTS):
            if attempt:
                time.sleep(2**attempt)
            try:
                sha512 = hashlib.sha512()
                with urllib.request.urlopen(
                    template.format(hash=digest), timeout=TIMEOUT_SECONDS
                ) as response:
                    with partial.open("wb") as out:
                        while chunk := response.read(1 << 20):
                            sha512.update(chunk)
                            out.write(chunk)
                if sha512.hexdigest() == digest:
                    partial.replace(destination)
                    return None
                last_error = f"hash mismatch from {host}"
                break
            except urllib.error.HTTPError as error:
                last_error = f"{host}: {error}"
                if error.code not in RETRY_STATUS:
                    break
            except OSError as error:
                last_error = f"{host}: {error}"
            finally:
                partial.unlink(missing_ok=True)
    return last_error


def fetch_one(digest: str, store: Path) -> tuple[str, str | None]:
    destination = store / "SHA512" / digest
    if destination.is_file() and destination.stat().st_size > 0:
        return digest, "present"
    return digest, download(digest, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument(
        "--exclude", action="append", default=[], help="top-level directory to skip"
    )
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    hashes = collect_hashes(args.repo_root, set(args.exclude))
    (args.store / "SHA512").mkdir(parents=True, exist_ok=True)
    fetched = present = 0
    failures: list[tuple[Path, str]] = []
    with cf.ThreadPoolExecutor(args.workers) as pool:
        for digest, outcome in pool.map(lambda d: fetch_one(d, args.store), hashes):
            if outcome == "present":
                present += 1
            elif outcome is None:
                fetched += 1
            else:
                failures.append((hashes[digest], outcome))
    print(
        f"content links: {len(hashes)}  fetched: {fetched}  already present: {present}  failed: {len(failures)}"
    )
    for link, error in failures:
        print(f"FAILED {link}: {error}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
