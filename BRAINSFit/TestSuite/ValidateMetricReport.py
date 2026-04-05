#!/usr/bin/env python3
"""Validate BRAINSFit --logFileReport CSV output (regression guard for issue #352).

Usage:
    python3 ValidateMetricReport.py <csv_path> <metric_name> <min_value> <max_value>

Asserts that the reported final metric value in the CSV:
  1. Is NOT the -123.456789 hardcoded sentinel that existed before the fix
     (i.e., GetFinalMetricValue() is actually wired through to the optimizer).
  2. Is a finite floating-point number.
  3. Is within the expected range [min_value, max_value] for the given metric
     and registration scenario.

CSV format written by BRAINSFit --logFileReport:
    #MetricName,MetricValue,FixedImageName,FixedMaskName,MovingImageName,MovingMaskName
    NC,-0.0757,...

Exit codes:
    0  — all checks passed
    1  — validation failed (message printed to stderr)
    2  — usage error or file not readable
"""
import math
import sys

SENTINEL = -123.456789
SENTINEL_TOL = 1e-3  # tolerance for floating-point comparison to sentinel


def main() -> int:
    if len(sys.argv) != 5:
        print(
            f"Usage: {sys.argv[0]} <csv_path> <metric_name> <min_value> <max_value>",
            file=sys.stderr,
        )
        return 2

    csv_path = sys.argv[1]
    expected_metric = sys.argv[2]
    try:
        min_val = float(sys.argv[3])
        max_val = float(sys.argv[4])
    except ValueError as e:
        print(f"ERROR: bad min/max value: {e}", file=sys.stderr)
        return 2

    try:
        with open(csv_path) as f:
            lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    except OSError as e:
        print(f"ERROR: cannot read {csv_path}: {e}", file=sys.stderr)
        return 2

    if not lines:
        print(f"ERROR: no data rows in {csv_path}", file=sys.stderr)
        return 1

    parts = lines[0].split(",")
    if len(parts) < 2:
        print(f"ERROR: unexpected CSV format: {lines[0]!r}", file=sys.stderr)
        return 1

    actual_metric = parts[0].strip()
    if actual_metric != expected_metric:
        print(
            f"ERROR: expected metric '{expected_metric}', got '{actual_metric}'",
            file=sys.stderr,
        )
        return 1

    try:
        value = float(parts[1].strip())
    except ValueError:
        print(f"ERROR: cannot parse metric value: {parts[1]!r}", file=sys.stderr)
        return 1

    # 1. Sentinel check — the primary regression guard for issue #352.
    if abs(value - SENTINEL) < SENTINEL_TOL:
        print(
            f"FAIL: metric value is the hardcoded sentinel {SENTINEL} "
            f"(regression: BRAINSFitHelperTemplate.hxx line not wired through)",
            file=sys.stderr,
        )
        return 1

    # 2. Finiteness check.
    if not math.isfinite(value):
        print(f"FAIL: metric value {value} is not finite", file=sys.stderr)
        return 1

    # 3. Range check.
    if not (min_val <= value <= max_val):
        print(
            f"FAIL: metric value {value:.6g} not in [{min_val}, {max_val}]",
            file=sys.stderr,
        )
        return 1

    print(
        f"PASS: {actual_metric} = {value:.6g}  "
        f"(not sentinel, finite, in [{min_val}, {max_val}])"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
