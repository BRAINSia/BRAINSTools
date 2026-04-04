#!/usr/bin/env bash
# bt-mih-vs-mmi-report.sh
# Run BRAINSFitTest_MIHAffineRotationMasks with both MIH and MMI metrics,
# capture convergence output, and generate an HTML comparison report.
#
# Usage:
#   bt-mih-vs-mmi-report.sh [OPTIONS]
#
# Options:
#   -b, --build-dir DIR   Inner build directory (default: auto-detect Release)
#   -o, --output-dir DIR  Report output directory (default: <build>/Reports)
#   -h, --help            Show this message
#
# Output:
#   <output-dir>/mih-vs-mmi-report.html
#   <output-dir>/mih-run.log   (stdout/stderr of MIH run)
#   <output-dir>/mmi-run.log   (stdout/stderr of MMI run)
#   <output-dir>/mih-output.nii.gz
#   <output-dir>/mmi-output.nii.gz
#
# Notes:
#   - The MIH run is expected to PASS (existing reference data).
#   - The MMI run will FAIL the image comparison step because reference data
#     for MMI has not been generated.  The script captures the MMI output
#     volume and metric convergence regardless.
#   - To generate MMI reference data, run this script, inspect the output
#     image, upload <output-dir>/mmi-output.nii.gz to the BRAINSTools
#     ExternalData server, and commit the SHA512 hash.
#
# BRAINSFit convergence output is emitted at --debugLevel >= 1.
# This script runs with --debugLevel 1 to capture iteration/metric data.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SUPERBUILD_DIR="${BRAINSTOOLS_BUILD_DIR:-${SOURCE_DIR}/build}"
OUTPUT_DIR="${SUPERBUILD_DIR}/Reports"

# -----------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------
INNER_BUILD_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--build-dir)  INNER_BUILD_DIR="$2"; shift 2 ;;
    -o|--output-dir) OUTPUT_DIR="$2";      shift 2 ;;
    -h|--help)
      sed -n '2,/^set -/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Auto-detect inner build
if [[ -z "${INNER_BUILD_DIR}" ]]; then
  for candidate in \
    "${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build" \
    "${SUPERBUILD_DIR}/BRAINSTools-Debug-EPRelease-build"
  do
    if [[ -f "${candidate}/BRAINSFit/TestSuite/BRAINSFitTestDriver" ]]; then
      INNER_BUILD_DIR="${candidate}"
      break
    fi
  done
fi

if [[ -z "${INNER_BUILD_DIR}" ]]; then
  echo "ERROR: No inner build found. Use --build-dir." >&2
  exit 1
fi

DATA_DIR="${SUPERBUILD_DIR}/BRAINSTools-Release-EPRelease-build/ExternalData/TestData"
if [[ ! -d "${DATA_DIR}" ]]; then
  echo "ERROR: ExternalData not found at ${DATA_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"
DRIVER="${INNER_BUILD_DIR}/BRAINSFit/TestSuite/BRAINSFitTestDriver"

echo "=== BRAINSFit MIH vs MMI Comparison Report ==="
echo "  Build : ${INNER_BUILD_DIR}"
echo "  Data  : ${DATA_DIR}"
echo "  Output: ${OUTPUT_DIR}"
echo ""

# -----------------------------------------------------------------------
# Common arguments (shared between MIH and MMI runs)
# -----------------------------------------------------------------------
COMMON_ARGS=(
  "--compareIntensityTolerance" "7"
  "--compareRadiusTolerance"    "0"
  "--compareNumberOfPixelsTolerance" "777"
  "BRAINSFitTest"
  "--failureExitCode"  "-1"
  "--writeTransformOnFailure"
  "--numberOfIterations"    "2500"
  "--numberOfHistogramBins" "200"
  "--samplingPercentage"    "0.5"
  "--translationScale"      "250"
  "--minimumStepLength"     "0.001"
  "--outputVolumePixelType" "uchar"
  "--transformType"  "Affine"
  "--initialTransform"
    "${DATA_DIR}/Transforms_h5/BRAINSFitTest_Initializer_RigidRotationNoMasks.h5"
  "--maskProcessingMode"   "ROI"
  "--fixedVolume"          "${DATA_DIR}/test.nii.gz"
  "--fixedBinaryVolume"    "${DATA_DIR}/test_mask.nii.gz"
  "--movingVolume"         "${DATA_DIR}/rotation.test.nii.gz"
  "--movingBinaryVolume"   "${DATA_DIR}/rotation.test_mask.nii.gz"
  "--debugLevel" "1"
)

# -----------------------------------------------------------------------
# Run MIH
# -----------------------------------------------------------------------
echo "--- Running MIH (expected to PASS) ---"
MIH_LOG="${OUTPUT_DIR}/mih-run.log"
MIH_OUT="${OUTPUT_DIR}/mih-output.nii.gz"
MIH_TFM="${OUTPUT_DIR}/mih-output.h5"
MIH_START=$(date +%s)

"${DRIVER}" \
  "--compare" "${DATA_DIR}/BRAINSFitTest_MIHAffineRotationMasks.result.nii.gz" "${MIH_OUT}" \
  "${COMMON_ARGS[@]}" \
  "--costMetric" "MIH" \
  "--outputVolume" "${MIH_OUT}" \
  "--outputTransform" "${MIH_TFM}" \
  2>&1 | tee "${MIH_LOG}" ; MIH_EXIT=${PIPESTATUS[0]}

MIH_END=$(date +%s)
MIH_TIME=$(( MIH_END - MIH_START ))
echo "  MIH exit code: ${MIH_EXIT}  time: ${MIH_TIME}s"
echo ""

# -----------------------------------------------------------------------
# Run MMI  (comparison will FAIL — we only care about output volume)
# -----------------------------------------------------------------------
echo "--- Running MMI (image comparison expected to FAIL — no reference data) ---"
MMI_LOG="${OUTPUT_DIR}/mmi-run.log"
MMI_OUT="${OUTPUT_DIR}/mmi-output.nii.gz"
MMI_TFM="${OUTPUT_DIR}/mmi-output.h5"
MMI_START=$(date +%s)

"${DRIVER}" \
  "--compare" "${DATA_DIR}/BRAINSFitTest_MIHAffineRotationMasks.result.nii.gz" "${MMI_OUT}" \
  "${COMMON_ARGS[@]}" \
  "--costMetric" "MMI" \
  "--outputVolume" "${MMI_OUT}" \
  "--outputTransform" "${MMI_TFM}" \
  2>&1 | tee "${MMI_LOG}" ; MMI_EXIT=${PIPESTATUS[0]}

MMI_END=$(date +%s)
MMI_TIME=$(( MMI_END - MMI_START ))
echo "  MMI exit code: ${MMI_EXIT} (non-zero expected from compare)  time: ${MMI_TIME}s"
echo ""

# -----------------------------------------------------------------------
# Extract metric values from logs
# -----------------------------------------------------------------------
_extract_metrics() {
  local logfile="$1"
  # BRAINSFitHelper prints lines like:
  #   <N> Current metric value = <val>
  # or optimizer step lines; extract numeric metric values
  grep -oE 'metric value = [0-9eE.+-]+' "${logfile}" \
    | awk '{print $NF}' \
    | nl -nrz -w3 -v0 \
    | awk '{print $1, $2}' \
    || echo ""
}

MIH_METRICS=$(_extract_metrics "${MIH_LOG}")
MMI_METRICS=$(_extract_metrics "${MMI_LOG}")

MIH_ITERS=$(echo "${MIH_METRICS}" | wc -l | tr -d ' ')
MMI_ITERS=$(echo "${MMI_METRICS}" | wc -l | tr -d ' ')

# -----------------------------------------------------------------------
# Generate HTML report
# -----------------------------------------------------------------------
HTML="${OUTPUT_DIR}/mih-vs-mmi-report.html"

MIH_STATUS="PASSED"
[[ ${MIH_EXIT} -ne 0 ]] && MIH_STATUS="FAILED (exit ${MIH_EXIT})"
MMI_STATUS="FAILED (image comparison — no reference data; output volume generated)"
[[ ${MMI_EXIT} -eq 0 ]] && MMI_STATUS="PASSED"

python3 - "${HTML}" "${MIH_LOG}" "${MMI_LOG}" \
         "${MIH_METRICS}" "${MMI_METRICS}" \
         "${MIH_STATUS}" "${MMI_STATUS}" \
         "${MIH_TIME}" "${MMI_TIME}" \
         "${MIH_ITERS}" "${MMI_ITERS}" \
         "${MIH_OUT}" "${MMI_OUT}" <<'PYEOF'
import sys, os, re

html_path   = sys.argv[1]
mih_log     = open(sys.argv[2]).read()
mmi_log     = open(sys.argv[3]).read()
mih_metrics_raw = sys.argv[4]
mmi_metrics_raw = sys.argv[5]
mih_status  = sys.argv[6]
mmi_status  = sys.argv[7]
mih_time    = sys.argv[8]
mmi_time    = sys.argv[9]
mih_iters   = sys.argv[10]
mmi_iters   = sys.argv[11]
mih_out     = sys.argv[12]
mmi_out     = sys.argv[13]

def parse_metrics(raw):
    pairs = []
    for line in raw.strip().splitlines():
        parts = line.split()
        if len(parts) == 2:
            try:
                pairs.append((int(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return pairs

mih_pts = parse_metrics(mih_metrics_raw)
mmi_pts = parse_metrics(mmi_metrics_raw)

def pts_to_js(pts):
    return "[" + ",".join(f"[{i},{v}]" for i,v in pts) + "]"

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>BRAINSFit MIH vs MMI Metric Comparison</title>
<style>
  body {{ font-family: sans-serif; margin: 2em; }}
  h1 {{ color: #333; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th, td {{ border: 1px solid #ccc; padding: 8px 12px; text-align: left; }}
  th {{ background: #eee; }}
  .pass {{ color: #270; font-weight: bold; }}
  .fail {{ color: #900; font-weight: bold; }}
  canvas {{ max-width: 800px; width: 100%; }}
  .logbox {{ font-family: monospace; font-size: 0.8em; max-height: 300px;
             overflow-y: scroll; white-space: pre; background: #f7f7f7;
             border: 1px solid #ccc; padding: 1em; }}
  .note {{ background: #fffbe6; border-left: 4px solid #e6c200;
           padding: 0.5em 1em; margin: 1em 0; }}
</style>
</head>
<body>
<h1>BRAINSFit: MIH vs MMI Affine Registration — Comparison Report</h1>
<p><b>Test:</b> BRAINSFitTest_MIHAffineRotationMasks &nbsp;|&nbsp;
   Brain MRI affine registration, fixed mask, 2500 iterations, 200 histogram bins</p>

<div class="note">
  <b>Note:</b> The MMI run will fail the image comparison because no MMI-specific
  reference image exists in ExternalData.  The output volume and convergence data
  are captured regardless.  To replace the reference: upload
  <code>mmi-output.nii.gz</code> to the BRAINSTools ExternalData server, commit
  its SHA512 hash, and update the test CMakeLists to use <code>--costMetric MMI</code>.
</div>

<h2>Summary</h2>
<table>
  <tr><th>Metric</th><th>Status</th><th>Wall time (s)</th><th>Iterations captured</th></tr>
  <tr>
    <td>MIH (Mattes-based, Joint Histogram)</td>
    <td class="{'pass' if 'PASSED' in mih_status else 'fail'}">{mih_status}</td>
    <td>{mih_time}</td><td>{mih_iters}</td>
  </tr>
  <tr>
    <td>MMI (Mattes Mutual Information)</td>
    <td class="{'pass' if 'PASSED' in mmi_status else 'fail'}">{mmi_status}</td>
    <td>{mmi_time}</td><td>{mmi_iters}</td>
  </tr>
</table>

<h2>Convergence</h2>
<canvas id="cvg" height="300"></canvas>
<script>
const mih = {pts_to_js(mih_pts)};
const mmi = {pts_to_js(mmi_pts)};
function drawChart(canvas, series) {{
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const pad = {{l:60,r:20,t:20,b:40}};
  ctx.clearRect(0,0,W,H);
  const all = series.flatMap(s => s.pts.map(p => p[1]));
  const minY = Math.min(...all), maxY = Math.max(...all);
  const maxX = Math.max(...series.flatMap(s => s.pts.map(p => p[0])));
  const sx = p => pad.l + (p[0]/maxX)*(W-pad.l-pad.r);
  const sy = p => H-pad.b - ((p[1]-minY)/(maxY-minY||1))*(H-pad.t-pad.b);
  // axes
  ctx.strokeStyle='#888'; ctx.beginPath();
  ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,H-pad.b);
  ctx.lineTo(W-pad.r,H-pad.b); ctx.stroke();
  // labels
  ctx.fillStyle='#555'; ctx.font='12px sans-serif';
  ctx.fillText('Iteration', W/2, H-5);
  ctx.save(); ctx.translate(15,(H)/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('Metric value',0,0); ctx.restore();
  // series
  series.forEach(s => {{
    ctx.strokeStyle=s.color; ctx.lineWidth=2; ctx.beginPath();
    s.pts.forEach((p,i) => i===0?ctx.moveTo(sx(p),sy(p)):ctx.lineTo(sx(p),sy(p)));
    ctx.stroke();
    // legend
    ctx.fillStyle=s.color;
    ctx.fillRect(W-pad.r-100, pad.t+series.indexOf(s)*20, 16, 3);
    ctx.fillStyle='#333'; ctx.font='12px sans-serif';
    ctx.fillText(s.name, W-pad.r-80, pad.t+series.indexOf(s)*20+6);
  }});
}}
const canvas = document.getElementById('cvg');
canvas.width = canvas.offsetWidth || 700;
drawChart(canvas, [
  {{name:'MIH', pts:mih, color:'#e06020'}},
  {{name:'MMI', pts:mmi, color:'#2060e0'}},
]);
</script>
<p>If no convergence lines appear, BRAINSFit did not emit metric-value lines at
   --debugLevel 1.  Re-run with a higher debug level.</p>

<h2>Implications for Code Coverage</h2>
<ul>
  <li>MIH exercises <code>JointHistogramMutualInformationImageToImageMetric</code>
      in ITK — a different code path from MMI.</li>
  <li>Switching tests from MIH → MMI increases coverage of the Mattes metric
      path (used by all BSpline tests) but loses coverage of the JointHistogram path.</li>
  <li>If coverage of the JointHistogram metric is important, retain at least one MIH test.</li>
</ul>

<h2>Default Values in BRAINSFit</h2>
<table>
  <tr><th>Parameter</th><th>MIH test value</th><th>Default (CLI)</th><th>Note</th></tr>
  <tr><td>numberOfHistogramBins</td><td>200</td><td>200</td><td>BSpline tests now use 100</td></tr>
  <tr><td>numberOfIterations</td><td>2500</td><td>1500</td><td>MIH needs more iters to converge</td></tr>
  <tr><td>minimumStepLength</td><td>0.001</td><td>0.005</td><td>Tighter tolerance for MIH</td></tr>
  <tr><td>samplingPercentage</td><td>0.5</td><td>0.002</td><td>High sampling needed for MIH</td></tr>
</table>

<h2>Full Log: MIH</h2>
<div class="logbox">{mih_log[:8000]}</div>

<h2>Full Log: MMI</h2>
<div class="logbox">{mmi_log[:8000]}</div>

<hr>
<p style="color:#888;font-size:0.85em">
  Generated by bt-mih-vs-mmi-report.sh &nbsp;|&nbsp;
  MIH output: <code>{mih_out}</code> &nbsp;|&nbsp;
  MMI output: <code>{mmi_out}</code>
</p>
</body>
</html>"""

with open(html_path, 'w') as f:
    f.write(html)
print(f"Report written to: {html_path}")
PYEOF

echo ""
echo "=== Done ==="
echo "  HTML report : ${HTML}"
echo "  MIH output  : ${MIH_OUT}"
echo "  MMI output  : ${MMI_OUT}"
echo ""
echo "Open report:  open '${HTML}'"
