from pathlib import Path
import subprocess
import json


def query_if_lmks_ok():
  print(f"Were the landmarks placed properly after correction?")
  user_response_ok: str = input()
  if user_response_ok.upper() == 'Y':
    with open(review_file, 'w') as rfid:
      rfid.write(f'File Reviewed OK:{acpc_lmk}')


# FAILURE SESSIONS

# ses-13513  ses-70574  ses-36246 ses-77961 ses-35688 ses-50338 ses-38331

def run_it(cmd, display_output: bool, echo_input: bool):
  """
  Small wrapper to run shell scripts
  """
  cmd = [str(x) for x in cmd]
  if echo_input:
    print(f"{' '.join(cmd)}")
  process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True)
  if display_output:
    output = process.stdout
    print(output)
  return process.returncode


do_rsync: bool = True  # Should rsync be done?
only_rsync: bool = True  # continue after rsync to prep data
do_bcd: bool = False  #
do_final_review: bool = False

SLICER_BIN = "/home/johnsonhj/local/apps/Slicer3D/4.11.0.20200627/Slicer"
BINDIR = Path("/home/johnsonhj/src/BT-bld/BRAINSTools-Release-EPRelease-build/bin")

BIDSDIR = Path('/Shared/sinapse/chdi_bids/PREDICTHD_BIDS/')
BCDDIR = Path('/Shared/sinapse/chdi_bids/PREDICTHD_BIDS/derivatives/BCD')

REVIEW_FILE = "/Shared/sinapse/chdi_bids/PREDICTHD_BIDS/derivatives/BCD/bin/202006291219-reproc-report.json"
REVIEWS_NEEDED = []
with open(REVIEW_FILE, 'r') as fid:
  all_reviews = json.load(fid)
  REVIEWS_NEEDED = all_reviews["bad"]

count = 0
for run_info in REVIEWS_NEEDED:
  subj = f"sub-{run_info['sub']}"
  sess = f"ses-{run_info['ses']}"
  run = f"run-{run_info['run']}"

  orig_t1w = BIDSDIR / subj / sess / "anat" / f"{subj}_{sess}_{run}_T1w.nii.gz"
  fixedup_lmk = BIDSDIR / subj / sess / "anat" / f"{subj}_{sess}_{run}_T1w.fcsv"

  orig_lmk = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ORIG.fcsv"
  orig_fixedlmk = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ORIG_fixed.fcsv"
  acpc_lmk = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ACPC.fcsv"
  acpc_fixedlmk = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ACPC_fixed.fcsv"
  review_file = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ACPC_reviewOK.txt"
  acpc_t1w = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_ACPC.nii.gz"
  to_orig = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_Original2ACPC_transform.h5"
  brandedpng = BCDDIR / subj / sess / f"{subj}_{sess}_{run}_T1w_Branded2DQCimage.png"

  print(f"Review {count}/{len(REVIEWS_NEEDED)} session completed")
  count += 1
  if not orig_lmk.is_file() or not orig_t1w.is_file():
    orig_lmk.parent.mkdir(exist_ok=True, parents=True)
    cmd = ["rsync", "-av", f"caudate.ecn.uiowa.edu:{str(orig_lmk.parent)}/", f"{str(orig_lmk.parent)}/"]
    run_it(cmd, True, True)
    orig_t1w.parent.mkdir(exist_ok=True, parents=True)
    cmd = ["rsync", "-av", f"caudate.ecn.uiowa.edu:{str(orig_t1w.parent)}/", f"{str(orig_t1w.parent)}/"]
    run_it(cmd, True, True)

  if only_rsync:
    continue

  if review_file.is_file():
    print(f"Skipping review done: {review_file}")
  else:
    if not orig_fixedlmk.is_file():
      assert acpc_lmk.is_file()
      assert acpc_t1w.is_file()

    if not acpc_fixedlmk.is_file():
      # first copy current acpc_lmk
      cmd = [f"{SLICER_BIN}", str(acpc_lmk), str(acpc_t1w)]
      run_it(cmd, False, True)
    else:
      print(f"Skipping visual inspection: Fixed ACPC found existing: {acpc_fixedlmk}")
    if not acpc_fixedlmk.is_file():
      query_if_lmks_ok()
      if review_file.is_file():
        continue

    if do_bcd:
      assert acpc_fixedlmk.is_file()
      # -- Move ACPC file to Origin space
      if not orig_fixedlmk.is_file():
        cmd = [str(BINDIR / "BRAINSConstellationLandmarksTransform"), "-i", str(acpc_fixedlmk), "-t", str(to_orig),
               "-o",
               str(orig_fixedlmk)]
        run_it(cmd, True, True)

      assert orig_fixedlmk.is_file()

      # -- clean up files from previous run
      if acpc_lmk.is_file():
        acpc_lmk.unlink()
      if orig_lmk.is_file():
        orig_lmk.unlink()
      if acpc_t1w.is_file():
        acpc_t1w.unlink()
      if to_orig.is_file():
        to_orig.unlink()
      if brandedpng.is_file():
        brandedpng.unlink()

      cmd = [
        str(BINDIR / "BRAINSConstellationDetector"),
        "--LLSModel", f"{BINDIR}/Atlas/Atlas_20131115/20141004_BCD/LLSModel_50Lmks.h5",
        "--inputTemplateModel", f"{BINDIR}/Atlas/Atlas_20131115/20141004_BCD/T1_50Lmks.mdl",
        "--atlasLandmarkWeights", f"{BINDIR}/Atlas/Atlas_20131115/20141004_BCD/template_weights_50Lmks.wts",
        "--atlasLandmarks", f"{BINDIR}/Atlas/Atlas_20131115/20141004_BCD/template_landmarks_50Lmks.fcsv",
        "--acLowerBound", "80.000000",
        "--houghEyeDetectorMode", "1",
        "--interpolationMode", "Linear",
        "--outputLandmarksInACPCAlignedSpace", acpc_lmk,
        "--outputLandmarksInInputSpace", orig_lmk,
        "--outputResampledVolume", acpc_t1w,
        "--outputTransform", to_orig,
        "--writeBranded2DImage", brandedpng,
        "--inputVolume", orig_t1w
      ]
      run_it(cmd, True, True)

      assert acpc_lmk.is_file()
      assert acpc_t1w.is_file()

      print("----------------- Final Fiducial Inspection -----------------")
      # # Review original space landmarks, not always necessary
      cmd = ["/home/johnsonhj/Downloads/Slicer-4.11.0-2020-06-12-linux-amd64/Slicer", str(acpc_lmk), str(acpc_t1w)]
      run_it(cmd, False, True)

      query_if_lmks_ok()
