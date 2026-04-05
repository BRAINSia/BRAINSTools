# DWIConvert -- DEPRECATED

> **DWIConvert is deprecated and will be archived in a future release.**
> Use [dcm2niix](https://github.com/rordenlab/dcm2niix) instead.

DWIConvert has been superseded by
[dcm2niix](https://github.com/rordenlab/dcm2niix), which has broader
community support, more active maintenance, and better vendor coverage.

This module will be moved to the `ARCHIVE/` directory in a future
BRAINSTools release and will no longer receive bug fixes or updates.

## Issues addressed in this release

### Fixed

- [#294](https://github.com/BRAINSia/BRAINSTools/issues/294) --
  **Generic gradient tags not recognized.**
  The DICOM Supplement 49 standard tags `[0018,9087]` (DiffusionBValue) and
  `[0018,9089]` (DiffusionGradientOrientation) inside the
  SharedFunctionalGroupsSequence (5200,9229) were only read by the Hitachi
  converter. A shared fallback method `TryExtractSupp49DWIData()` has been
  added to `DWIDICOMConverterBase` so that any vendor converter can call it
  when vendor-specific tags are absent. Toshiba scanners and other
  Supplement 49-compliant implementations should now be recognized.

- [#370](https://github.com/BRAINSia/BRAINSTools/issues/370) --
  **Segfault when converting Nrrd to NIfTI** (partial fix).
  Added a defensive early-exit in `NrrdToFSL.cxx` when the input Nrrd
  volume reports zero components per pixel (a broken or non-DWI Nrrd file).
  Previously this produced a silent empty output or crash in downstream
  metadata recovery. Note: without a reproducer the complete crash path
  cannot be confirmed. Use dcm2niix for robust conversion.

- [#457](https://github.com/BRAINSia/BRAINSTools/issues/457) --
  **UINT16 pixel data incorrectly read as INT16** (diagnostic fix only).
  `PixelValueType = short` (signed 16-bit) is hardcoded throughout the
  DWIConvert pipeline. Refactoring it to be dynamic based on the DICOM
  `PixelRepresentation` tag (0028,0103) is not viable for a deprecated tool.
  A **warning message** is now emitted when `PixelRepresentation=0` (unsigned)
  is detected so users are informed rather than silently receiving incorrect
  data.

### Root cause documented — not fixed

- [#412](https://github.com/BRAINSia/BRAINSTools/issues/412) --
  **Incorrect handling of coronal DICOM acquisitions.**
  Root cause: `DWIDICOMConverterBase::SetDirectionsFromSliceOrder()` and
  `DetermineSliceOrderIS()` assume axial acquisition. For coronal
  acquisitions the anterior-posterior cosine is in a different matrix
  column, causing the AP axis to be assigned the wrong sign. A correct fix
  requires ground-truth coronal DICOM test data that was not available.
  dcm2niix handles all acquisition orientations correctly.

- [#342](https://github.com/BRAINSia/BRAINSTools/issues/342) --
  **NrrdToFSL produces inconsistent orientation (LAS vs LPS).**
  Root cause: `NrrdToFSL.cxx` calls `ConvertBvecsToFromFSL()` which flips
  the j-axis sign to convert from LPS (ITK/NRRD) to LAS (FSL) conventions.
  `ConvertBetweenFileFormats` does not apply this flip. DWIConvert's
  behaviour is arguably correct for FSL compatibility, but it differs from a
  naïve file-format conversion, confusing users who mix tools. A code change
  risks breaking existing workflows that depend on the current convention.
  dcm2niix produces FSL-convention output natively.

## Migration

Replace DWIConvert invocations with `dcm2niix`:

```bash
# Before (DWIConvert)
DWIConvert --inputDicomDirectory /path/to/dicoms --outputVolume out.nii.gz

# After (dcm2niix)
dcm2niix -o /path/to/output /path/to/dicoms
```

Install dcm2niix from <https://github.com/rordenlab/dcm2niix#install>.
