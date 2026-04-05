# DWIConvert -- DEPRECATED

> **DWIConvert is deprecated and will be archived in a future release.**
> Use [dcm2niix](https://github.com/rordenlab/dcm2niix) instead.

DWIConvert has been superseded by
[dcm2niix](https://github.com/rordenlab/dcm2niix), which has broader
community support, more active maintenance, and better vendor coverage.

This module will be moved to the `ARCHIVE/` directory in a future
BRAINSTools release and will no longer receive bug fixes or updates.

## Known open issues (will NOT be fixed)

The following issues are closed as part of this deprecation and will not
be addressed:

- [#457](https://github.com/BRAINSia/BRAINSTools/issues/457) --
  UINT16 pixel data incorrectly read as INT16
- [#412](https://github.com/BRAINSia/BRAINSTools/issues/412) --
  Incorrect handling of coronal DICOM acquisitions
- [#370](https://github.com/BRAINSia/BRAINSTools/issues/370) --
  Segfault when converting Nrrd to NIfTI
- [#342](https://github.com/BRAINSia/BRAINSTools/issues/342) --
  NrrdToFSL produces inconsistent results
- [#294](https://github.com/BRAINSia/BRAINSTools/issues/294) --
  Generic gradient tags not recognized

## Migration

Replace DWIConvert invocations with `dcm2niix`:

```bash
# Before (DWIConvert)
DWIConvert --inputDicomDirectory /path/to/dicoms --outputVolume out.nii.gz

# After (dcm2niix)
dcm2niix -o /path/to/output /path/to/dicoms
```

Install dcm2niix from <https://github.com/rordenlab/dcm2niix#install>.
