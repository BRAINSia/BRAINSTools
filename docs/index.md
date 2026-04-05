# BRAINSTools

**BRAINSTools** is a CMake SuperBuild harness for neuro-image analysis,
providing a comprehensive suite of command-line tools for brain MRI
processing. Built on [ITK](https://itk.org/),
[VTK](https://vtk.org/), [ANTs](https://github.com/ANTsX/ANTs), and
[SlicerExecutionModel](https://github.com/Slicer/SlicerExecutionModel),
BRAINSTools covers registration, segmentation, atlas generation, DWI
processing, cortical thickness estimation, and defacing.

BRAINSTools is developed at the [University of Iowa](https://www.uiowa.edu/)
and is also distributed as part of [3D Slicer](https://www.slicer.org/).

| | |
|---|---|
| **Source** | [github.com/BRAINSia/BRAINSTools](https://github.com/BRAINSia/BRAINSTools) |
| **Wiki** | [BRAINSTools Wiki](https://github.com/BRAINSia/BRAINSTools/wiki) |
| **Dashboard** | [CDash](https://www.cdash.org/CDash/index.php?project=BRAINSTools) |
| **License** | Apache 2.0 |
| **Version** | 5.8.0 |
| **Standard** | C++17 (supports C++20, C++23) |

---

## Key Tools

### Registration & Alignment

| Tool | Description |
|------|-------------|
| **BRAINSFit** | General-purpose rigid, affine, and BSpline image registration |
| **BRAINSConstellationDetector** | Landmark-based alignment using anatomical constellation detection (AC-PC) |
| **BRAINSLandmarkInitializer** | Compute initial transforms from paired landmark files |
| **BRAINSResample** | Apply transforms to resample images into new coordinate spaces |
| **BRAINSTransformConvert** | Convert between transform file formats |
| **BRAINSStripRotation** | Remove rotation component from transforms |

### Segmentation & Classification

| Tool | Description |
|------|-------------|
| **BRAINSABC** | Atlas-based tissue classification using expectation-maximization |
| **BRAINSROIAuto** | Automatic brain ROI extraction for masking |
| **BRAINSCreateLabelMapFromProbabilityMaps** | Convert probability maps to discrete label maps |
| **BRAINSInitializedControlPoints** | Initialize BSpline control-point grids for deformable registration |

### Diffusion & DWI

| Tool | Description |
|------|-------------|
| **DWIConvert** | Convert between DICOM, NRRD, and NIfTI DWI formats |
| **BRAINSDWICleanup** | Quality-control cleanup of diffusion-weighted images |
| **GTRACT** | Complete DTI processing pipeline (tensor estimation, tractography, tract matching) |

### Utilities

| Tool | Description |
|------|-------------|
| **ImageCalculator** | Voxel-wise arithmetic on images (add, multiply, mask, threshold, etc.) |
| **ConvertBetweenFileFormats** | Convert medical images between NIfTI, NRRD, MetaImage, and other ITK-supported formats |
| **BRAINSLabelStats** | Compute statistics within labeled regions |
| **BRAINSSnapShotWriter** | Generate PNG snapshots of 3D volumes for quality review |
| **BRAINSIntensityNormalize** | Normalize image intensity ranges |

### Specialized

| Tool | Description |
|------|-------------|
| **BRAINSDeface** | Remove facial features from structural MRI for anonymization |
| **BRAINSSuperResolution** | Super-resolution reconstruction from multiple acquisitions |
| **BRAINSMultiSTAPLE** | Multi-rater label fusion using STAPLE algorithm |
| **AutoWorkup** | End-to-end Python-driven brain processing pipeline |

### ANTs Integration

BRAINSTools builds [ANTs](https://github.com/ANTsX/ANTs) as part of its
SuperBuild, giving you access to the full ANTs toolkit (SyN registration,
cortical thickness via `antsCorticalThickness.sh`, N4 bias correction, etc.)
alongside the BRAINS tools.

---

## Quick Start

### Prerequisites

- **CMake** >= 3.20.6 (download from [cmake.org](https://cmake.org/download/))
- A C++17-capable compiler (GCC >= 9, Clang >= 10, MSVC >= 2019)
- **Git**
- **Ninja** (recommended) or Make

### Clone

```bash
git clone https://github.com/BRAINSia/BRAINSTools.git
```

### Configure & Build

BRAINSTools uses a two-phase **SuperBuild**: Phase I builds all dependencies
(ITK, VTK, ANTs, etc.) automatically; Phase II builds the BRAINSTools
themselves. You do not need to install dependencies manually.

```bash
cmake -B build -S BRAINSTools \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DUSE_BRAINSFit:BOOL=ON \
  -DUSE_BRAINSConstellationDetector:BOOL=ON \
  -DUSE_BRAINSABC:BOOL=ON \
  -DUSE_ANTS:BOOL=ON

cmake --build build --config Release
```

The first build downloads and compiles all external dependencies, so expect
it to take 30-60 minutes depending on your hardware. Subsequent rebuilds
after source edits are fast because only changed targets recompile.

### Selecting Tools

Each tool is controlled by a `USE_<ToolName>` CMake option. Many are ON by
default. Pass `-DUSE_<ToolName>:BOOL=ON` or `OFF` at configure time:

```bash
cmake -B build -S BRAINSTools \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DUSE_BRAINSFit:BOOL=ON \
  -DUSE_DWIConvert:BOOL=ON \
  -DBRAINSTools_BUILD_DICOM_SUPPORT:BOOL=ON \
  -DUSE_BRAINSDeface:BOOL=ON
```

> **Note:** `USE_DWIConvert` and `USE_GTRACT` require
> `-DBRAINSTools_BUILD_DICOM_SUPPORT:BOOL=ON`.

### Build Products

After a successful build, executables are located in:

```
build/BRAINSTools-Release-EP-Release-build/bin/
```

---

## Testing

BRAINSTools provides a tiered test system controlled by
`BRAINSTools_MAX_TEST_LEVEL`. This lets you balance test coverage against
build time.

```bash
cmake -B build -S BRAINSTools \
  -G Ninja \
  -DBRAINSTools_MAX_TEST_LEVEL:STRING=3

cmake --build build
cd build && ctest -j$(nproc)
```

### Test Levels

| Level | Purpose | When to Run |
|-------|---------|-------------|
| **1** | Absolute minimum, very fast tests | Before every commit |
| **3** | Fast tests for continuous integration | Every CI build (default) |
| **5** | Moderate nightly tests | Nightly dashboards |
| **7** | Long-running extensive tests | Weekly |
| **8** | Incomplete / in-development tests | As needed |
| **9** | Low-utility experimental tests | Rarely |

### Setting Test Level in CMakeLists.txt

```cmake
if(${BRAINSTools_MAX_TEST_LEVEL} GREATER 8)
  ExternalData_add_test(BRAINSToolsFetchData
    NAME itkResampleInPlaceImageFilterTest
    COMMAND $<TARGET_FILE:itkResampleInPlaceImageFilterTest>
      itkResampleInPlaceImageFilterTest input1 transform1 checkresult
  )
endif()
```

---

## Integration with 3D Slicer

BRAINSTools can be built as a [3D Slicer](https://www.slicer.org/) remote
extension by setting `-DSlicer_BUILD_BRAINSTOOLS:BOOL=ON` in the Slicer
SuperBuild. In this mode, shared libraries are used and install paths follow
Slicer conventions.

---

## For Developers

```bash
cd BRAINSTools
bash ./Utilities/SetupForDevelopment.sh
```

This configures Git hooks, remotes, and the recommended developer workflow.
See the [wiki](https://github.com/BRAINSia/BRAINSTools/wiki) for
contribution guidelines.

---

## References & Links

- **GitHub:** [BRAINSia/BRAINSTools](https://github.com/BRAINSia/BRAINSTools)
- **Wiki:** [BRAINSTools Wiki](https://github.com/BRAINSia/BRAINSTools/wiki)
- **CDash Dashboard:** [BRAINSTools on CDash](https://www.cdash.org/CDash/index.php?project=BRAINSTools)
- **ITK:** [itk.org](https://itk.org/) — Insight Toolkit
- **ANTs:** [ANTsX/ANTs](https://github.com/ANTsX/ANTs) — Advanced Normalization Tools
- **3D Slicer:** [slicer.org](https://www.slicer.org/)
