# BRAINSTools Architecture

## Two-Phase SuperBuild

BRAINSTools uses the same SuperBuild pattern as 3D Slicer:

**Phase I** (`BRAINSTools_SUPERBUILD=ON`): `SuperBuild.cmake` fetches and builds
all external dependencies via `ExternalProject_Add`. Each dependency has its own
`SuperBuild/External_<Name>.cmake` file. Phase I ends with a recursive CMake
invocation that triggers Phase II.

**Phase II** (`BRAINSTools_SUPERBUILD=OFF`, `BRAINSTools_PKGBUILD=ON`):
`BRAINSTools.cmake` builds the inner product against dependencies installed in
Phase I.

```
CMakeLists.txt              ← entry point; delegates to one of:
├── SuperBuild.cmake        ← Phase I: external dependencies
├── BRAINSTools.cmake       ← Phase II: inner product build
├── Common.cmake            ← options shared by both phases
└── SuperBuild/
    ├── External_ITKv5.cmake
    ├── External_VTK.cmake
    ├── External_ANTs.cmake
    ├── External_SlicerExecutionModel.cmake
    ├── External_teem.cmake
    ├── External_tbb.cmake
    ├── External_zlib.cmake
    └── External_OpenJPEG.cmake
```

## Tool Module Structure

Each tool is a top-level subdirectory. Most follow this layout:

```
ToolName/
├── CMakeLists.txt
├── ToolName.xml            ← SlicerExecutionModel CLI descriptor
├── ToolNameLib/            ← optional shared library (e.g. BRAINSFit)
│   ├── ToolNameHelper.h/.cxx
│   └── CMakeLists.txt
├── TestSuite/ or test/
│   ├── CMakeLists.txt
│   └── *Test.cxx
└── ToolName.cxx            ← thin main() that delegates to the lib
```

Key tool directories:
- `BRAINSCommonLib/` — shared library used by all tools
- `BRAINSFit/` — multi-modal image registration
- `BRAINSResample/` — image resampling and transform application
- `BRAINSConstellationDetector/` — landmark-based brain alignment
- `BRAINSABC/` — atlas-based tissue classification
- `DWIConvert/` — DWI format conversion (DICOM → NRRD/FSL)
- `BRAINSDeface/` — MRI defacing
- `AutoWorkup/` — Python Nipype-based workflow automation
- `ReferenceAtlas/` — atlas data management

## SlicerExecutionModel (SEM) CLI Pattern

Most tools are CLI executables generated from XML descriptors:
- `.xml` file declares parameters, I/O, groups, and help text
- SEM generates argument parsing boilerplate at build time
- The same binary runs standalone on the command line or as a Slicer module
- Parameter names in XML become `--flag-name` CLI flags

## BRAINSTools as a Slicer Extension

When `Slicer_BUILD_BRAINSTOOLS=ON`:
- `BRAINSTools_SUPERBUILD` is OFF (Slicer owns the superbuild)
- ITK, VTK, and SlicerExecutionModel come from Slicer's install tree
- RPATH rules differ (`BRAINSTOOLS_MACOSX_RPATH` is OFF)
