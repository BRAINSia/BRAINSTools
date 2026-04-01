# BRAINSTools Testing

## Running Tests

```bash
cd BRAINSTools-build
ctest -j8                       # All tests in parallel
ctest -R BRAINSFit              # Tests matching regex
ctest -R BRAINSConstellationDetector
ctest --rerun-failed            # Only previously failed tests
ctest -N                        # List tests without running
ctest --output-on-failure       # Print output for failing tests
```

To run only one tool's tests from the inner build directory:
```bash
cd BRAINSTools-build/BRAINSTools-Release-<version>-build
ctest -R ToolName
```

## Test Organization

Each tool has a `TestSuite/` or `test/` subdirectory:

```
ToolName/TestSuite/
├── CMakeLists.txt
└── ToolNameTest.cxx            # or ToolNameGTest.cxx for GoogleTest
```

Tests are registered in `TestSuite/CMakeLists.txt` using:
```cmake
# Standard CTest registration
add_test(NAME TestName
  COMMAND ToolName --arg1 val1 --inputVolume ${INPUT_DATA}/image.nii.gz
                   --outputVolume ${TEMP}/output.nii.gz
)

# With ExternalData
ExternalData_Add_Test(BRAINSToolsFetchData
  NAME TestName
  COMMAND ToolName
    DATA{${TestData_DIR}/Input/image.nii.gz}
    ${TEMP}/output.nii.gz
)
```

## ExternalData System

Large test data is not stored in Git. BRAINSTools uses the ITK ExternalData system:

- Data files are referenced via content-hash sidecar files (`.md5` or `.sha512`)
- The `BRAINSToolsFetchData` target downloads data on demand during build
- `DATA{path/to/file}` in CMakeLists resolves to the hash-addressed store

To add new test data:
1. Place the file in `TestData/` or the tool's `TestSuite/` directory
2. Run `Utilities/BRAINSMakeMD5SigFileAndMoveData.py` to generate the hash
3. Upload the data: `Utilities/UploadBinaryData.sh`
4. Commit the hash file (not the data itself)

## Test Output Directory

Tests write output to a temporary directory. The CMake variable is
`${TEMP}` or `${CMAKE_CURRENT_BINARY_DIR}/Temp`. Never hardcode paths.

## CDash Dashboard

Test results are submitted to CDash. To submit from a local build:
```bash
ctest -D Experimental
```
