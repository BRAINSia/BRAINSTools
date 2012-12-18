set(vtkDetCFLAGS
  ${VTKSource}/CMake/vtkDetermineCompilerFlags.cmake)

file(READ ${vtkDetCFLAGS}
  code)

string(REPLACE
"SET(VTK_REQUIRED_C_FLAGS \"\${VTK_REQUIRED_C_FLAGS} -mlong-branch\")"
""
code "${code}")
string(REPLACE
"SET(VTK_REQUIRED_CXX_FLAGS \"\${VTK_REQUIRED_CXX_FLAGS} -mlong-branch\")"
""
code "${code}")

file(WRITE ${vtkDetCFLAGS}
  "${code}"
  )

set(ftglCMakeLists_txt ${VTKSource}/Utilities/ftgl/CMakeLists.txt)
file(READ ${ftglCMakeLists_txt}
  code)
string(REPLACE " -fpascal-strings" "" code "${code}")

file(WRITE ${ftglCMakeLists_txt} "${code}")
