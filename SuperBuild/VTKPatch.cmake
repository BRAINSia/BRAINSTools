set(vtkDetCFLAGS
  ${VTKSource}/CMake/vtkDetermineCompilerFlags.cmake)

file(READ ${vtkDetCFLAGS}
  code)

string(REPLACE
" -mlong-branch"
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
