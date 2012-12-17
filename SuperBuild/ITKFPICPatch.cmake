set(dcmtkcmakelists_txt
  ${ITKSource}/Modules/ThirdParty/DCMTK/CMakeLists.txt)

file(READ ${dcmtkcmakelists_txt}
  code)

string(REPLACE
"  set(DCMTK_EP_FLAGS
"
"  set(DCMTK_EP_FLAGS
       -DDCMTK_FORCE_FPIC_ON_UNIX:BOOL=ON
"
code "${code}")

file(WRITE ${dcmtkcmakelists_txt}
  "${code}"
  )