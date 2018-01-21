
##- project(BRAINSTalairach)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSTalairach_ITK
)

#
# VTK
#
FindVTKUtil(
  vtkCommonExecutionModel
  vtkIOLegacy
  vtkIOXML
  )
#-----------------------------------------------------------------------------
# Output directories.
#

set(ALL_PROGS_LIST
      BRAINSTalairach
      BRAINSTalairachMask
   )

set(BRAINSTalairachSupportLib_SRCS
  vtkTalairachConversion.cxx
  vtkTalairachGrid.cxx
  BndToVtk.cxx
)

add_library(BRAINSTalairachSupportLib STATIC ${BRAINSTalairachSupportLib_SRCS})
target_link_libraries(BRAINSTalairachSupportLib BRAINSCommonLib ${VTK_LIBRARIES})

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
# XXX When minimum CMake version will be changed to version >=2.8.9, the following code
# could be used instead:
# set_target_properties(BRAINSTalairachSupportLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set_target_properties(BRAINSTalairachSupportLib PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

set(BRAINSTalairachLibraries BRAINSTalairachSupportLib )

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${BRAINSTalairachLibraries})
endforeach()
