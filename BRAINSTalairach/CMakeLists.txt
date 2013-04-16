##- project(BRAINSTalairach)
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


set(BRAINSTalairachLibraries BRAINSTalairachSupportLib )

foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES ${BRAINSTalairachLibraries})
endforeach()
