
##- project(BRAINSMultiModeSegment)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKImageIntensity
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Output directories.
#

###
set(ALL_PROGS_LIST
  BRAINSMultiModeSegment
  )
foreach(prog ${ALL_PROGS_LIST})
  StandardBRAINSBuildMacro(NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()

if(BUILD_TESTING AND NOT Slicer_BUILD_BRAINSTOOLS)
  #  add_subdirectory(TestSuite)
endif()
