#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKCommon
  ITKIOImageBase
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

StandardBRAINSBuildMacro(NAME BRAINSDWICleanup TARGET_LIBRARIES ${ITK_LIBRARIES})
