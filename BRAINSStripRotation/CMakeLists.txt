
project(BRAINSStripRotation)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
find_package(ITK COMPONENTS
  ITKIOImageBase
  ITKIOTransformBase
  ITKTransform
  ${ITK_IO_MODULES_USED}
  REQUIRED
)
if(Slicer_BUILD_BRAINSTOOLS)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1)
endif()
include(${ITK_USE_FILE})

set(prog BRAINSStripRotation)

StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
