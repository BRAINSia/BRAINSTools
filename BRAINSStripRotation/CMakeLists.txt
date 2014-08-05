
project(BRAINSStripRotation)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil(ITKTransform
  ITKImageCompare)

set(prog BRAINSStripRotation)

StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
