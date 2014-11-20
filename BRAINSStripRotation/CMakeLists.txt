
project(BRAINSStripRotation)

#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSStripRotation_ITK
  ITKTransform
  ITKImageCompare)

set(prog BRAINSStripRotation)

StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
