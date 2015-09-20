#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSDWICleanup_ITK
ITKNrrdIO
)

StandardBRAINSBuildMacro(NAME BRAINSDWICleanup TARGET_LIBRARIES ${BRAINSDWICleanup_ITK_LIBRARIES} BRAINSCommonLib)
