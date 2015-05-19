#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#

FindITKUtil(BRAINSDWICleanup_ITK
)

StandardBRAINSBuildMacro(NAME BRAINSDWICleanup TARGET_LIBRARIES ${BRAINSDWICleanup_ITK_LIBRARIES} BRAINSCommonLib)
