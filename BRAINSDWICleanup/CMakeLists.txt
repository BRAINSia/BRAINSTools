#-----------------------------------------------------------------------------
# Dependencies.
#

#
# ITK
#
FindITKUtil()

StandardBRAINSBuildMacro(NAME BRAINSDWICleanup TARGET_LIBRARIES ${ITK_LIBRARIES})
