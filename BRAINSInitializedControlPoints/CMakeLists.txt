project(BRAINSInitializedControlPoints)
set(LOCAL_PROJECT_NAME BRAINSInitializedControlPoints)

#-----------------------------------------------------------------------------
# Output directories.
#

set(CLP BRAINSInitializedControlPoints)
set(${CLP}_SOURCE BRAINSInitializedControlPoints.cxx)
CONFIGUREBRAINSORSLICERPROPERTIES(${CLP} BRAINSInitializedControlPoints.xml "${${CLP}_SOURCE}" "" "main" "")
target_link_libraries(BRAINSInitializedControlPoints ${ITK_LIBRARIES} crypto)

