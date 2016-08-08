## NOTE: On mavericks you need to edit mexopts.sh to match required compiler options.
##  vim /Applications/Matlab_R2013b.app/bin/mexopts.sh
find_package(Matlab COMPONENTS MX_LIBRARY ENG_LIBRARY REQUIRED)
include_directories(${Matlab_INCLUDE_DIRS})

find_package(VTK 7 REQUIRED)
include(${VTK_USE_FILE})

find_package(ITK 4.5 REQUIRED)
include(${ITK_USE_FILE})

#####################
#####################
# Provide a macro to build the mex files from
# within CMake
#####################
# BuildMex.cmake
# \author Kent Williams norman-k-williams@uiowa.edu
# \author Hans J. Johnson hans-johnson@uiowa.edu
include(CMakeParseArguments)
include_directories(${Matlab_INCLUDE_DIRS})
#
# BuildMex -- arguments
# MEXNAME = root of mex library name
# TARGETDIR = location for the mex library files to be created
# SOURCE = list of source files
# LIBRARIES = libraries needed to link mex library
macro(BuildMex)
  set(oneValueArgs MEXNAME TARGETDIR)
  set(multiValueArgs SOURCE LIBRARIES)
  cmake_parse_arguments(BuildMex "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # message("MEXNAME=${BuildMex_MEXNAME} SOURCE=${BuildMex_SOURCE} LIBRARIES=${BuildMex_LIBRARIES}")
#--  set_source_files_properties(${BuildMex_SOURCE} COMPILE_DEFINITIONS -DMatlab_MEX_FILE )
  matlab_add_mex(
      NAME ${BuildMex_MEXNAME}
      SRC ${BuildMex_SOURCE}
#--      [OUTPUT_NAME output_name]
#--      [DOCUMENTATION file.txt]
      LINK_TO ${BuildMex_LIBRARIES}
#--      [...]
  )
endmacro(BuildMex)

matlab_add_mex(
      NAME itkLoadWithMetadata
      SRC itkLoadWithMetadata.cxx nrrdCommon.cxx
#--      [DOCUMENTATION file.txt]
      LINK_TO ${ITK_LIBRARIES} ${VTK_LIBRARIES}
)

matlab_add_mex(
      NAME itkSaveWithMetadata
      SRC itkSaveWithMetadata.cxx nrrdCommon.cxx
#--      [DOCUMENTATION file.txt]
      LINK_TO ${ITK_LIBRARIES} ${VTK_LIBRARIES}
)

matlab_add_mex(
      NAME vtkLoadPolyData
      SRC vtkLoadPolyData.cxx
#--      [DOCUMENTATION file.txt]
      LINK_TO ${VTK_LIBRARIES}
)

matlab_add_mex(
      NAME vtkSavePolyData
      SRC vtkSavePolyData.cxx
#--      [DOCUMENTATION file.txt]
      LINK_TO ${VTK_LIBRARIES}
)


if(0)
BuildMex(MEXNAME itkLoadWithMetadata
  TARGETDIR "${${PRIMARY_PROJECT_NAME}_SOURCE_DIR}/MATLAB_SCRIPTS"
  SOURCE itkLoadWithMetadata.cxx nrrdCommon.cxx
  LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})
BuildMex(MEXNAME itkSaveWithMetadata
  TARGETDIR "${${PRIMARY_PROJECT_NAME}_SOURCE_DIR}/MATLAB_SCRIPTS"
  SOURCE itkSaveWithMetadata.cxx nrrdCommon.cxx
  LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES})

BuildMex(MEXNAME vtkLoadPolyData
  TARGETDIR "${${PRIMARY_PROJECT_NAME}_SOURCE_DIR}/MATLAB_SCRIPTS"
  SOURCE vtkLoadPolyData.cxx
  LIBRARIES ${VTK_LIBRARIES} )

BuildMex(MEXNAME vtkSavePolyData
  TARGETDIR "${${PRIMARY_PROJECT_NAME}_SOURCE_DIR}/MATLAB_SCRIPTS"
  SOURCE vtkSavePolyData.cxx
  LIBRARIES ${VTK_LIBRARIES} )

endif()
