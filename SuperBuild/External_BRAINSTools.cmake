## The purpose of this file is to provide a common way
## to build all the BRAINSTools for distribution in
## various build environments.


## This macro builds each of the BRAINSTools Components.
macro(BuildExtPackage PackageName PACKAGE_DEPENDANCIES)
  ## The following assumes that each PackageName
  ## is a git submodule of the current ${CMAKE_CURRENT_SOURCE_DIR}
  ## that was created with a command similar to:
  ## cd ${CMAKE_CURRENT_SOURCE_DIR} 
  ## git submodule add -b master -- git://github.com/BRAINSia/BRAINSABC.git
  ## git submodule init
  ## git submodule update
  set(SRCCMDS SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${PackageName} )
  if(0) ## TODO:  if PACKAGE_DEPENDANCIES contains VTK, then add -DVTK_DIR:
    set(VTK_BUILD_FLAGS
                 -DVTK_DIR:PATH=${VTK_DIR}
                 )
  elseif()
    set(VTK_BUILD_FLAGS "")
  endif()
  if(0) ## TODO:  if PACKAGE_DEPENDANCIES contains OpenCV, then add -DVTK_DIR:
    set(OpenCV_BUILD_FLAGS
                 -DOpenCV_DIR:PATH=${OpenCV_DIR}
                 )
  elseif()
    set(OpenCV_BUILD_FLAGS "")
  endif()
  if(0) ## TODO:  if PACKAGE_DEPENDANCIES contains Qt, then add -DVTK_DIR:
    set(QT_BUILD_FLAGS
                -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
                -D${CMAKE_PROJECT_NAME}_USE_QT:BOOL=${${CMAKE_PROJECT_NAME}_USE_QT}
                )
  elseif()
    set(QT_BUILD_FLAGS "")
  endif()

  if(NOT INTEGRATE_WITH_SLICER)
    set(INTEGRATE_WITH_SLICER OFF)
  endif()

  message(STATUS "                                    CONFIGUREING PROJECT ${PackageName} WITH : ${SlicerExecutionModel_DIR}")
  message(STATUS " -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}")
## Configure and build the package.
  ExternalProject_Add(${PackageName}
    ${SRCCMDS}
    BINARY_DIR ${PackageName}-build
    CMAKE_GENERATOR ${gen}
    DEPENDS ${PACKAGE_DEPENDANCIES}
    CMAKE_ARGS
    ${ep_common_args}
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
    -DITK_DIR:PATH=${ITK_DIR}
    -DINTEGRATE_WITH_SLICER:BOOL=${INTEGRATE_WITH_SLICER}
    -DSlicer_SOURCE_DIR:PATH=${Slicer_SOURCE_DIR}
    -DBRAINSCommonLib_DIR:PATH=${BRAINSCommonLib_DIR}
    -DBUILD_TESTING:BOOL=ON
    -D${CMAKE_PROJECT_NAME}_USE_ITK4:BOOL=ON
    ${VTK_BUILD_FLAGS}
    ${OpenCV_BUILD_FLAGS}
    ${QT_BUILD_FLAGS}
    #INSTALL_COMMAND ""
    #INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )

## Force building the package every time (i.e. look for code changes).i
## Without this step, after the first successful build, the
## the package is never built again.  This means that source code
## changes do not trigger rebuilds
  ExternalProject_Add_Step(${PackageName} forcebuild
    COMMAND ${CMAKE_COMMAND} -E remove
    ${CMAKE_CURRENT_BUILD_DIR}/${PackageName}-prefix/src/${PackageName}-stamp/${PackageName}-build
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
    )

## Setup some common variables for each package.
  set(${PackageName}_DEPEND "${proj}")
  set(${PackageName}_DIR ${CMAKE_INSTALL_PREFIX}/lib/${PackageName})
  set(${PackageName}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${PackageName} )
  message(STATUS "${CMAKE_PROJECT_NAME}")
  message(STATUS "${PackageName}_DIR = ${CMAKE_INSTALL_PREFIX}/lib/${PackageName})
  message(STATUS "${PackageName}_DIR = ${${PackageName}_DIR}})
  list(APPEND BRAINSTools_PROVIDES ${PackageName})
endmacro(BuildExtPackage)


option(USE_BRAINSCommonLib             "Build BRAINSCommonLib"             ON)
option(USE_BRAINSFit                   "Build BRAINSFit"                   OFF)
#option(USE_BRAINSCut                   "Build BRAINSCut"                   OFF)
option(USE_BRAINSABC                   "Build BRAINSABC"                   OFF)
option(USE_BRAINSConstellationDetector "Build BRAINSConstellationDetector" OFF)

#-----------------------------------------------------------------------------
# BRAINSCommonLib
#-----------------------------------------------------------------------------
if(USE_BRAINSCommonLib)
  ## This is required option(USE_BRAINSCommonLib "Build BRAINSCommonLib" ON)
  BuildExtPackage(BRAINSCommonLib "SlicerExecutionModel")
endif()

#-----------------------------------------------------------------------------
# BRAINSFit
#-----------------------------------------------------------------------------
if(USE_BRAINSFit)
  BuildExtPackage(BRAINSFit "BRAINSCommonLib" )
endif()

#-----------------------------------------------------------------------------
# BRAINSConstellationDetector
#-----------------------------------------------------------------------------
if(USE_BRAINSConstellationDetector)
  BuildExtPackage(BRAINSConstellationDetector "BRAINSCommonLib" )
endif()

#-----------------------------------------------------------------------------
# ReferenceAtlas
#-----------------------------------------------------------------------------
if(USE_BRAINSABC) # OR USE_BRAINSCut)
  # Define the atlas subdirectory in one place
  set(${CMAKE_PROJECT_NAME}_RUNTIME_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/bin)
  include(External_ReferenceAtlas)
  list(APPEND ${CMAKE_PROJECT_NAME}_DEPENDENCIES ${ReferenceAtlas_DEPEND})
endif()

#-----------------------------------------------------------------------------
# BRAINSABC
#-----------------------------------------------------------------------------
if(USE_BRAINSABC)
  BuildExtPackage(BRAINSABC "BRAINSCommonLib;${ReferenceAtlas_DEPEND}" )
endif()


if(0)
  BuildExtPackage(BRAINSResample "BRAINSCommonLib" )
  BuildExtPackage(BRAINSROIAuto "BRAINSCommonLib" )
  BuildExtPackage(GTRACT "BRAINSCommonLib" )
  BuildExtPackage(BRAINSCut "BRAINSCommonLib;${OpenCV_DEPEND}" )
  BuildExtPackage(BRAINSMush "BRAINSCommonLib" )
  BuildExtPackage(BRAINSDemonWarp "BRAINSCommonLib" )
  BuildExtPackage(BRAINSMultiModeSegment "BRAINSCommonLib" )
  BuildExtPackage(BRAINSInitializedControlPoints "BRAINSCommonLib" )
endif() ## Comment out things that are not yet ready


