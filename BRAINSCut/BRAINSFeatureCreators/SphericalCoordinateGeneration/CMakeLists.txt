##- project(SphericalCoordinateGeneration)

## Build all the programs
##
set(ALL_PROGS_LIST
  SphericalCoordinateGeneration
  )

include_directories( ${BRAINSTools_SOURCE_DIR}/BRAINSCut )
include_directories( ${BRAINSTools_SOURCE_DIR}/BRAINSCut/BRAINSCutConfiguration )

foreach( prog ${ALL_PROGS_LIST} )
  GENERATECLP( ${prog}.cxx ${prog}.xml )
  # Create the shared library
  add_library( ${prog}Lib SHARED ${prog}.cxx )
  set_target_properties( ${prog}Lib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint" )
  target_link_libraries( ${prog}Lib BRAINSCutCOMMONLIB NetConfigurationCOMMONLIB )
  # Create the executable
  add_executable( ${prog} SphericalCoordinateGeneration.cxx )
  target_link_libraries( ${prog} BRAINSCutCOMMONLIB NetConfigurationCOMMONLIB )
endforeach()
