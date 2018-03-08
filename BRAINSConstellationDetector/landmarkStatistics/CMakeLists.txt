##- project(landmarkStatistics)

## Build all the programs
##
## Not CLP programs
set(PROGS_LIST_1
  LmkStatistics
  )
foreach(prog ${PROGS_LIST_1})
  add_executable(${prog} ${prog}.cxx)
  target_link_libraries(${prog} BRAINSCommonLib ${BRAINSConstellationDetector_ITK_LIBRARIES} ${VTK_LIBRARIES})
  set_target_properties(${prog} PROPERTIES FOLDER ${MODULE_FOLDER})
endforeach()

## Simple CLP programs
set(PROGS_LIST_2
  GenerateAverageLmkFile
  )

foreach(prog ${PROGS_LIST_2})
    StandardBRAINSBuildMacro(NAME ${prog}
     TARGET_LIBRARIES landmarksConstellationCOMMONLIB BRAINSCommonLib)
endforeach()
