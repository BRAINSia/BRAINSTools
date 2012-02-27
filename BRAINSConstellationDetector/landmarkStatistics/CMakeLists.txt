project(landmarkStatistics)

## Build all the programs
##
set(ALL_PROGS_LIST
  LmkStatistics
  GenerateAverageLmkFile
  )
foreach(prog ${ALL_PROGS_LIST})
  add_executable(${prog} ${prog}.cxx)
  target_link_libraries(${prog} BRAINSCommonLib ${ITK_LIBRARIES})
endforeach()
