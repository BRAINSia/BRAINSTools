##- project(GenerateTestImage)

## Build all the programs
##
set(ALL_PROGS_LIST
  GenerateTestImage
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES ${ITK_LIBRARIES} )
endforeach()
