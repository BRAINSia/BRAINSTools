##- project(GenerateSummedGradientImage)

## Build all the programs
##
set(ALL_PROGS_LIST
  GenerateSummedGradientImage
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()
