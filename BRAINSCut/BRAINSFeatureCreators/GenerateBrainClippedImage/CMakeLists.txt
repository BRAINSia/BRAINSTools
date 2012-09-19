##- project(GenerateBrainClippedImage)

## Build all the programs
##
set(ALL_PROGS_LIST
  GenerateBrainClippedImage
  )
foreach( prog ${ALL_PROGS_LIST} )
  StandardBRAINSBuildMacro( NAME ${prog} TARGET_LIBRARIES BRAINSCommonLib )
endforeach()
