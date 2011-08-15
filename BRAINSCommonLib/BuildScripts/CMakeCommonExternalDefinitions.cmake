
#-----------------------------------------------------------------------------
## The FORCE_BUILD_CHECK macro adds a forecebuild step that will cause the
## external project build process to be checked for updates each time
## a dependent project is built.  It MUST be called AFTER the ExternalProject_Add
## step for the project that you want to force building on.
macro(FORCE_BUILD_CHECK  proj)
    ExternalProject_Add_Step(${proj} forcebuild
      COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BUILD_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
      DEPENDEES configure
      DEPENDERS build
      ALWAYS 1
    )
endmacro()

