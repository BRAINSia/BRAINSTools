# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

option(USE_SYSTEM_UNCRUSTIFY  "Use already-installed uncrustify" OFF)
set(UNCRUSTIFY_EXE "" CACHE PATH "absolute path for uncrustify utility")

if(USE_SYSTEM_UNCRUSTIFY)
  find_program(UNCRUSTIFY_EXE uncrustify
    DOC "path of uncrustify program"
    )
  if("${UNCRUSTIFY_EXE}" STREQUAL "")
    message(WARNING "To use the system uncrustify, set UNCRUSTIFY_EXE")
  endif()

endif()

if(NOT USE_SYSTEM_UNCRUSITY)
  ExternalProject_add(uncrustify
    GIT_REPOSITORY git://uncrustify.git.sourceforge.net/gitroot/uncrustify/uncrustify
    GIT_TAG 60f3681da60462eda539b78e0c6c3eea823481e5
    CONFIGURE_COMMAND ../uncrustify/configure
    --prefix=${CMAKE_BINARY_DIR}/Utils
    )
  set(UNCRUSTIFY_EXE "${CMAKE_BINARY_DIR}/Utils/bin/uncrustify" CACHE PATH
    "Absolute path for uncrustify")
endif()

option(USE_SYSTEM_KWSTYLE  "Use already-installed KWStyle" OFF)
set(KWSTYLE_EXE "" CACHE PATH "absolute path for KWStyle utility")

if(USE_SYSTEM_KWSTYLE)
  find_program(KWSTYLE_EXE KWStyle
    DOC "path of KWStyle program"
    )
  if("${KWSTYLE_EXE}" STREQUAL "")
    message(WARNING "To use the system KWStyle, set KWSTYLE_EXE")
  endif()

endif()

if(NOT USE_SYSTEM_UNCRUSITY)
  ExternalProject_add(KWStyle
    CVS_REPOSITORY :pserver:anoncvs@public.kitware.com:/cvsroot/KWStyle
    CVS_MODULE KWStyle
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Utils
    -DCMAKE_BUILD_TYPE:STRING=Release
    )
  set(KWSTYLE_EXE "${CMAKE_BINARY_DIR}/Utils/bin/KWStyle" CACHE PATH
    "Absolute path for KWStyle")
endif()

option(USE_SYSTEM_CPPCHECK  "Use already-installed cppcheck" OFF)
set(CPPCHECK_EXE "" CACHE PATH "absolute path for cppcheck utility")

if(USE_SYSTEM_CPPCHECK)
  find_program(CPPCHECK_EXE cppcheck
    DOC "path of cppcheck program"
    )
  if("${CPPCHECK_EXE}" STREQUAL "")
    message(WARNING "To use the system cppcheck, set CPPCHECK_EXE")
  endif()

endif()

if(NOT USE_SYSTEM_CPPCHECK)
  ExternalProject_add(cppcheck
    GIT_REPOSITORY git://github.com/danmar/cppcheck.git
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} DESTDIR=${CMAKE_BINARY_DIR} HAVE_RULES=no
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} DESTDIR=${CMAKE_BINARY_DIR} HAVE_RULES=no
    PREFIX=Utils
    )
  set(CPPCHECK_EXE "${CMAKE_BINARY_DIR}/Utils/bin/cppcheck" CACHE PATH
    "Absolute path for cppcheck")
  option(USE_SYSTEM_UNCRUSTIFY  "Use already-installed uncrustify" OFF)
  set(UNCRUSTIFY_EXE "" CACHE PATH "absolute path for uncrustify utility")

endif()
