#
# Follow the boost suggestions
#

if(WIN32)

  execute_process(COMMAND ./b2 install --prefix=${BOOST_INSTALL_DIR}
    --without-atomic --without-chrono --without-context --without-date_time --without-exception --without-filesystem --without-graph --without-graph_parallel --without-iostreams --without-locale --without-math --without-mpi --without-program_options --without-python --without-random --without-regex --without-serialization --without-signals
    --without-test --without-timer --without-wave
# Needed for UKF  --without-system --without-thread
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE build_result)

else(WIN32)
if(CMAKE_COMPILER_IS_CLANGXX)
  set(CLANG_ARGS "toolset=clang cxxflags=\"-std=c++11 -stdlib=libc++\" linkflags=\"-stdlib=libc++\"")
endif()
  execute_process(COMMAND ./b2 ${CLANG_ARGS} install
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE build_result)

endif(WIN32)

return(${build_result})
