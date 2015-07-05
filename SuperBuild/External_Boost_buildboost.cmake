#
# Follow the boost suggestions
#

if(WIN32)
set(XXX "
# Needed for UKF
        --with-system
        --with-thread
# Needed for DTIProcess
        --with-program_options
        --with-locale
        --with-format
        --with-BLAH
# Not Needed
        --without-atomic
        --without-chrono
        --without-context
        --without-date_time
        --without-exception
        --without-filesystem
        --without-graph
        --without-graph_parallel
        --without-iostreams
        --without-log
        --without-math
        --without-mpi
        --without-python
        --without-random
        --without-regex
        --without-serialization
        --without-signals
        --without-test
        --without-timer
        --without-wave
"
)
  execute_process(COMMAND ./b2 install --prefix=${BOOST_INSTALL_DIR}
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE build_result)

else(WIN32)

  execute_process(COMMAND ./b2 install
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE build_result)

endif(WIN32)

return(${build_result})
