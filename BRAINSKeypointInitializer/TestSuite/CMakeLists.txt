add_executable(BRAINSKeypointRecoveryTest BRAINSKeypointRecoveryTest.cxx)
target_link_libraries(BRAINSKeypointRecoveryTest PRIVATE BRAINSKeypointLib)

foreach(estimator RANSAC MSAC PROSAC)
  add_test(NAME BRAINSKeypointRecoveryTest_${estimator}
    COMMAND $<TARGET_FILE:BRAINSKeypointRecoveryTest> ${estimator})
  set_tests_properties(BRAINSKeypointRecoveryTest_${estimator} PROPERTIES
    PROCESSORS 3
    LABELS "BRAINSKeypointInitializer")
endforeach()

add_test(NAME BRAINSKeypointRecoveryTest_Model
  COMMAND $<TARGET_FILE:BRAINSKeypointRecoveryTest> Model)
set_tests_properties(BRAINSKeypointRecoveryTest_Model PROPERTIES
  PROCESSORS 3
  LABELS "BRAINSKeypointInitializer")
