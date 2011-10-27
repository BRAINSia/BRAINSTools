
# --------------------------------------------------------- #
# BRAINSCut Testing
# --------------------------------------------------------- #
##
#if ( 1 ) ## HACK:  Just silencing failing test until they can be properly addressed

## BRAINSCut NetConfiguration(XML) file for testing
##
set( InputNetConfigurationForTesting ${BRAINSCut_SOURCE_DIR}/TestSuite/NetConfigurations/BRAINSCutTest.xml.in)
set( OutputNetConfigurationForTesting ${CMAKE_CURRENT_BINARY_DIR}/BRAINSCutTest.xml )

configure_file( ${InputNetConfigurationForTesting} ${OutputNetConfigurationForTesting}  IMMEDIATE @ONLY)

##
## BRAINSCut Test Run
add_test( GenerateProbabilityTest $<TARGET_FILE:BRAINSCut>
                                  --generateProbability --netConfiguration ${OutputNetConfigurationForTesting})
add_test( CreateVectorsTest $<TARGET_FILE:BRAINSCut>
                                  --createVectors --netConfiguration ${OutputNetConfigurationForTesting})
add_test( TrainModelTest $<TARGET_FILE:BRAINSCut>
                                  --trainModel --netConfiguration ${OutputNetConfigurationForTesting})
add_test( ApplyModelTest $<TARGET_FILE:BRAINSCut>
                                  --applyModel --netConfiguration ${OutputNetConfigurationForTesting})
add_test( AllInOneTest $<TARGET_FILE:BRAINSCut>
                                  --generateProbability  --createVectors  --trainModel  --applyModel
                                  --netConfiguration ${OutputNetConfigurationForTesting})
#endif()







