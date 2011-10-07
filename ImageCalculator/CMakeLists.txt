PROJECT(ImageCalculator)

ADD_EXECUTABLE(ImageCalculator
  ImageCalculator.cxx
  ImageCalculatorProcess2D.cxx
  ImageCalculatorProcess3D.cxx
  ImageCalculatorUtils.cxx)
TARGET_LINK_LIBRARIES(ImageCalculator ${ITK_LIBRARIES} )

ADD_EXECUTABLE(ImageGenerate ImageGenerate.cxx )
TARGET_LINK_LIBRARIES(ImageGenerate ${ITK_LIBRARIES})

IF( ____BUILD_TESTING )

  SET(IMAGE_CALCULATOR_TESTS ${CXX_TEST_PATH}/ImageCalculatorTests)
  SET(IC_TSTDATA ${ImageCalculator_SOURCE_DIR}/TestImages)
  SET(IC_STATCMDS --statNPX --statMIN --statAVG --statMAX --statVAR --statSUM)
  SET(IC_BIN "${ImageCalculator_BINARY_DIR}")
  SET(IC_OUT --out "${IC_BIN}/Test.png")
  #Test to show arithmetic operation and printing the statistics.
  ADD_TEST(ImageCalculator2DTest1 ${IMAGE_CALCULATOR_TESTS}
    --compare
       "${IC_BIN}/Test.png"
       "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
       --in "${IC_TSTDATA}/UpperHalfHundreds.png"
           "${IC_TSTDATA}/LowerHalfHundreds.png"
       --intype UCHAR --outtype UCHAR --add
       ${IC_OUT}
       -d 2 ${IC_STATCMDS}  )

  #Test gaussian filtering on input and output
  ADD_TEST(ImageCalculator2DGaussianTest ${IMAGE_CALCULATOR_TESTS}
    --compare
       "${IC_BIN}/Test.png"
       "${IC_TSTDATA}/GaussianIOFilteringTest.png"
    ImageCalculatorTest
       --in "${IC_TSTDATA}/UpperHalfHundreds.png"
       --intype UCHAR --outtype UCHAR --add --ifgaussiansigma 1 --ofgaussiansigma 1
       ${IC_OUT}
       -d 2 ${IC_STATCMDS}  )
  #Test for more than 2 images different data types for Input image and output image
  #
  # Test add
  ADD_TEST(ImageCalculator2DTest2 ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllFiveHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
        --intype UCHAR --outtype USHORT --add
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )

  #
  # Test SUB
  ADD_TEST(ImageCalculator2DTest3 ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllZeroes.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
        --intype UCHAR --outtype USHORT --sub
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )
  #
  # Test Mult
  ADD_TEST(ImageCalc2DMult ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllTens.png"
            "${IC_TSTDATA}/AllTens.png"
        --intype UCHAR --outtype UCHAR --mul
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )
  #
  # Test Mult
  ADD_TEST(ImageCalc2DDiv ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllTens.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllTens.png"
        --intype UCHAR --outtype UCHAR --div
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )

  #
  # Test Avg
  ADD_TEST(ImageCalc2DAvg ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllFifties.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/LowerHalfHundreds.png"
            "${IC_TSTDATA}/UpperHalfHundreds.png"
        --intype UCHAR --outtype UCHAR --avg
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )


  #
  # Test Input Square
  ADD_TEST(ImageCalc2DTestSquare ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllTens.png"
            "${IC_TSTDATA}/AllZeroes.png"
        --ifsqr --intype UCHAR --outtype UCHAR --sub
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )
  #
  # Test Input Square Root
  ADD_TEST(ImageCalc2DTestSquareRoot ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllTens.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllZeroes.png"
        --ifsqrt --intype UCHAR --outtype USHORT --sub
        ${IC_OUT}
        -d 2 ${IC_STATCMDS}  )

  # Test for 3d Images and division
  #
  # test Divide
  ADD_TEST(ImageCalculator3DTest1 ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.hdr"
        "${IC_TSTDATA}/AllTwos.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
            "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR --div
        --out "${IC_BIN}/Test.hdr"
        -d 3 ${IC_STATCMDS}  )

  #A simple scalar multiplcation test
  # test multiply
  ADD_TEST(ImageCalc3DTestOMult ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.hdr"
        "${IC_TSTDATA}/TwoHundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/Test.hdr"
        -d 3 ${IC_STATCMDS} --ofmulc 2  )

  #A simple scalar multiplcation test
  # test multiply
  ADD_TEST(ImageCalc3DTestIMult ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.hdr"
        "${IC_TSTDATA}/TwoHundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/Test.hdr"
        -d 3 ${IC_STATCMDS} --ifmulc 2  )

  #A simple scalar division test
  # test
  ADD_TEST(ImageCalc3DTestODiv ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.hdr"
        "${IC_TSTDATA}/Hundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/Test.hdr"
        -d 3 ${IC_STATCMDS} --ofdivc 2  )

  #A simple scalar divide test
  # test diviply
  ADD_TEST(ImageCalc3DTestIDiv ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.hdr"
        "${IC_TSTDATA}/Hundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/Test.hdr"
        -d 3 ${IC_STATCMDS} --ifdivc 2  )

  #A simple scalar addition test
  # test
  ADD_TEST(ImageCalc2DTestIAdd ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllFifties.png"
        --intype UCHAR --outtype UCHAR
        ${IC_OUT}
        -d 2 ${IC_STATCMDS} --ifaddc 50 )
  ADD_TEST(ImageCalc2DTestOAdd ${IMAGE_CALCULATOR_TESTS}
    --compare
        "${IC_BIN}/Test.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllFifties.png"
        --intype UCHAR --outtype UCHAR
        ${IC_OUT}
        -d 2 ${IC_STATCMDS} --ofaddc 50 )


  ADD_EXECUTABLE(ImageCalculatorTests ImageCalculatorTests.cxx ImageCalculatorUtils.cxx ImageCalculatorProcess2D.cxx ImageCalculatorProcess3D.cxx)
  TARGET_LINK_LIBRARIES(ImageCalculatorTests ${ITK_LIBRARIES})

ENDIF( BUILD_TESTING )

