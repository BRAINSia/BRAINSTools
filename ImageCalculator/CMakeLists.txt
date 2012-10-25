##- project(ImageCalculator)

add_executable(ImageCalculator
  ImageCalculator.cxx
  ImageCalculatorProcess2D.cxx
  ImageCalculatorProcess3D.cxx
  ImageCalculatorUtils.cxx)
target_link_libraries(ImageCalculator ${ITK_LIBRARIES} )

add_executable(ImageGenerate ImageGenerate.cxx )
target_link_libraries(ImageGenerate ${ITK_LIBRARIES})

if(BUILD_TESTING)

  add_executable(ImageCalculatorTests
    ImageCalculatorTests.cxx ImageCalculatorUtils.cxx
    ImageCalculatorProcess2D.cxx
    ImageCalculatorProcess3D.cxx)
  target_link_libraries(ImageCalculatorTests ${ITK_LIBRARIES})
  set_target_properties(ImageCalculatorTests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)

  set(IC_TSTDATA ${CMAKE_CURRENT_SOURCE_DIR}/TestImages)
  set(IC_STATCMDS --statNPX --statMIN --statAVG --statMAX --statVAR --statSUM)
  set(IC_BIN "${CMAKE_CURRENT_BINARY_DIR}/testbin")
  set(IC_OUT --out "${IC_BIN}/Test.png")
  #Test to show arithmetic operation and printing the statistics.
  add_test(NAME ImageCalculator2DTest1 COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
       "${IC_BIN}/ImageCalculator2DTest1.png"
       "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
       --in "${IC_TSTDATA}/UpperHalfHundreds.png"
           "${IC_TSTDATA}/LowerHalfHundreds.png"
       --intype UCHAR --outtype UCHAR --add
       --out ${IC_BIN}/ImageCalculator2DTest1.png
       -d 2 ${IC_STATCMDS}  )

  #Test gaussian filtering on input and output
  add_test(NAME ImageCalculator2DGaussianTest COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
       "${IC_BIN}/ImageCalculator2DGaussianTest.png"
       "${IC_TSTDATA}/GaussianIOFilteringTest.png"
    ImageCalculatorTest
       --in "${IC_TSTDATA}/UpperHalfHundreds.png"
       --intype UCHAR --outtype UCHAR --add --ifgaussiansigma 1 --ofgaussiansigma 1
       --out ${IC_BIN}/ImageCalculator2DGaussianTest.png
       -d 2 ${IC_STATCMDS}  )
  
  #Test for more than 2 images different data types for Input image and output image
  #
  # Test add
  add_test(NAME ImageCalculator2DTest2 COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalculator2DTest2.png"
        "${IC_TSTDATA}/AllFiveHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
        --intype UCHAR --outtype USHORT --add
        --out ${IC_BIN}/ImageCalculator2DTest2.png
        -d 2 ${IC_STATCMDS}  )

  #
  # Test SUB
  add_test(NAME ImageCalculator2DTest3 COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalculator2DTest3.png"
        "${IC_TSTDATA}/AllZeroes.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
        --intype UCHAR --outtype USHORT --sub
        --out ${IC_BIN}/ImageCalculator2DTest3.png
        -d 2 ${IC_STATCMDS}  )
  #  
  # Test histogram equalization on inputs before making the difference image
  add_test(NAME ImageCalculator2DHisteqTest COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
       "${IC_BIN}/ImageCalculator2DHisteqTest.png"
       "${IC_TSTDATA}/Histeq2DTest.png"
    ImageCalculatorTest
       --in "${IC_TSTDATA}/AllFiveHundreds.png"
            "${IC_TSTDATA}/AllHundreds.png"
       --intype UCHAR --outtype UCHAR --sub --ifhisteq 32
       --out ${IC_BIN}/ImageCalculator2DHisteqTest.png
       -d 2 ${IC_STATCMDS}  )
      
  #
  # Test Mult
  add_test(NAME ImageCalc2DMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DMult.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllTens.png"
            "${IC_TSTDATA}/AllTens.png"
        --intype UCHAR --outtype UCHAR --mul
        --out ${IC_BIN}/ImageCalc2DMult.png
        -d 2 ${IC_STATCMDS}  )
  #
  # Test Mult
  add_test(NAME ImageCalc2DDiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DDiv.png"
        "${IC_TSTDATA}/AllTens.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllTens.png"
        --intype UCHAR --outtype UCHAR --div
        --out ${IC_BIN}/ImageCalc2DDiv.png
        -d 2 ${IC_STATCMDS}  )

  #
  # Test Avg
  add_test(NAME ImageCalc2DAvg COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DAvg.png"
        "${IC_TSTDATA}/AllFifties.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/LowerHalfHundreds.png"
            "${IC_TSTDATA}/UpperHalfHundreds.png"
        --intype UCHAR --outtype UCHAR --avg
        --out ${IC_BIN}/ImageCalc2DAvg.png
        -d 2 ${IC_STATCMDS}  )


  #
  # Test Input Square
  add_test(NAME ImageCalc2DTestSquare COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DTestSquare.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllTens.png"
            "${IC_TSTDATA}/AllZeroes.png"
        --ifsqr --intype UCHAR --outtype UCHAR --sub
        --out ${IC_BIN}/ImageCalc2DTestSquare.png
        -d 2 ${IC_STATCMDS}  )
  #
  # Test Input Square Root
  add_test(NAME ImageCalc2DTestSquareRoot COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DTestSquareRoot.png"
        "${IC_TSTDATA}/AllTens.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllHundreds.png"
            "${IC_TSTDATA}/AllZeroes.png"
        --ifsqrt --intype UCHAR --outtype USHORT --sub
        --out ${IC_BIN}/ImageCalc2DTestSquareRoot.png
        -d 2 ${IC_STATCMDS}  )

  # Test for 3d Images and division
  #
  # test Divide
  add_test(NAME ImageCalculator3DTest1 COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalculator3DTest1.hdr"
        "${IC_TSTDATA}/AllTwos.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
            "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR --div
        --out "${IC_BIN}/ImageCalculator3DTest1.hdr"
        -d 3 ${IC_STATCMDS}  )

  #A simple scalar multiplcation test
  # test multiply
  add_test(NAME ImageCalc3DTestOMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc3DTestOMult.hdr"
        "${IC_TSTDATA}/TwoHundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/ImageCalc3DTestOMult.hdr"
        -d 3 ${IC_STATCMDS} --ofmulc 2  )

  #A simple scalar multiplcation test
  # test multiply
  add_test(NAME ImageCalc3DTestIMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc3DTestIMult.hdr"
        "${IC_TSTDATA}/TwoHundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/Hundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/ImageCalc3DTestIMult.hdr"
        -d 3 ${IC_STATCMDS} --ifmulc 2  )

  #A simple scalar division test
  # test
  add_test(NAME ImageCalc3DTestODiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc3DTestODiv.hdr"
        "${IC_TSTDATA}/Hundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/ImageCalc3DTestODiv.hdr"
        -d 3 ${IC_STATCMDS} --ofdivc 2  )

  #A simple scalar divide test
  # test diviply
  add_test(NAME ImageCalc3DTestIDiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc3DTestIDiv.hdr"
        "${IC_TSTDATA}/Hundred.hdr"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/TwoHundred.hdr"
        --intype UCHAR --outtype UCHAR
        --out "${IC_BIN}/ImageCalc3DTestIDiv.hdr"
        -d 3 ${IC_STATCMDS} --ifdivc 2  )

  #A simple scalar addition test
  # test
  add_test(NAME ImageCalc2DTestIAdd COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        ${IC_BIN}/ImageCalculator3DTest1.png
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllFifties.png"
        --intype UCHAR --outtype UCHAR
        --out ${IC_BIN}/ImageCalculator3DTest1.png
        -d 2 ${IC_STATCMDS} --ifaddc 50 )
  add_test(NAME ImageCalc2DTestOAdd COMMAND $<TARGET_FILE:ImageCalculatorTests>
    --compare
        "${IC_BIN}/ImageCalc2DTestOAdd.png"
        "${IC_TSTDATA}/AllHundreds.png"
    ImageCalculatorTest
        --in "${IC_TSTDATA}/AllFifties.png"
        --intype UCHAR --outtype UCHAR
        --out ${IC_BIN}/ImageCalc2DTestOAdd.png
        -d 2 ${IC_STATCMDS} --ofaddc 50 )


endif()

