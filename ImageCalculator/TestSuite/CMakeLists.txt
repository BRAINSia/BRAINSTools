#########
add_executable(ImageCalculatorTests
  ../ImageCalculatorTests.cxx ../ImageCalculatorUtils.cxx
  ../ImageCalculatorProcess2D.cxx
  ../ImageCalculatorProcess3D.cxx)
target_link_libraries(ImageCalculatorTests ${ITK_LIBRARIES})
set_target_properties(ImageCalculatorTests PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testbin)

set(IC_STATCMDS --statNPX --statMIN --statAVG --statMAX --statVAR --statSUM)
set(IC_BIN "${CMAKE_CURRENT_BINARY_DIR}/testbin")
set(IC_OUT --out "${IC_BIN}/Test.png")
#Test to show arithmetic operation and printing the statistics.
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator2DTest1 COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
     "${IC_BIN}/ImageCalculator2DTest1.png"
     "DATA{${TestData_DIR}/AllHundreds.png}"
  ImageCalculatorTest
     --in "DATA{${TestData_DIR}/UpperHalfHundreds.png}"
         "DATA{${TestData_DIR}/LowerHalfHundreds.png}"
     --intype UCHAR --outtype UCHAR --add
     --out ${IC_BIN}/ImageCalculator2DTest1.png
     -d 2 ${IC_STATCMDS}  )

#Test gaussian filtering on input and output
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator2DGaussianTest COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
     "${IC_BIN}/ImageCalculator2DGaussianTest.png"
     "DATA{${TestData_DIR}/GaussianIOFilteringTest.png}"
  ImageCalculatorTest
     --in "DATA{${TestData_DIR}/UpperHalfHundreds.png}"
     --intype UCHAR --outtype UCHAR --add --ifgaussiansigma 1 --ofgaussiansigma 1
     --out ${IC_BIN}/ImageCalculator2DGaussianTest.png
     -d 2 ${IC_STATCMDS}  )

#Test for more than 2 images different data types for Input image and output image
#
# Test add
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator2DTest2 COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalculator2DTest2.png"
      "DATA{${TestData_DIR}/AllFiveHundreds.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllHundreds.png}"
      --intype UCHAR --outtype USHORT --add
      --out ${IC_BIN}/ImageCalculator2DTest2.png
      -d 2 ${IC_STATCMDS}  )

#
# Test SUB
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator2DTest3 COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalculator2DTest3.png"
      "DATA{${TestData_DIR}/AllZeroes.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllHundreds.png}"
      --intype UCHAR --outtype USHORT --sub
      --out ${IC_BIN}/ImageCalculator2DTest3.png
      -d 2 ${IC_STATCMDS}  )
#
# Test histogram equalization on inputs before making the difference image
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator2DHisteqTest COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
     "${IC_BIN}/ImageCalculator2DHisteqTest.png"
     "DATA{${TestData_DIR}/histeq2D.png}"
  ImageCalculatorTest
     --in "DATA{${TestData_DIR}/histeq_input1.png}"
          "DATA{${TestData_DIR}/histeq_input2.png}"
          --intype USHORT --outtype USHORT --sub --ifhisteq 32
     --out ${IC_BIN}/ImageCalculator2DHisteqTest.png
     -d 2 ${IC_STATCMDS}  )

#
# Test Mult
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DMult.png"
      "DATA{${TestData_DIR}/AllHundreds.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllTens.png}"
          "DATA{${TestData_DIR}/AllTens.png}"
      --intype UCHAR --outtype UCHAR --mul
      --out ${IC_BIN}/ImageCalc2DMult.png
      -d 2 ${IC_STATCMDS}  )
#
# Test Mult
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DDiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DDiv.png"
      "DATA{${TestData_DIR}/AllTens.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllTens.png}"
      --intype UCHAR --outtype UCHAR --div
      --out ${IC_BIN}/ImageCalc2DDiv.png
      -d 2 ${IC_STATCMDS}  )

#
# Test Avg
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DAvg COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DAvg.png"
      "DATA{${TestData_DIR}/AllFifties.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/LowerHalfHundreds.png}"
          "DATA{${TestData_DIR}/UpperHalfHundreds.png}"
      --intype UCHAR --outtype UCHAR --avg
      --out ${IC_BIN}/ImageCalc2DAvg.png
      -d 2 ${IC_STATCMDS}  )


#
# Test Input Square
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DTestSquare COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DTestSquare.png"
      "DATA{${TestData_DIR}/AllHundreds.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllTens.png}"
          "DATA{${TestData_DIR}/AllZeroes.png}"
      --ifsqr --intype UCHAR --outtype UCHAR --sub
      --out ${IC_BIN}/ImageCalc2DTestSquare.png
      -d 2 ${IC_STATCMDS}  )
#
# Test Input Square Root
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DTestSquareRoot COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DTestSquareRoot.png"
      "DATA{${TestData_DIR}/AllTens.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllHundreds.png}"
          "DATA{${TestData_DIR}/AllZeroes.png}"
      --ifsqrt --intype UCHAR --outtype USHORT --sub
      --out ${IC_BIN}/ImageCalc2DTestSquareRoot.png
      -d 2 ${IC_STATCMDS}  )

# Test for 3d Images and division
#
# test Divide
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalculator3DTest1 COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalculator3DTest1.hdr"
      "DATA{${TestData_DIR}/AllTwos.hdr,AllTwos.img}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/TwoHundred.hdr,TwoHundred.img}"
          "DATA{${TestData_DIR}/Hundred.hdr,Hundred.img}"
      --intype UCHAR --outtype UCHAR --div
      --out "${IC_BIN}/ImageCalculator3DTest1.hdr"
      -d 3 ${IC_STATCMDS}  )

#A simple scalar multiplcation test
# test multiply
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc3DTestOMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc3DTestOMult.hdr"
      "DATA{${TestData_DIR}/TwoHundred.hdr,TwoHundred.img}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/Hundred.hdr,Hundred.img}"
      --intype UCHAR --outtype UCHAR
      --out "${IC_BIN}/ImageCalc3DTestOMult.hdr"
      -d 3 ${IC_STATCMDS} --ofmulc 2  )

#A simple scalar multiplcation test
# test multiply
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc3DTestIMult COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc3DTestIMult.hdr"
      "DATA{${TestData_DIR}/TwoHundred.hdr,TwoHundred.img}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/Hundred.hdr,Hundred.img}"
      --intype UCHAR --outtype UCHAR
      --out "${IC_BIN}/ImageCalc3DTestIMult.hdr"
      -d 3 ${IC_STATCMDS} --ifmulc 2  )

#A simple scalar division test
# test
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc3DTestODiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc3DTestODiv.hdr"
      "DATA{${TestData_DIR}/Hundred.hdr,Hundred.img}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/TwoHundred.hdr,Hundred.img}"
      --intype UCHAR --outtype UCHAR
      --out "${IC_BIN}/ImageCalc3DTestODiv.hdr"
      -d 3 ${IC_STATCMDS} --ofdivc 2  )

#A simple scalar divide test
# test diviply
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc3DTestIDiv COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc3DTestIDiv.hdr"
      "DATA{${TestData_DIR}/Hundred.hdr,Hundred.img}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/TwoHundred.hdr,TwoHundred.img}"
      --intype UCHAR --outtype UCHAR
      --out "${IC_BIN}/ImageCalc3DTestIDiv.hdr"
      -d 3 ${IC_STATCMDS} --ifdivc 2  )

#A simple scalar addition test
# test
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DTestIAdd COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      ${IC_BIN}/ImageCalculator3DTest1.png
      "DATA{${TestData_DIR}/AllHundreds.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllFifties.png}"
      --intype UCHAR --outtype UCHAR
      --out ${IC_BIN}/ImageCalculator3DTest1.png
      -d 2 ${IC_STATCMDS} --ifaddc 50 )
ExternalData_add_test(${PROJECT_NAME}FetchData NAME ImageCalc2DTestOAdd COMMAND $<TARGET_FILE:ImageCalculatorTests>
  --compare
      "${IC_BIN}/ImageCalc2DTestOAdd.png"
      "DATA{${TestData_DIR}/AllHundreds.png}"
  ImageCalculatorTest
      --in "DATA{${TestData_DIR}/AllFifties.png}"
      --intype UCHAR --outtype UCHAR
      --out ${IC_BIN}/ImageCalc2DTestOAdd.png
        -d 2 ${IC_STATCMDS} --ofaddc 50 )
