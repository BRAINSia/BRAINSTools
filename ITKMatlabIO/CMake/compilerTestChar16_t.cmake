#
# Figure out how the compiler prefers to declare templated functions as friends
# of a templated class.
#
set(TEST_CHAR16_T_TEST_OUTPUT)
try_compile(SUPPORTS_CHAR16_T
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake/compilerTestChar16_t.cxx
      OUTPUT_VARIABLE TEST_CHAR16_T_TEST_OUTPUT
      )


if(SUPPORTS_CHAR16_T)
  if(TEST_CHAR16_T_TEST_OUTPUT MATCHES "[Ww]arning")
    if(NOT TEST_CHAR16_T_TEST_OUTPUT MATCHES " 0 [Ww]arning")
      set(SUPPORTS_CHAR16_T CACHE INTERNAL FALSE)
    endif()
  endif()
endif()

if(NOT SUPPORTS_CHAR16_T)
  message("SUPPORTS: ${SUPPORTS_CHAR16_T}")
  message("OUTPUT\n${TEST_CHAR16_T_TEST_OUTPUT}")

  message(FATAL_ERROR "ERROR: Your compiler settings do not support char16_t\nFor clang try to add CXXFLAGS=\"-stdlib=libc++ -std=c++11\" when running cmake")
endif()

