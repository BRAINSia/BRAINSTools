
include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/SparseLibMVIml)
link_directories(${SPHARMP_BINARY_DIR}/bin)

# here prepared only for ISO C++ compatibles compilers

# IBM xlC  v. 1.1
# CCCFLAGS        = -DCOMPLEX=complex
# LDFLAGS         = -lm -lcomplex

# Sun C++ 4.0.1
# CCCFLAGS        = -DMV_VECTOR_BOUNDS_CHECK -g -DCOMPLEX_OSTREAM -DCOMPLEX=complex
# LDFLAGS         =   -lm -lcomplex

# g++ v. 2.6.3
# CCCFLAGS                =   -O5 '-DCOMPLEX=std::complex<double>'
# LDFLAGS                 =   -lm


if(WIN32)
add_definitions("-DCOMPLEX=std::complex<double>")
else()
add_definitions('-DCOMPLEX=std::complex<double>')
endif()

link_libraries(SparseMatrixLib)

