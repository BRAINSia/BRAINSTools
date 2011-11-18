project( liblinear )

cmake_minimum_required(VERSION 2.8)

add_library( tron tron.cpp tron.h )
add_library( liblinear linear.cpp linear.h blas/dnrm2.c blas/daxpy.c blas/ddot.c blas/dscal.c blas/blas.h blas/blas.h blas/blasp.h )
target_link_libraries(liblinear tron)

