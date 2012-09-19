##- project( liblinear )

cmake_minimum_required(VERSION 2.8)

add_library( blas blas/dnrm2.c blas/daxpy.c blas/ddot.c blas/dscal.c blas/blas.h blas/blas.h blas/blasp.h )
add_library( tron tron.cpp tron.h )
target_link_libraries(tron blas)
add_library( liblinear linear.cpp linear.h )
target_link_libraries(liblinear tron)

