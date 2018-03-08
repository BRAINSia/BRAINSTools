##- project( liblinear )

add_library( blas blas/dnrm2.c blas/daxpy.c blas/ddot.c blas/dscal.c blas/blas.h blas/blas.h blas/blasp.h )
set_target_properties(blas PROPERTIES FOLDER ${MODULE_FOLDER})

add_library( tron tron.cpp tron.h )
target_link_libraries(tron blas)
set_target_properties(tron PROPERTIES FOLDER ${MODULE_FOLDER})

add_library( liblinear linear.cpp linear.h )
target_link_libraries(liblinear tron ${BRAINSContinuousClass_ITK_LIBRARIES})
set_target_properties(liblinear PROPERTIES FOLDER ${MODULE_FOLDER})

