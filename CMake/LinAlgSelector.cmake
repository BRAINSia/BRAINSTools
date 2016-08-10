##==============================================================
## Switch based on the linear algebra optimized library to
## use.  Note that the first library successfully found
## will be used
##
## -*-*- Try OpenBLAS first
if(NOT LINALG_VENDOR OR LINALG_VENDOR MATCHES "OpenBLAS")
    include( ${CMAKE_CURRENT_LIST_DIR}/FindOpenBLAS.cmake )
    if(OpenBLAS_FOUND)
        set(LINALG_VENDOR_FOUND TRUE)
        set(LINALG_VENDOR "OpenBLAS")
        set(LINALG_INCLUDE_DIRS ${OpenBLAS_INCLUDE_DIRS})
        if(OpenBLAS_HAS_PARALLEL_LIBRARIES)
            # HACK!! set(LINALG_LIBRARIES ${OpenBLAS_PARALLEL_LIBRARIES})
            set(LINALG_LIBRARIES ${OpenBLAS_LIBRARIES})
        else()
            set(LINALG_LIBRARIES ${OpenBLAS_LIBRARIES})
        endif()
    endif()
endif()

##
## -*-*- Try ATLAS version next
if(NOT LINALG_VENDOR OR LINALG_VENDOR MATCHES "ATLAS")
    include( ${CMAKE_CURRENT_LIST_DIR}/FindATLAS.cmake )
    if(ATLAS_FOUND)
        set(LINALG_VENDOR_FOUND TRUE)
        set(LINALG_VENDOR "ATLAS")
        set(LINALG_INCLUDE_DIRS ${ATLAS_INCLUDE_DIRS})
        set(LINALG_LIBRARIES ${ATLAS_LIBRARIES})
    endif()
endif()

##
## -*-*- Try Generic LAPACKE version Last
if(NOT LINALG_VENDOR OR LINALG_VENDOR MATCHES "LAPACKE")
    #NOTE: By specifying Fortran here, linking to lapack becomes easier
    # See https://blog.kitware.com/fortran-for-cc-developers-made-easier-with-cmake/
    enable_language(Fortran)
    ## Only very new versions of LAPACK (> 3.5.0) have built in support
    ## for cblas and lapacke.  This method is not very robust to older
    ## versions of lapack that might be able to be supported.
    ## It is know to work local builds
    include( ${CMAKE_CURRENT_LIST_DIR}/FindLAPACKE.cmake )
    if(LAPACKE_FOUND)
        set(LINALG_VENDOR_FOUND TRUE)
        set(LINALG_VENDOR "LAPACKE")
        set(LINALG_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIRS})
        set(LINALG_LIBRARIES ${LAPACKE_LIBRARIES})
    endif()
endif()

##
## -*-*- Finally, set include_directories -*-*-*
if(NOT LINALG_VENDOR_FOUND)
    message(FATAL_ERROR "No valid linear algebra libraries found!")
endif()
include_directories(${LINALG_INCLUDE_DIRS})