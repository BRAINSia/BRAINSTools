NAMICExternalProjects
=====================

A superbuild infrastructure to build all superbuild structures


This project acts as a template to model building all other SuperBuild projects around.

NAMICExternalProjects/
        CMakeLists.txt       | A file that acts as a top level switch, it runs the 
                               building of pre-requisite ExternalProjects Defined in
                               SuperBuild.cmake (where ${MYProject} is last ExternalProject
                               to be built).  This file is processed by cmake 2 times, the
                               first time it processes SuperBuild.cmake, the second time
                               it processes ${MyProject}.cmake
        Common.cmake         | Included in SuperBuild.cmake and ${MYProject}.cmake for common
                               conditional compilation options to be set.  These options
                               must be passed to all ExternalProjects (including ${MYProject})
        SuperBuild.cmake     | Infrastructure to build all the pre-requisite packages where
                               ${MYProject} is the last one listed
        ${MYProject}.cmake   | Standard cmake build instructions for ${MYProject}
        SuperBuild/          | A directory full of External_${extProjName}.cmake files defining
                               how to build external dependancies.
        CMake/               | A directory of support files.


