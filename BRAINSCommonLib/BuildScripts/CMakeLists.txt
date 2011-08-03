### This is just the absolute minimum script necessary
### to define a separate project that contains all the
### goodie scripts needed to support the biulding
### the BRAINS* tools

project(BRAINSBuildScripts)
cmake_minimum_required(VERSION 2.8)

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

install( FILES License.txt DESTINATION ${CMAKE_INSTALL_PREFIX} )

