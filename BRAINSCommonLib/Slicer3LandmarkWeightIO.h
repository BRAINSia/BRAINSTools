/*
 * Author: Ali Ghayoor, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

#ifndef __Slicer3LandmarkWeightIO__h
#define __Slicer3LandmarkWeightIO__h

#include <itksys/SystemTools.hxx>

#include "itkPoint.h"

#include <fstream>
#include <map>
#include <string>

/*
 * This IO utility program write and read the ITK landmarks weight.
 */

typedef std::map<std::string, double> LandmarksWeightMapType;

/*
 * Write ITK landmark weights to a Slicer3 landmark list file (.wts)
 * Input:
 * landmarksWeightFilename  - the filename of the output Slicer landmark list file
 * landmarks          - a map of landmark weights (double values) to be written into file
 *
 * Output:
 * NONE
 */
extern void WriteITKtoSlicer3LmkWts( const std::string & landmarksWeightFilename,
                                     const LandmarksWeightMapType & landmarks );

/*
 * Read Slicer3 landmark list file (.wts) into a map of doubles
 * Input:
 * landmarksWeightFilename  - the filename of the input Slicer landmark list file
 *
 * Output:
 * landmark weights          - a map of doubles to save the landmark wights in ITK
 */
extern LandmarksWeightMapType ReadSlicer3toITKLmkWts( const std::string & landmarksWeightFilename );

#endif
