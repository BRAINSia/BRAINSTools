/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
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
#include <BRAINSTypes.h>

/*
 * This IO utility program write and read the ITK landmarks weight.
 */


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
