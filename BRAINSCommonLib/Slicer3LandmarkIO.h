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
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __Slicer3LandmarkIO__h
#define __Slicer3LandmarkIO__h

#include "itkPoint.h"
#include <itksys/SystemTools.hxx>

#include <fstream>
#include <map>
#include <string>

#include <BRAINSTypes.h>

/*
 * This IO utility program transforms between ITK landmarks (in LPS coordinate
 * system) to Slicer3 landmarks (in RAS coordinate system).
 */


/*
 * Read lmk weights
 * weightFilename  -
 * Output: A map of weights
 */
LandmarkWeightMapType ReadLandmarkWeights( const std::string & weightFilename );

/*
 * Write ITK landmarks to a Slicer3 landmark list file (.fcsv)
 * Input:
 * landmarksFilename  - the filename of the output Slicer landmark list file
 * landmarks          - a map of landmarks (itkPoints) to be written into file
 *
 * Output:
 * NONE
 */
extern void WriteITKtoSlicer3Lmk( const std::string & landmarksFilename, const LandmarksMapType & landmarks );

/*
 * Read Slicer3 landmark list file (.fcsv) into a map of ITK points
 * Input:
 * landmarksFilename  - the filename of the input Slicer landmark list file
 *
 * Output:
 * landmarks          - a map of itkPoints to save the landmarks in ITK
 */
extern LandmarksMapType ReadSlicer3toITKLmk( const std::string & landmarksFilename );

#endif
