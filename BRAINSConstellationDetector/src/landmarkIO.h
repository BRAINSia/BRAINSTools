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
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __landmarkIO__h
#define __landmarkIO__h

#include "landmarksConstellationCommon.h"
#include "landmarksConstellationDetector.h"
#include "landmarksConstellationModelIO.h"
#include "itkVersorRigid3DTransform.h"

#include <cstring>
#include <cstdio>
#include <map>
#include <vector>

namespace
{
typedef vnl_matrix<double>                           MatrixType;
typedef vnl_vector<double>                           VectorType;
typedef itk::VersorRigid3DTransform<double>          VersorTransformType;
typedef VersorTransformType::Pointer                 VersorTransformPointer;
typedef VersorTransformType::MatrixType              VersorTransformMatrixType;
typedef itk::ImageDuplicator<SImageType>             DuplicatorType;
}

extern void MakeBrandeddebugImage(SImageType::ConstPointer in, const landmarksConstellationModelIO & mDef,
                                  const SImageType::PointType & RP, const SImageType::PointType & AC,
                                  const SImageType::PointType & PC, const SImageType::PointType & VN4,
                                  const std::string & fname, const SImageType::PointType & RP2,
                                  const SImageType::PointType & AC2, const SImageType::PointType & PC2,
                                  const SImageType::PointType & VN42);

extern void MakePointBranded3DImage(SImageType::ConstPointer in, const SImageType::PointType & CenterPoint,
                                    const std::string & fname);

extern void MakeBranded2DImage(SImageType::ConstPointer in, landmarksConstellationDetector & myDetector,
                               const SImageType::PointType & RP, const SImageType::PointType & AC,
                               const SImageType::PointType & PC, const SImageType::PointType & VN4,
                               const SImageType::PointType & CM, const std::string & fname);

// Write Slicer scene file (.mrml)
extern void WriteMRMLFile(std::string outputMRML, std::string outputLandmarksInInputSpace,
                          std::string outputLandmarksInOutputSpace, std::string inputVolume, std::string outputVolume,
                          std::string outputTransform, const LandmarksMapType & outputLandmarksInInputSpaceMap,
                          const LandmarksMapType & outputLandmarksInOutputSpaceMap,
                          VersorTransformType::ConstPointer versorTransform);

// load linear least squares model for selected landmarks
// .load from txt file
extern void loadLLSModel(std::string llsModel, std::map<std::string, std::vector<double> > & llsMeans,
                         std::map<std::string, MatrixType> & llsMatrices, std::map<std::string, double> & searchRadii);

/*
// load from .mat file
extern void loadLLSModelMat(std::string llsModel, std::string processingList,
                            std::map<std::string,
                                     std::vector<double> > & llsMeans,
                            std::map<std::string, MatrixType> & llsMatrices, std::map<std::string,
                                                                                      double>
                            & searchRadii);
*/
extern void writeLLSModel(const std::string & modelName, const std::map<std::string, std::vector<double> > & llsMeans,
                          const std::map<std::string, MatrixType> & llsMatrices, const std::map<std::string,
                                                                                                double> & searchRadii);

extern void readLLSModel(const std::string & modelName, std::map<std::string, std::vector<double> > & llsMeans,
                         std::map<std::string, MatrixType> & llsMatrices, std::map<std::string, double> & searchRadii);

// write out verification script
extern void writeVerificationScript(std::string outputVerificationScript, std::string outputVolume,
                                    std::string saveOutputLandmarks);

#endif // __landmarkIO__h
