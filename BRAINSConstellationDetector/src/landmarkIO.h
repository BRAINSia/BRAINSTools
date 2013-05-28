/*
 * Author: Han J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef __landmarkIO__h
#define __landmarkIO__h

#include "landmarksConstellationCommon.h"
#include "landmarksConstellationDetector.h"
#include "itkOtsuThresholdImageFilter.h"
#include "landmarksConstellationModelIO.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkImageFileReader.h"
#include "itkImageDuplicator.h"

#include <cstring>
#include <map>
#include <cstdio>
#include <vector>
#include <vnl/vnl_matlab_read.h>

typedef std::map<std::string, SImageType::PointType> LandmarksMapType;
typedef vnl_matrix<double>                           MatrixType;
typedef vnl_vector<double>                           VectorType;
typedef itk::VersorRigid3DTransform<double>          VersorTransformType;
typedef VersorTransformType::Pointer                 VersorTransformPointer;
typedef VersorTransformType::MatrixType              VersorTransformMatrixType;
typedef itk::ImageDuplicator<SImageType>             DuplicatorType;

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
