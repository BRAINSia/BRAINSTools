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
#ifndef __itkSimpleDiffeomorphicRegistration_h
#define __itkSimpleDiffeomorphicRegistration_h

#include <string>
#include "itkImage.h"
#include "itkObject.h"
// #include "BRAINSDemonWarp.h"
#include "itkVector.h"
#include "DemonsPreprocessor.h"
#include "DemonsRegistrator.h"

#define DIM 3
#define MaxStepLength 2
#define SmoothDisplacementFieldSigma 1.5
#define NumberOfLevels 5
#define NumberOfIteration0 300
#define NumberOfIteration1 100 // 100
#define NumberOfIteration2 30  // 30
#define NumberOfIteration3 20  // 20
#define NumberOfIteration4 15  // 15
#define FixedPyramid 16
#define NumberOfMatchPoints 7
#define NumberOfHistogramLevels 1024

/** INFO:  Need to document this class
 */
class itkSimpleDiffeomorphicRegistration : public itk::Object
{
public:
  using TRealImage = itk::Image<float, DIM>;
  using DemonsPreprocessorType = itk::DemonsPreprocessor<TRealImage, TRealImage>;
  using DemonsRegistratorType = itk::DemonsRegistrator<TRealImage, TRealImage, float>;
  using TDisplacementField = itk::Image<itk::Vector<float, DIM>, DIM>;

  itkSimpleDiffeomorphicRegistration();
  itkSetObjectMacro(FixedImage, TRealImage);
  itkSetObjectMacro(MovingImage, TRealImage);

  itkSetStringMacro(DisplacementFieldName);
  itkSetStringMacro(DeformedImageName);
  itkGetStringMacro(DeformedImageName);
  itkGetConstObjectMacro(DisplacementField, TDisplacementField);
  void
  Update();

protected:
  // std::string GetFixedImage(void);
  // std::string GetMovingImage(void);
  // std::string GetDeformedImageName(void);
  // std::string GetDisplacementPrefixName(void);
  void
  InitializePreprocessor();

  void
  Initialization(void);

private:
  DemonsPreprocessorType::Pointer m_DemonsPreprocessor;
  DemonsRegistratorType::Pointer  m_DemonsRegistrator;
  TDisplacementField::Pointer     m_DisplacementField;
  TRealImage::Pointer             m_FixedImage;
  TRealImage::Pointer             m_MovingImage;
  std::string                     m_DeformedImageName;
  std::string                     m_DisplacementFieldName;
};

#endif
