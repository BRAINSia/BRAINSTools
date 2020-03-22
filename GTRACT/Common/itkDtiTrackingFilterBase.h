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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDtiTrackingFilterBase_h
#define __itkDtiTrackingFilterBase_h
#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
// #include <vtkXMLImageDataWriter.h>
#include <vtkAppendPolyData.h>

#include "itkObject.h"
#include "itkImage.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
// #include "itkIOCommon.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkPoint.h"
// #include "itkBlobSpatialObject.h"

// #include <metaCommand.h>
#include <itkImage.h>
#include <itkVector.h>
// #include <itkBlobSpatialObject.h>
// #include <itkSceneSpatialObject.h>
// #include <itkSpatialObjectReader.h>
// #include <itkEuclideanDistancePointMetric.h>
// #include <itkIdentityTransform.h>
#include <itkDiffusionTensor3D.h>

#include "algo.h"
#include "GtractTypes.h"
#include "itkTensorLinearInterpolateImageFunction.h"

#include <map>
#include <string>

// ////////////////////////////////////////////////////////////////////////

/*
using PointSetType = itk::PointSet<vtkFloatingPointType,3>;
//using PointType = PointSetType::PointType;
using FiberType = std::list<FiberPointType>;
using FiberListType = std::list<FiberType>;
using VTKFiberListType = std::list<vtkPolyData *>;
using TensorElementType = double;
using TensorPixelType = itk::DiffusionTensor3D<TensorElementType>;
using TensorImageType = itk::Image<TensorPixelType,3>;
using MaskPixelType = unsigned char;
using MaskImageType = itk::Image<MaskPixelType,3>;
*/

namespace itk
{
/** \class DtiTrackingFilterBase
 *  \brief base class for
 *         DtiFreeTrackingFilter,
 *         DtiStreamlineTrackingFilter,
 *         DtiGraphSearchTrackingFilter,
 *         DtiGuidedTrackingFilter,
 */

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
class DtiTrackingFilterBase : public itk::Object
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiTrackingFilterBase);

  /** Standard class type alias. */
  using Self = DtiTrackingFilterBase;
  using Superclass = itk::Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Some convenient type alias. */
  using TensorImageType = TTensorImageType;
  using TensorImagePointer = typename TensorImageType::Pointer;
  using TensorImageConstPointer = typename TensorImageType::ConstPointer;
  using TensorImageRegionType = typename TensorImageType::RegionType;
  using TensorImageSizeType = typename TensorImageType::SizeType;
  using TensorImageSpacingType = typename TensorImageType::SpacingType;
  using TensorImagePointType = typename TensorImageType::PointType;
  using TensorImagePixelType = typename TensorImageType::PixelType;
  using TensorImageDirectionType = typename TensorImageType::DirectionType;

  using AnisotropyImageType = TAnisotropyImageType;
  using AnisotropyImagePointer = typename AnisotropyImageType::Pointer;
  using AnisotropyImageConstPointer = typename AnisotropyImageType::ConstPointer;
  using AnisotropyImageRegionType = typename AnisotropyImageType::RegionType;
  using AnisotropyImageSizeType = typename AnisotropyImageType::SizeType;
  using AnisotropyImageSpacingType = typename AnisotropyImageType::SpacingType;
  using AnisotropyImagePointType = typename AnisotropyImageType::PointType;
  using AnisotropyImagePixelType = typename AnisotropyImageType::PixelType;
  using AnisotropyImageDirectionType = typename AnisotropyImageType::DirectionType;

  using MaskImageType = TMaskImageType;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using MaskImageConstPointer = typename MaskImageType::ConstPointer;
  using MaskImageRegionType = typename MaskImageType::RegionType;
  using MaskImageSizeType = typename MaskImageType::SizeType;
  using MaskImageSpacingType = typename MaskImageType::SpacingType;
  using MaskImagePointType = typename MaskImageType::PointType;
  using MaskImagePixelType = typename MaskImageType::PixelType;
  using MaskImageDirectionType = typename MaskImageType::DirectionType;

  using MaskIPType = itk::LinearInterpolateImageFunction<MaskImageType, double>;
  using ScalarIPType = itk::LinearInterpolateImageFunction<AnisotropyImageType, double>;
  using VectorIPType = itk::TensorLinearInterpolateImageFunction<TensorImageType, double>;

  using ContinuousIndexType = typename itk::ContinuousIndex<double, 3>;

  using PointType = itk::Point<double, 3>;
  using SeedListType = std::list<ContinuousIndexType>;
  using BranchListType = std::list<BranchPointType>;
  using DirectionListType = std::list<TVector>;

  using PointSetType = itk::PointSet<double, 3>;
  using DtiFiberType = vtkPolyData *;

  /** ImageDimension constants * /
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  / ** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<Self::InputImageDimension,Self::OutputImageDimension>));

  / ** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<Self::InputImageDimension,4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DtiTrackingFilterBase, itk::Object);

  itkSetObjectMacro(TensorImage, TensorImageType);
  itkSetObjectMacro(AnisotropyImage, AnisotropyImageType);
  itkSetObjectMacro(StartingRegion, MaskImageType);
  itkSetObjectMacro(EndingRegion, MaskImageType);

  itkSetMacro(SeedThreshold, float);
  itkSetMacro(AnisotropyThreshold, float);
  itkSetMacro(MaximumLength, float);
  itkSetMacro(MinimumLength, float);
  itkSetMacro(StepSize, float);
  itkSetMacro(UseLoopDetection, bool);
  itkSetMacro(UseTend, bool);
  itkSetMacro(TendG, float);
  itkSetMacro(TendF, float);

  DtiFiberType
  GetOutput();

  // void Update();
protected:
  DtiTrackingFilterBase();
  ~DtiTrackingFilterBase() override = default;

private:
protected:
  bool
  IsLoop(vtkPoints * fiber, double tolerance = 0.001);

  void
  InitializeSeeds();

  void
  ContinuousIndexToMM(ContinuousIndexType & index, PointType & p);

  void
  MMToContinuousIndex(PointType & p, ContinuousIndexType & index);

  void
  MMToContinuousIndex(double * pt, ContinuousIndexType & index);

  void
  StepIndexInPointSpace(ContinuousIndexType & newIndex, ContinuousIndexType & oldIndex, TVector & vec);

  void
  StepIndex(ContinuousIndexType & newIndex, ContinuousIndexType & oldIndex, TVector & vec);

  void
  ApplyTensorDeflection(TVector & vin, TMatrix & fullTensorPixel, TVector & e2, TVector & vout);

  void
  AddFiberToOutput(vtkPoints * currentFiber, vtkFloatArray * fiberTensors);

  DirectionListType m_TrackingDirections;

  // Input and Output Image
  TensorImagePointer     m_TensorImage;
  AnisotropyImagePointer m_AnisotropyImage;
  DtiFiberType           m_Output;
  SeedListType           m_Seeds;
  MaskImagePointer       m_StartingRegion;
  MaskImagePointer       m_EndingRegion;

  // Interpolation data:  The Vector is the Tensor image, the Scalar is the
  // Anisotropy image.
  typename ScalarIPType::Pointer m_ScalarIP;
  typename VectorIPType::Pointer m_VectorIP;
  typename MaskIPType::Pointer   m_StartIP;
  typename MaskIPType::Pointer   m_EndIP;

  float m_SeedThreshold;
  float m_AnisotropyThreshold;
  float m_MaximumLength;
  float m_MinimumLength;
  float m_StepSize;
  bool  m_UseLoopDetection;
  bool  m_UseTend;
  float m_TendG;
  float m_TendF;

  float pi;
}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDtiTrackingFilterBase.hxx"
#endif

#endif
