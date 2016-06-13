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
typedef itk::PointSet<vtkFloatingPointType,3>          PointSetType;
//typedef PointSetType::PointType                        PointType;
typedef std::list<FiberPointType>                      FiberType;
typedef std::list<FiberType>                           FiberListType;
typedef std::list<vtkPolyData *>                       VTKFiberListType;
typedef double                                     TensorElementType;
typedef itk::DiffusionTensor3D<TensorElementType>   TensorPixelType;
typedef itk::Image<TensorPixelType,3>               TensorImageType;
typedef unsigned char                             MaskPixelType;
typedef itk::Image<MaskPixelType,3>                 MaskImageType;
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

template <class TTensorImageType, class TAnisotropyImageType, class TMaskImageType>
class DtiTrackingFilterBase : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef DtiTrackingFilterBase    Self;
  typedef itk::Object              Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TTensorImageType                        TensorImageType;
  typedef typename TensorImageType::Pointer       TensorImagePointer;
  typedef typename TensorImageType::ConstPointer  TensorImageConstPointer;
  typedef typename TensorImageType::RegionType    TensorImageRegionType;
  typedef typename TensorImageType::SizeType      TensorImageSizeType;
  typedef typename TensorImageType::SpacingType   TensorImageSpacingType;
  typedef typename TensorImageType::PointType     TensorImagePointType;
  typedef typename TensorImageType::PixelType     TensorImagePixelType;
  typedef typename TensorImageType::DirectionType TensorImageDirectionType;

  typedef TAnisotropyImageType                        AnisotropyImageType;
  typedef typename AnisotropyImageType::Pointer       AnisotropyImagePointer;
  typedef typename AnisotropyImageType::ConstPointer  AnisotropyImageConstPointer;
  typedef typename AnisotropyImageType::RegionType    AnisotropyImageRegionType;
  typedef typename AnisotropyImageType::SizeType      AnisotropyImageSizeType;
  typedef typename AnisotropyImageType::SpacingType   AnisotropyImageSpacingType;
  typedef typename AnisotropyImageType::PointType     AnisotropyImagePointType;
  typedef typename AnisotropyImageType::PixelType     AnisotropyImagePixelType;
  typedef typename AnisotropyImageType::DirectionType AnisotropyImageDirectionType;

  typedef TMaskImageType                        MaskImageType;
  typedef typename MaskImageType::Pointer       MaskImagePointer;
  typedef typename MaskImageType::ConstPointer  MaskImageConstPointer;
  typedef typename MaskImageType::RegionType    MaskImageRegionType;
  typedef typename MaskImageType::SizeType      MaskImageSizeType;
  typedef typename MaskImageType::SpacingType   MaskImageSpacingType;
  typedef typename MaskImageType::PointType     MaskImagePointType;
  typedef typename MaskImageType::PixelType     MaskImagePixelType;
  typedef typename MaskImageType::DirectionType MaskImageDirectionType;

  typedef itk::LinearInterpolateImageFunction<MaskImageType, double>         MaskIPType;
  typedef itk::LinearInterpolateImageFunction<AnisotropyImageType, double>   ScalarIPType;
  typedef itk::TensorLinearInterpolateImageFunction<TensorImageType, double> VectorIPType;

  typedef typename itk::ContinuousIndex<double, 3> ContinuousIndexType;

  typedef itk::Point<double, 3>          PointType;
  typedef std::list<ContinuousIndexType> SeedListType;
  typedef std::list<BranchPointType>     BranchListType;
  typedef std::list<TVector>             DirectionListType;

  typedef itk::PointSet<double, 3>       PointSetType;
  typedef vtkPolyData *                  DtiFiberType;

  /** ImageDimension constants * /
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  / ** The dimensions of the input image must equal those of the
      output image. * /
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),itkGetStaticConstMacro(OutputImageDimension)>));

  / ** The dimension of the input image must be 4. * /
  itkConceptMacro(DimensionShouldBe4,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),4>));
*/
  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(DtiTrackingFilterBase, itk::Object);

  itkSetObjectMacro(TensorImage,  TensorImageType);
  itkSetObjectMacro(AnisotropyImage,  AnisotropyImageType);
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

  DtiFiberType GetOutput();

  // void Update();
protected:
  DtiTrackingFilterBase();
  ~DtiTrackingFilterBase()
  {
  }

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiTrackingFilterBase);

protected:
  bool IsLoop(vtkPoints *fiber, double tolerance = 0.001);

  void InitializeSeeds();

  void ContinuousIndexToMM(ContinuousIndexType & index, PointType & p);

  void MMToContinuousIndex(PointType & p, ContinuousIndexType & index);

  void MMToContinuousIndex(double *p, ContinuousIndexType & index);

  void StepIndexInPointSpace(ContinuousIndexType & newIndex, ContinuousIndexType & oldIndex, TVector & vec);

  void StepIndex(ContinuousIndexType & newIndex, ContinuousIndexType & oldIndex, TVector & vec);

  void ApplyTensorDeflection(TVector & vin, TMatrix & fullTensorPixel, TVector & e2, TVector & vout);

  void AddFiberToOutput( vtkPoints *currentFiber, vtkFloatArray *fiberTensors );

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
  typename MaskIPType::Pointer m_StartIP;
  typename MaskIPType::Pointer m_EndIP;

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
};  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiTrackingFilterBase.hxx"
#endif

#endif
