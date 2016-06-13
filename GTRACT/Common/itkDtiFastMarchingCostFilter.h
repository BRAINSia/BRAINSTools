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

#ifndef __itkDtiFastMarchingCostFilter_h
#define __itkDtiFastMarchingCostFilter_h

#include <iostream>
#include <fstream>
#include <itkImage.h>
#include <itkObject.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include "itkProcessObject.h"
#include <itkDiffusionTensor3D.h>
#include <itkConstNeighborhoodIterator.h>

#include "GtractTypes.h"
#include <map>
#include <string>

#include <itkLevelSet.h>
#include <itkIndex.h>
#include <vnl/vnl_math.h>

#include <functional>
#include <queue>

namespace itk
{
template <
  class TLevelSet,
  class TTensorImage = Image<itk::DiffusionTensor3D<float>, 3> >
// class TTensorImage=
//
// Image<itk::DiffusionTensor3D<float>,TLevelSet::ImageDimension>
// >

class DtiFastMarchingCostFilter :
  public ImageToImageFilter<TTensorImage, TLevelSet>
{
public:

  /** Standard class typdedefs. */
  typedef DtiFastMarchingCostFilter Self;
  typedef ImageSource<TLevelSet>    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DtiFastMarchingCostFilter, ImageToImageFilter);

  /** Typedef support of level set method types for output image. */
  typedef LevelSetTypeDefault<TLevelSet>              LevelSetType;
  typedef typename LevelSetType::LevelSetImageType    LevelSetImageType;
  typedef typename LevelSetType::LevelSetPointer      LevelSetPointer;
  typedef typename LevelSetType::PixelType            PixelType;
  typedef typename LevelSetType::NodeType             NodeType;
  typedef typename LevelSetType::NodeContainer        NodeContainer;
  typedef typename LevelSetType::NodeContainerPointer NodeContainerPointer;
  typedef typename LevelSetImageType::SizeType        OutputSizeType;
  typedef typename LevelSetImageType::RegionType      OutputRegionType;
  typedef typename LevelSetImageType::SpacingType     OutputSpacingType;
  typedef typename LevelSetImageType::PointType       OutputPointType;
  typedef typename LevelSetImageType::DirectionType   OutputDirectionType;

  class AxisNodeType : public NodeType
  {
public:
    int GetAxis() const
    {
      return m_Axis;
    }

    void SetAxis( int axis )
    {
      m_Axis = axis;
    }

    const AxisNodeType & operator=(const NodeType & node)
    {
      this->NodeType::operator=(node); return *this;
    }

private:
    int m_Axis;
  };

  enum { dimension = 3 };

  typedef vnl_vector_fixed<float, dimension> TVector;

  /** Index typedef support. */
  typedef Index<dimension> IndexType;

  /** Typedef support of input Tensor Image Type */
  typedef TTensorImage                            TensorImageType;
  typedef typename TensorImageType::Pointer       TensorImagePointer;
  typedef typename TensorImageType::ConstPointer  TensorImageConstPointer;
  typedef typename TensorImageType::RegionType    TensorImageRegionType;
  typedef typename TensorImageType::SizeType      TensorImageSizeType;
  typedef typename TensorImageType::SpacingType   TensorImageSpacingType;
  typedef typename TensorImageType::PointType     TensorImagePointType;
  typedef typename TensorImageType::PixelType     TensorImagePixelType;
  typedef typename TensorImageType::IndexType     TensorImageIndexType;
  typedef typename TensorImageType::DirectionType TensorImageDirectionType;

  /** Typedef support of input Anisotropy Image Type */

  typedef itk::Image<float, dimension>                 AnisotropyImageType;
  typedef typename  AnisotropyImageType::Pointer       AnisotropyImagePointer;
  typedef typename  AnisotropyImageType::ConstPointer  AnisotropyImageConstPointer;
  typedef typename  AnisotropyImageType::RegionType    AnisotropyImageRegionType;
  typedef typename  AnisotropyImageType::SizeType      AnisotropyImageSizeType;
  typedef typename  AnisotropyImageType::SpacingType   AnisotropyImageSpacingType;
  typedef typename  AnisotropyImageType::PointType     AnisotropyImagePointType;
  typedef typename  AnisotropyImageType::PixelType     AnisotropyImagePixelType;
  typedef typename  AnisotropyImageType::DirectionType AnisotropyImageDirectionType;

  /** Enum of Fast Marching algorithm point types. FarPoints represent far
   * away points; TrialPoints represent points within a narrowband of the
   * propagating front; and AlivePoints represent points which have already
   * been processed. */
  enum LabelType { FarPoint, AlivePoint, TrialPoint };

  /** LabelImage typedef support. */
  typedef Image<unsigned char, dimension> LabelImageType;

  /** LabelImagePointer typedef support. */
  typedef typename LabelImageType::Pointer LabelImagePointer;

  /** OutputSpeedImage typedef support. */
  typedef itk::Image<float, dimension>             OutputSpeedImageType;
  typedef typename OutputSpeedImageType::Pointer   OutputSpeedImagePointer;
  typedef typename OutputSpeedImageType::PixelType SpeedImagePixelType;

  typedef itk::Vector<float, 3>               EigenvectorPixelType;
  typedef itk::Image<EigenvectorPixelType, 3> EigenvectorImageType;
  // typedef itk::Image<EigenvectorPixelType, dimension> EigenvectorImageType;
  typedef typename EigenvectorImageType::Pointer               EigenvectorImagePointer;
  typedef itk::ConstNeighborhoodIterator<EigenvectorImageType> ConstNeighborhoodIteratorType;

  /** Set the container of Alive Points representing the initial front.
   * Alive points are represented as a VectorContainer of LevelSetNodes. */
  void SetAlivePoints( NodeContainer *points )
  {
    m_AlivePoints = points;
    this->Modified();
  }

  /** Get the container of Alive Points representing the initial front. */
  NodeContainerPointer GetAlivePoints()
  {
    return m_AlivePoints;
  }

  /** Set the container of Trial Points representing the initial front.
   * Trial points are represented as a VectorContainer of LevelSetNodes. */
  void SetTrialPoints( NodeContainer *points )
  {
    m_TrialPoints = points;
    this->Modified();
  }

  /** Get the container of Trial Points representing the initial front. */
  NodeContainerPointer GetTrialPoints()
  {
    return m_TrialPoints;
  }

  /** Get the point type label image. */
  LabelImagePointer GetLabelImage() const
  {
    return m_LabelImage;
  }

  /** Get the eigenvector image. */
  EigenvectorImagePointer GetEigenvectorImage() const
  {
    return m_EigenvectorImage;
  }

  /** Get the output speed image. */
  OutputSpeedImagePointer GetOutputSpeedImage() const
  {
    return m_OutputSpeedImage;
  }

  itkSetObjectMacro(OutputImage,  LevelSetImageType);
  itkGetConstObjectMacro(OutputImage,  LevelSetImageType);
  itkSetObjectMacro(TensorImage,  TensorImageType);
  itkGetConstObjectMacro(TensorImage,  TensorImageType);
  itkSetObjectMacro(AnisotropyImage,  AnisotropyImageType);
  itkGetConstObjectMacro(AnisotropyImage,  AnisotropyImageType);

  /** Set/Get the Normalization Factor for the Speed Image.
      The values in the Speed Image is divided by this
      factor. This allows the use of images with
      integer pixel types to represent the speed. */

  itkSetMacro( NormalizationFactor, double );
  itkGetMacro( NormalizationFactor, double );

  /** Set the Fast Marching algorithm Stopping Value. The Fast Marching
   * algorithm is terminated when the value of the smallest trial point
   * is greater than the stopping value. */
  itkSetMacro( StoppingValue, double );

  /** Set Weight for Anisotropy Image*/
  itkSetMacro( AnisotropyWeight, float);

  /** Get the Fast Marching algorithm Stopping Value. */
  itkGetConstReferenceMacro( StoppingValue, double );

  /** Set the Collect Points flag. Instrument the algorithm to collect
   * a container of all nodes which it has visited. Useful for
   * creating Narrowbands for level set algorithms that supports
   * narrow banding. */
  itkSetMacro( CollectPoints, bool );

  /** Get thConste Collect Points flag. */
  itkGetConstReferenceMacro( CollectPoints, bool );
  itkBooleanMacro( CollectPoints );

  /** Get the container of Processed Points. If the CollectPoints flag
   * is set, the algorithm collects a container of all processed nodes.
   * This is useful for defining creating Narrowbands for level
   * set algorithms that supports narrow banding. */
  NodeContainerPointer GetProcessedPoints() const
  {
    return m_ProcessedPoints;
  }

  // DON'T NEED TO MANUALLY SET OUTPUT IMAGE!!!!
  /** The output largeset possible, spacing and origin is computed as follows.
   * If the speed image is NULL or if the OverrideOutputInformation is true,
   * the output information is set from user specified parameters. These
   * parameters can be specified using methods SetOutputRegion(), SetOutputSpacing()
   * and SetOutputOrigin(). Else if the speed image is not NULL, the output information
   * is copied from the input speed image. */

  virtual void SetOutputSize( const OutputSizeType & size )
  {
    m_OutputRegion = size;
  }

  virtual OutputSizeType GetOutputSize() const
  {
    return m_OutputRegion.GetSize();
  }

  itkSetMacro( OutputRegion, OutputRegionType );
  itkGetConstReferenceMacro( OutputRegion, OutputRegionType );
  itkSetMacro( OutputSpacing, OutputSpacingType );
  itkGetConstReferenceMacro( OutputSpacing, OutputSpacingType );
  itkSetMacro( OutputOrigin, OutputPointType );
  itkGetConstReferenceMacro( OutputOrigin, OutputPointType );
  itkSetMacro( OutputDirection, OutputDirectionType);
  itkGetConstReferenceMacro( OutputDirection, OutputDirectionType );

  itkSetMacro( OverrideOutputInformation, bool );
  itkGetConstReferenceMacro( OverrideOutputInformation, bool );
  itkBooleanMacro( OverrideOutputInformation );

  void GenerateData() ITK_OVERRIDE;

protected:
  DtiFastMarchingCostFilter();
  ~DtiFastMarchingCostFilter()
  {
  }

  void PrintSelf( std::ostream & os, Indent indent ) const ITK_OVERRIDE;

  virtual void UpdateFront( // const TensorImageType *,
    LevelSetImageType * );

  virtual double InitializeTrialPoints( // const TensorImageType *,
    const IndexType & index, LevelSetImageType * );

  virtual void UpdateNeighbors( // const TensorImageType *,
    const IndexType & index, LevelSetImageType * );

  virtual double UpdateValue( // const TensorImageType *,
    const IndexType & index, LevelSetImageType * );

  virtual double ModifiedUpdateValue( const TensorImageType *, const IndexType & index, LevelSetImageType * );

  /** Set Image NRRD Meta Data - Used for Image I/O */
  void SetMetaDataHeader();

  const AxisNodeType & GetNodeUsedInCalculation(unsigned int idx) const
  {
    return m_NodesUsed[idx];
  }

  /** Generate the output image meta information. */
  virtual void GenerateOutputInformation() ITK_OVERRIDE;

  virtual void EnlargeOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  /** Get Large Value. This value is used to
      represent the concept of infinity for the time assigned to pixels that
      have not been visited. This value is set by default to half the
      max() of the pixel type used to represent the time-crossing map. */
  itkGetConstReferenceMacro( LargeValue, PixelType );

  OutputRegionType m_BufferedRegion;
  typedef typename LevelSetImageType::IndexType LevelSetIndexType;
  LevelSetIndexType m_StartIndex;
  LevelSetIndexType m_LastIndex;

  itkGetConstReferenceMacro( StartIndex, LevelSetIndexType );
  itkGetConstReferenceMacro( LastIndex, LevelSetIndexType );
private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DtiFastMarchingCostFilter);

  NodeContainerPointer m_AlivePoints;
  NodeContainerPointer m_TrialPoints;

  LevelSetPointer        m_OutputImage;
  TensorImagePointer     m_TensorImage;
  AnisotropyImagePointer m_AnisotropyImage;

  LabelImagePointer       m_LabelImage;
  OutputSpeedImagePointer m_OutputSpeedImage;
  EigenvectorImagePointer m_EigenvectorImage;

  double m_StoppingValue;
  float  m_AnisotropyWeight;

  bool                 m_CollectPoints;
  NodeContainerPointer m_ProcessedPoints;

  OutputRegionType    m_OutputRegion;
  OutputDirectionType m_OutputDirection;
  OutputSpacingType   m_OutputSpacing;
  OutputPointType     m_OutputOrigin;
  bool                m_OverrideOutputInformation;

  typename LevelSetImageType::PixelType m_LargeValue;
  AxisNodeType m_NodesUsed[dimension];

  /** Trial points are stored in a min-heap. This allow efficient access
   * to the trial point with minimum value which is the next grid point
   * the algorithm processes. */
  typedef std::vector<AxisNodeType>                                      HeapContainer;
  typedef std::greater<AxisNodeType>                                     NodeComparer;
  typedef std::priority_queue<AxisNodeType, HeapContainer, NodeComparer> HeapType;

  HeapType m_TrialHeap;

  double m_NormalizationFactor;
}; // end class
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDtiFastMarchingCostFilter.hxx"
#endif

#endif
