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
#ifndef __itkHistogramMatchingQuadEdgeMeshFilter_h
#define __itkHistogramMatchingQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkHistogram.h"
#include "vnl/vnl_matrix.h"

namespace itk
{
/**
 * \class HistogramMatchingQuadEdgeMeshFilter
 * \brief brief Normalize the scalar values between two meshes by histogram
 * matching.
 *
 * HistogramMatchingQuadEdgeMeshFilter normalizes the scalar values of a source
 * mesh based on the scalar values of a reference mesh.
 * This filter uses a histogram matching technique where the histograms of the
 * two meshes are matched only at a specified number of quantile values.
 *
 * This filter was originally designed to normalized scalar values of two
 * brain surfaces. No mean scalar value is used to cut out the background
 * as the itkHistogramMatchingImageFilter does.
 *
 * The source mesh can be set via either SetInput() or SetSourceMesh().
 * The reference image can be set via SetReferenceMesh().
 *
 * SetNumberOfHistogramLevels() sets the number of bins used when
 * creating histograms of the source and reference meshes.
 * SetNumberOfMatchPoints() governs the number of quantile values to be
 * matched.
 *
 * This filter assumes that both the source and reference are of the same
 * type and have the same number of scalars, and that the input and output
 * mesh type have the same geometry and have scalar pixel types.
 *
 * \ingroup MeshFilters
 *
 */
/* THistogramMeasurement -- The precision level for which to do HistogramMeasurmenets */
template <class TInputMesh, class TOutputMesh, class THistogramMeasurement = typename TInputMesh::PixelType>
class HistogramMatchingQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef HistogramMatchingQuadEdgeMeshFilter Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>                     Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( HistogramMatchingQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef typename  Superclass::InputMeshType                     InputMeshType;
  typedef typename  InputMeshType::PixelType                      InputPixelType;
  typedef typename  InputMeshType::ConstPointer                   InputMeshConstPointer;
  typedef typename  InputMeshType::PointDataContainer             InputPointDataContainer;
  typedef typename  InputMeshType::PointDataContainerConstPointer InputPointDataContainerConstPointer;

  typedef typename  Superclass::OutputMeshType           OutputMeshType;
  typedef typename  OutputMeshType::PixelType            OutputPixelType;
  typedef typename  OutputMeshType::Pointer              OutputMeshPointer;
  typedef typename  Superclass::OutputPointDataContainer OutputPointDataContainer;
  typedef typename  OutputPointDataContainer::Pointer    OutputPointDataContainerPointer;

  /** Histogram related typedefs. */
  typedef Statistics::Histogram<THistogramMeasurement> HistogramType;

  typedef typename HistogramType::Pointer HistogramPointer;

  /** Set/Get the source mesh. */
  void SetSourceMesh( const InputMeshType * source );

  const InputMeshType * GetSourceMesh( void ) const;

  /** Set/Get the reference mesh. */
  void SetReferenceMesh( const InputMeshType * reference );

  const InputMeshType * GetReferenceMesh( void ) const;

  /** Set/Get the number of histogram levels used. */
  itkSetMacro( NumberOfHistogramLevels, unsigned long );
  itkGetConstMacro( NumberOfHistogramLevels, unsigned long );

  /** Set/Get the number of match points used. */
  itkSetMacro( NumberOfMatchPoints, unsigned long );
  itkGetConstMacro( NumberOfMatchPoints, unsigned long );

  /** Methods to get the histograms of the source, reference, and
   * output. Objects are only valid after Update() has been called
   * on this filter. */
  itkGetConstObjectMacro(SourceHistogram, HistogramType);
  itkGetConstObjectMacro(ReferenceHistogram, HistogramType);
  itkGetConstObjectMacro(OutputHistogram, HistogramType);
protected:
  HistogramMatchingQuadEdgeMeshFilter();
  ~HistogramMatchingQuadEdgeMeshFilter();

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  void BeforeTransform();

  void Transform();

  void GenerateData() ITK_OVERRIDE;

  /** Compute min, max and mean of an image. */
  void ComputeMinMax( const InputMeshType * mesh, THistogramMeasurement& minValue, THistogramMeasurement& maxValue );

  /** Construct a histogram from a mesh. */
  void ConstructHistogram( const InputMeshType * mesh, HistogramType * histogram, const THistogramMeasurement minValue,
                           const THistogramMeasurement maxValue );

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(HistogramMatchingQuadEdgeMeshFilter);

  unsigned long m_NumberOfHistogramLevels;
  unsigned long m_NumberOfMatchPoints;

  THistogramMeasurement m_SourceMinValue;
  THistogramMeasurement m_SourceMaxValue;
  THistogramMeasurement m_ReferenceMinValue;
  THistogramMeasurement m_ReferenceMaxValue;
  THistogramMeasurement m_OutputMinValue;
  THistogramMeasurement m_OutputMaxValue;

  HistogramPointer m_SourceHistogram;
  HistogramPointer m_ReferenceHistogram;
  HistogramPointer m_OutputHistogram;

  typedef vnl_matrix<double> TableType;
  TableType m_QuantileTable;

  typedef vnl_vector<double> GradientArrayType;
  GradientArrayType m_Gradients;
  double            m_LowerGradient;
  double            m_UpperGradient;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHistogramMatchingQuadEdgeMeshFilter.hxx"
#endif

#endif
