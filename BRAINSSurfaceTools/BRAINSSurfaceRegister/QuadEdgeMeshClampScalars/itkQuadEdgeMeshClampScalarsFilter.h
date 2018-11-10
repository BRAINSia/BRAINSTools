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
#ifndef __itkQuadEdgeMeshClampScalarsFilter_h
#define __itkQuadEdgeMeshClampScalarsFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class QuadEdgeMeshClampScalarsFilter
 * \brief This filter clamp scalar values on input mesh at either end (min or max).
 *
 *
 * This filter takes input and manipulate the scalar values.
 * If ClampMin is true, OutputMinimum has to be given by the user,
 * and the input scalar values lower than OutputMinimum is set to be OutputMinimum.
 * If ClampMax is true, OutputMaximum has to be given by the user,
 * and the input scalar values higher than OutputMaximum is set to be OutputMaximum.
 * \ingroup MeshFilters
 *
 */
template <typename TInputMesh, typename TOutputMesh>
class QuadEdgeMeshClampScalarsFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(QuadEdgeMeshClampScalarsFilter);

  using Self = QuadEdgeMeshClampScalarsFilter;
  using Superclass = QuadEdgeMeshToQuadEdgeMeshFilter<
      TInputMesh, TOutputMesh>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( QuadEdgeMeshClampScalarsFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  using OutputMeshType = typename  Superclass::OutputMeshType;
  using OutputMeshPointer = typename  OutputMeshType::Pointer;
  using OutputPixelType = typename  OutputMeshType::PixelType;
  using OutputPointDataContainer = typename  Superclass::OutputPointDataContainer;
  using OutputPointDataContainerPointer = typename  OutputPointDataContainer::Pointer;
  using OutputPointDataContainerIterator = typename  OutputPointDataContainer::Iterator;

  using InputMeshType = typename  Superclass::InputMeshType;
  using InputMeshPointer = typename  InputMeshType::Pointer;
  using InputPixelType = typename  InputMeshType::PixelType;
  using InputPointDataContainer = typename  Superclass::InputPointDataContainer;
  using InputPointDataContainerConstPointer = typename  InputPointDataContainer::ConstPointer;
  using InputPointDataContainerConstIterator = typename  InputPointDataContainer::ConstIterator;

  /** Set/Get ClampMin. */
  itkSetMacro( ClampMin, bool );
  itkGetMacro( ClampMin, bool );
  itkBooleanMacro( ClampMin );

  /** Set/Get ClampMax. */
  itkSetMacro( ClampMax, bool );
  itkGetMacro( ClampMax, bool );
  itkBooleanMacro( ClampMax );

  /** Set/Get OutputMinimum. */
  itkSetMacro( OutputMinimum, OutputPixelType );
  itkGetMacro( OutputMinimum, OutputPixelType );

  /** Set/Get OutputMaximum. */
  itkSetMacro( OutputMaximum, OutputPixelType );
  itkGetMacro( OutputMaximum, OutputPixelType );
protected:
  QuadEdgeMeshClampScalarsFilter();
  ~QuadEdgeMeshClampScalarsFilter();

  void GenerateData() override;

private:

  bool m_ClampMin;
  bool m_ClampMax;

  OutputPixelType m_OutputMinimum;
  OutputPixelType m_OutputMaximum;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadEdgeMeshClampScalarsFilter.hxx"
#endif

#endif
