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
#ifndef __itkPiecewiseRescaleQuadEdgeMeshFilter_h
#define __itkPiecewiseRescaleQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class PiecewiseRescaleQuadEdgeMeshFilter
 * \brief This filter rescales the scalar values of a mesh.
 *
 * The output mesh will have the same geometry as the input.
 *
 * It maps the input value that is in [min_in,cValue) to [min_out,cValue);
 *                                    (cValue,max_in] to (cValue,max_out].
 * and keeps cValue unchanged.
 *
 * \ingroup MeshFilters
 *
 */
template <class TInputMesh, class TOutputMesh>
class PiecewiseRescaleQuadEdgeMeshFilter :
  public QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh>
{
public:
  typedef PiecewiseRescaleQuadEdgeMeshFilter                        Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter<TInputMesh, TOutputMesh> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( PiecewiseRescaleQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter );

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro( Self );

  typedef TInputMesh                                             InputMeshType;
  typedef typename InputMeshType::PixelType                      InputPixelType;
  typedef typename InputMeshType::PointDataContainer             InputPointDataContainer;
  typedef typename InputMeshType::PointDataContainerConstPointer InputPointDataContainerConstPointer;

  typedef TOutputMesh                                        OutputMeshType;
  typedef typename OutputMeshType::PixelType                 OutputPixelType;
  typedef typename OutputMeshType::PointDataContainer        OutputPointDataContainer;
  typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;

  /** Set/Get the mesh that will be deformed. */
  void SetInputMesh( const InputMeshType * mesh );

  const InputMeshType * GetInputMesh( void ) const;

  /** Set/Get min value of the output mesh scalars. */
  itkSetMacro( OutputMinimum, OutputPixelType );
  itkGetMacro( OutputMinimum, OutputPixelType );

  /** Set/Get max value of the output mesh scalars. */
  itkSetMacro( OutputMaximum, OutputPixelType );
  itkGetMacro( OutputMaximum, OutputPixelType );

  /** Get min value of the input mesh scalars. */
  itkGetMacro( InputMinimum, InputPixelType );

  /** Get max value of the input mesh scalars. */
  itkGetMacro( InputMaximum, InputPixelType );

  /** Set/Get value that we want to keep unchanged. */
  itkSetMacro( cValue, InputPixelType );
  itkGetMacro( cValue, InputPixelType );

  /** Get scale to transform mesh scalars in [min_in,cValue). */
  itkGetMacro( Scale_a, double );

  /** Get scale to transform mesh scalars in (cValue,max_in]. */
  itkGetMacro( Scale_b, double );
protected:
  PiecewiseRescaleQuadEdgeMeshFilter();
  ~PiecewiseRescaleQuadEdgeMeshFilter();

  void GenerateData() ITK_OVERRIDE;

private:

  ITK_DISALLOW_COPY_AND_ASSIGN(PiecewiseRescaleQuadEdgeMeshFilter);

  OutputPixelType m_OutputMinimum;
  OutputPixelType m_OutputMaximum;

  InputPixelType m_InputMinimum;
  InputPixelType m_InputMaximum;

  InputPixelType m_cValue;

  double m_Scale_a;
  double m_Scale_b;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPiecewiseRescaleQuadEdgeMeshFilter.hxx"
#endif

#endif
