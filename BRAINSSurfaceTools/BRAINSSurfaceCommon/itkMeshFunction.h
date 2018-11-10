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
#ifndef __itkMeshFunction_h
#define __itkMeshFunction_h

#include "itkFunctionBase.h"
#include "itkMesh.h"
#include "itkPoint.h"

namespace itk
{
/** \class MeshFunction
 * \brief Evaluates a function of a Mesh at specified position.
 *
 * MeshFunction is a baseclass for all objects that evaluates
 * a function of a mesh at a given point.
 *
 * This class is templated over the input mesh type, the type
 * of the function output and the coordinate representation type
 * (e.g. float or double).
 *
 * The input mesh is set via method SetInputMesh().
 * Methods Evaluate, EvaluateAtIndex and EvaluateAtContinuousIndex
 * respectively evaluates the function at an geometric point,
 * mesh index and continuous mesh index.
 *
 * \warning Mesh information is cached during SetInputMesh( mesh ). If the mesh
 * has changed one must call SetInputMesh( mesh ) again to update the cache to
 * the current values. The cached values are the ones of the KdTree point locator.
 *
 * \sa Point
 *
 * \ingroup MeshFunctions
 */
template <typename TInputMesh, typename TOutput>
class MeshFunction :
  public         FunctionBase<typename TInputMesh::PointType, TOutput>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MeshFunction);

  /** Dimension underlying input mesh. */
  static constexpr unsigned int MeshDimension = TInputMesh::PointDimension;

  /** Standard class type alias. */
  using Self = MeshFunction;
  using Superclass = FunctionBase<
      typename TInputMesh::PointType, TOutput>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshFunction, FunctionBase);

  /** InputMeshType type alias support. */
  using InputMeshType = TInputMesh;

  /** InputPixel type alias support */
  using InputPixelType = typename InputMeshType::PixelType;

  /** InputMeshPointer type alias support */
  using InputMeshConstPointer = typename InputMeshType::ConstPointer;

  /** OutputType type alias support. */
  using OutputType = TOutput;

  /** CoordRepType type alias support. */
  using CoordRepType = typename InputMeshType::CoordRepType;

  /** Point Type. */
  using PointType = typename InputMeshType::PointType;

  /** Set the input mesh.
   * \warning this method caches information.  If the Mesh has changed, user
   * must call SetInputMesh again to update cached values. */
  virtual void SetInputMesh( const InputMeshType * ptr );

  /** Get the input mesh. */
  const InputMeshType * GetInputMesh() const;

  /** Evaluate the function at specified Point position.
   * Subclasses must provide this method. */
  virtual TOutput Evaluate( constexpr PointType& point ) const override   = 0;

protected:
  MeshFunction();
  ~MeshFunction()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** Const pointer to the input image. */
  InputMeshConstPointer m_Mesh;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMeshFunction.hxx"
#endif

#endif
