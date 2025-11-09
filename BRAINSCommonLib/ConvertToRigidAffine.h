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
#ifndef __ConvertToRigidAffine_h
#define __ConvertToRigidAffine_h
#include <itkMatrix.h>
#include <itkVector.h>
#include <itkAffineTransform.h>
#include <itkVersorRigid3DTransform.h>
#include "itkScaleVersor3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkMacro.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>
#include <vcl_compiler.h>
#include <iostream>

// INFO:  Need to make return types an input template type.
namespace AssignRigid
{
using AffineTransformType = itk::AffineTransform<double, 3>;
using AffineTransformPointer = AffineTransformType::Pointer;

using VnlTransformMatrixType44 = vnl_matrix_fixed<double, 4, 4>;

using Matrix3D = itk::Matrix<double, 3, 3>;
using VersorType = itk::Versor<double>;

using MatrixType = AffineTransformType::MatrixType;
using PointType = AffineTransformType::InputPointType;
using VectorType = AffineTransformType::OutputVectorType;

using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
using VersorRigid3DTransformPointer = VersorRigid3DTransformType::Pointer;
using VersorRigid3DParametersType = VersorRigid3DTransformType::ParametersType;

using ScaleVersor3DTransformType = itk::ScaleVersor3DTransform<double>;
using ScaleVersor3DTransformPointer = ScaleVersor3DTransformType::Pointer;
using ScaleVersor3DParametersType = ScaleVersor3DTransformType::ParametersType;

using ScaleSkewVersor3DTransformType = itk::ScaleSkewVersor3DTransform<double>;
using ScaleSkewVersor3DTransformPointer = ScaleSkewVersor3DTransformType::Pointer;
using ScaleSkewVersor3DParametersType = ScaleSkewVersor3DTransformType::ParametersType;

using Similarity3DTransformType = itk::Similarity3DTransform<double>;
using Similarity3DTransformPointer = Similarity3DTransformType::Pointer;
using Similarity3DParametersType = Similarity3DTransformType::ParametersType;

/**
 * AffineTransformPointer  :=  AffineTransformPointer
 */
inline void
AssignConvertedTransform(AffineTransformPointer & result, const AffineTransformType::ConstPointer & affine)
{
  if (result.IsNotNull())
  {
    result->SetParameters(affine->GetParameters());
    result->SetFixedParameters(affine->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, while assigning "
              << "AffineTransformPointer := AffineTransformPointer." << std::endl;
    throw;
  }
}

/**
 * AffineTransformPointer  :=  VnlTransformMatrixType44
 */
inline void
AssignConvertedTransform(AffineTransformPointer & result, const VnlTransformMatrixType44 & matrix)
{
  if (result.IsNotNull())
  {
    MatrixType rotator; // can't do = conversion.
    rotator.operator=(matrix.extract(3, 3, 0, 0));

    VectorType offset;
    for (unsigned int i = 0; i < 3; ++i)
    {
      offset[i] = matrix.get(i, 3);
    }
    itk::Point<double, 3> ZeroCenter;
    ZeroCenter.Fill(0.0);
    result->SetIdentity();
    result->SetCenter(ZeroCenter); // Assume that rotation is about 0.0
    result->SetMatrix(rotator);
    result->SetOffset(offset); // It is offset in this case, and not
    // Translation.
  }
  else
  {
    std::cout << "Error missing Pointer data, while assigning "
              << "AffineTransformPointer := VnlTransformMatrixType44." << std::endl;
    throw;
  }
}

/**
 * VnlTransformMatrixType44  :=  AffineTransformPointer
 */
inline void
AssignConvertedTransform(VnlTransformMatrixType44 & result, const AffineTransformType::ConstPointer & affine)
{
  if (affine.IsNotNull())
  {
    MatrixType rotator = affine->GetMatrix();
    VectorType offset = affine->GetOffset(); // This needs to be offst
                                             // in
    // this case, and not
    // Translation.
    result.update(rotator.GetVnlMatrix().as_matrix(), 0, 0);
    for (unsigned int i = 0; i < 3; ++i)
    {
      result.put(i, 3, offset[i]);
      result.put(3, i, 0.0);
    }
    result.put(3, 3, 1.0);
  }
  else
  {
    std::cout << "Error missing Pointer data, while assigning"
              << " VnlTransformMatrixType44 := AffineTransformPointer." << std::endl;
    throw;
  }
}

/**
 * AffineTransformPointer  :=  ScaleSkewVersor3DTransformPointer
 */
inline void
AssignConvertedTransform(AffineTransformPointer & result, const ScaleSkewVersor3DTransformType::ConstPointer & scale)
{
  if (result.IsNotNull() && scale.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(scale->GetCenter());
    result->SetMatrix(scale->GetMatrix());
    result->SetTranslation(scale->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning AffineTransformPointer"
              << " := ScaleSkewVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * ScaleSkewVersor3DTransformPointer  :=  ScaleSkewVersor3DTransformPointer
 */
inline void
AssignConvertedTransform(ScaleSkewVersor3DTransformPointer &                  result,
                         const ScaleSkewVersor3DTransformType::ConstPointer & scale)
{
  if (result.IsNotNull() && scale.IsNotNull())
  {
    result->SetParameters(scale->GetParameters());
    result->SetFixedParameters(scale->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning AffineTransformPointer "
              << ":= ScaleSkewVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * AffineTransformPointer  :=  ScaleVersor3DTransformPointer
 */

inline void
AssignConvertedTransform(AffineTransformPointer & result, const ScaleVersor3DTransformType::ConstPointer & scale)
{
  if (result.IsNotNull() && scale.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(scale->GetCenter());
    result->SetMatrix(scale->GetMatrix()); // NOTE:  This matrix has
                                           // both
    // rotation ans scale components.
    result->SetTranslation(scale->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning AffineTransformPointer"
              << " := ScaleVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * ScaleVersor3DTransformPointer  :=  ScaleVersor3DTransformPointer
 */

inline void
AssignConvertedTransform(ScaleVersor3DTransformPointer & result, const ScaleVersor3DTransformType::ConstPointer & scale)
{
  if (result.IsNotNull() && scale.IsNotNull())
  {
    result->SetParameters(scale->GetParameters());
    result->SetFixedParameters(scale->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning ScaleVersor3DTransform"
              << " := ScaleVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * AffineTransformPointer  :=  VersorRigid3DTransformPointer
 */
inline void
AssignConvertedTransform(AffineTransformPointer &                         result,
                         const VersorRigid3DTransformType::ConstPointer & versorTransform)
{
  if (result.IsNotNull() && versorTransform.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(versorTransform->GetCenter());
    result->SetMatrix(versorTransform->GetMatrix()); // We MUST
                                                     // SetMatrix
    // before the
    // SetOffset -- not
    // after!
    result->SetTranslation(versorTransform->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning AffineTransformPointer"
              << " := VersorRigid3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * *VersorRigid3DTransformPointer  :=  VersorRigid3DTransformPointer
 */

inline void
AssignConvertedTransform(VersorRigid3DTransformPointer &                  result,
                         const VersorRigid3DTransformType::ConstPointer & versorRigid)
{
  if (result.IsNotNull() && versorRigid.IsNotNull())
  {
    result->SetParameters(versorRigid->GetParameters());
    result->SetFixedParameters(versorRigid->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " VersorRigid3DTransformPointer := VersorRigid3DTTransformPointer." << std::endl;
    throw;
  }
}

/**
 * ScaleSkewVersor3DTransformPointer  :=  ScaleVersor3DTransformPointer
 */
inline void
AssignConvertedTransform(ScaleSkewVersor3DTransformPointer &              result,
                         const ScaleVersor3DTransformType::ConstPointer & scale)
{
  if (result.IsNotNull() && scale.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(scale->GetCenter());
    result->SetRotation(scale->GetVersor());
    result->SetScale(scale->GetScale());
    result->SetTranslation(scale->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " ScaleSkewVersor3DTransformPointer := ScaleVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * ScaleSkewVersor3DTransformPointer  :=  VersorRigid3DTransformPointer
 */
inline void
AssignConvertedTransform(ScaleSkewVersor3DTransformPointer &              result,
                         const VersorRigid3DTransformType::ConstPointer & versorRigid)
{
  if (result.IsNotNull() && versorRigid.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(versorRigid->GetCenter());
    result->SetRotation(versorRigid->GetVersor());
    result->SetTranslation(versorRigid->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " ScaleSkewVersor3DTransformPointer := VersorRigid3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * *VersorRigid3DTransformPointer  :=  VersorRigid3DTransformPointer
 */

inline void
AssignConvertedTransform(Similarity3DTransformPointer &                  result,
                         const Similarity3DTransformType::ConstPointer & similarity3D)
{
  if (result.IsNotNull() && similarity3D.IsNotNull())
  {
    result->SetParameters(similarity3D->GetParameters());
    result->SetFixedParameters(similarity3D->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " VersorRigid3DTransformPointer := VersorRigid3DTTransformPointer." << std::endl;
    throw;
  }
}

/**
 * ScaleVersor3DTransformPointer  :=  VersorRigid3DTransformPointer
 */
inline void
AssignConvertedTransform(ScaleVersor3DTransformPointer &                  result,
                         const VersorRigid3DTransformType::ConstPointer & versorRigid)
{
  if (result.IsNotNull() && versorRigid.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(versorRigid->GetCenter());
    result->SetRotation(versorRigid->GetVersor());
    result->SetTranslation(versorRigid->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " ScaleVersor3DTransformPointer := VersorRigid3DTransformPointer." << std::endl;
    throw;
  }
}

inline void
ExtractVersorRigid3DTransform(VersorRigid3DTransformPointer &                  result,
                              const ScaleVersor3DTransformType::ConstPointer & scaleVersorRigid)
{
  if (result.IsNotNull() && scaleVersorRigid.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(scaleVersorRigid->GetCenter());
    result->SetRotation(scaleVersorRigid->GetVersor());
    result->SetTranslation(scaleVersorRigid->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " VersorRigid3DTransformPointer := ScaleVersor3DTransformPointer." << std::endl;
    throw;
  }
}

inline void
ExtractVersorRigid3DTransform(VersorRigid3DTransformPointer &                      result,
                              const ScaleSkewVersor3DTransformType::ConstPointer & scaleSkewVersorRigid)
{
  if (result.IsNotNull() && scaleSkewVersorRigid.IsNotNull())
  {
    result->SetIdentity();
    result->SetCenter(scaleSkewVersorRigid->GetCenter());
    result->SetRotation(scaleSkewVersorRigid->GetVersor());
    result->SetTranslation(scaleSkewVersorRigid->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " VersorRigid3DTransformPointer := ScaleSkewVersor3DTransformPointer." << std::endl;
    throw;
  }
}

inline void
ExtractVersorRigid3DTransform(VersorRigid3DTransformPointer &                  result,
                              const VersorRigid3DTransformType::ConstPointer & versorRigid)
{
  if (result.IsNotNull() && versorRigid.IsNotNull())
  {
    result->SetParameters(versorRigid->GetParameters());
    result->SetFixedParameters(versorRigid->GetFixedParameters());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning"
              << " VersorRigid3DTransformPointer := ScaleVersor3DTransformPointer." << std::endl;
    throw;
  }
}

/**
 * VersorRigid3DTransformPointer  :=  AffineTransformPointer
 */

/**
 * Utility function in which we claim the singular-value decomposition (svd)
 * gives the orthogonalization of a matrix in its U component.  However we
 * must clip out the null subspace, if any.
 */
inline Matrix3D
orthogonalize(const Matrix3D & rotator)
{
  vnl_svd<double>                             decomposition(rotator.GetVnlMatrix().as_matrix(), -1E-6);
  vnl_diag_matrix<vnl_svd<double>::singval_t> Winverse(decomposition.Winverse());

  vnl_matrix<double> W(3, 3);
  W.fill(static_cast<double>(0));
  for (unsigned int i = 0; i < 3; ++i)
  {
    if (decomposition.Winverse()(i, i) != 0.0)
    {
      W(i, i) = 1.0;
    }
  }

  vnl_matrix<double> result(decomposition.U() * W * decomposition.V().conjugate_transpose());

  //    std::cout << " svd Orthonormalized Rotation: " << std::endl
  //      << result << std::endl;
  Matrix3D Orthog;
  Orthog.operator=(result);

  return Orthog;
}

inline void
ExtractVersorRigid3DTransform(VersorRigid3DTransformPointer & result, const AffineTransformType::ConstPointer & affine)
{
  if (result.IsNotNull() && affine.IsNotNull())
  {
    Matrix3D   NonOrthog = affine->GetMatrix();
    Matrix3D   Orthog(orthogonalize(NonOrthog));
    MatrixType rotator;
    rotator.operator=(Orthog);

    VersorType versor;
    versor.Set(rotator); //    --> controversial!  Is rotator
                         // orthogonal as
    // required?
    // versor.Normalize();

    result->SetIdentity();
    result->SetCenter(affine->GetCenter());
    result->SetRotation(versor);
    result->SetTranslation(affine->GetTranslation());
  }
  else
  {
    std::cout << "Error missing Pointer data, assigning "
              << "VersorRigid3DTransformPointer := AffineTransformPointer." << std::endl;
    throw;
  }
}
} // namespace AssignRigid
#endif // __RigidAffine_h
