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
#include "itkIO.h"
#include "GenericTransformImage.h"
#include "itkVersorRigid3DTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkAffineTransform.h"
#include "itkBSplineDeformableTransform.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_random.h"

#if !defined( _WIN32 )
#include <unistd.h>
#else
#include <process.h>
inline int getpid()
{
  return _getpid();
}

#endif

template <typename TTransform>
typename TTransform::Pointer
CreateTransform()
{
  typename TTransform::Pointer rval =
    TTransform::New();
  rval->SetIdentity();
  return rval;
}

int main(int argc, char * *argv)
{
  if( argc < 2 )
    {
    std::cerr << "Missing working directory argument" << std::endl;
    return EXIT_FAILURE;
    }

  using BSplineDeformableTransformType = itk::BSplineDeformableTransform<double, 3, 3>;
  using VersorRigid3DTransformType = itk::VersorRigid3DTransform<double>;
  using ScaleVersor3DTransformType = itk::ScaleVersor3DTransform<double>;
  using ScaleSkewVersor3DTransformType = itk::ScaleSkewVersor3DTransform<double>;
  VersorRigid3DTransformType::Pointer versorRigidTransform
    = CreateTransform<VersorRigid3DTransformType>();
  VersorRigid3DTransformType::InputPointType center;
  VersorRigid3DTransformType::AxisType       axis;

  center[0] = 0.25;  center[1] = -0.25;  center[2] = 0.333;
  versorRigidTransform->SetCenter(center);

  axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
  versorRigidTransform->SetRotation(axis, 1.5);

  std::string versorRigidName(argv[1]);
  versorRigidName += "/VersorRigidTransform.txt";
  itk::WriteTransformToDisk<double>(versorRigidTransform.GetPointer(), versorRigidName);

  ScaleVersor3DTransformType::Pointer scaleVersorTransform =
    CreateTransform<ScaleVersor3DTransformType>();
  ScaleVersor3DTransformType::OutputVectorType translation;

  translation[0] = 0.5; translation[1] = -0.6; translation[2] = 0.73;
  scaleVersorTransform->SetTranslation(translation);

  center[0] = -0.5; center[1] = 0.8; center[2] = 0.99;
  scaleVersorTransform->SetCenter(center);

  axis[0] = -1.0; axis[1] = 0.0; axis[2] = 0.0;
  scaleVersorTransform->SetRotation(axis, 0.5723);
  ScaleVersor3DTransformType::ScaleVectorType scale;

  scale[0] = .75; scale[1] = 1.1; scale[2] = 0.3333;
  scaleVersorTransform->SetScale(scale);

  std::string scaleVersorName(argv[1]);
  scaleVersorName += "/ScaleVersorTransform.txt";
  itk::WriteTransformToDisk<double>(scaleVersorTransform, scaleVersorName);

  ScaleSkewVersor3DTransformType::Pointer scaleSkewVersorTransform =
    CreateTransform<ScaleSkewVersor3DTransformType>();
  translation[0] = -0.5; translation[1] = 0.6; translation[2] = -0.73;
  scaleSkewVersorTransform->SetTranslation(translation);

  center[0] = 0.7; center[1] = -0.8; center[2] = 0.97;
  scaleSkewVersorTransform->SetCenter(center);

  axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
  scaleSkewVersorTransform->SetRotation(axis, 0.993);

  scale[0] = .33; scale[1] = .5; scale[2] = 0.666;
  scaleSkewVersorTransform->SetScale(scale);

  std::string scaleSkewVersorName(argv[1]);
  scaleSkewVersorName += "/ScaleSkewVersorTransform.txt";
  itk::WriteTransformToDisk<double>(scaleSkewVersorTransform, scaleSkewVersorName);

  using AffineTransformType = itk::AffineTransform<double, 3>;
  AffineTransformType::Pointer affineTransform =
    CreateTransform<AffineTransformType>();

  translation[0] = -1.5; translation[1] = -0.7; translation[2] = -0.02;
  affineTransform->Translate(translation);

  center[0] = -0.77; center[1] = 0.88; center[2] = -0.97;
  affineTransform->SetCenter(center);

  axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
  affineTransform->Rotate3D(axis, 0.993);

  scale[0] = .65; scale[1] = 1.3; scale[2] = 0.9;
  affineTransform->Scale(scale);

  affineTransform->Shear(0, 1, 0.25);
  affineTransform->Shear(1, 0, 0.75);
  affineTransform->Shear(1, 2, 0.35);

  std::string affineName(argv[1]);
  affineName += "/AffineTransform.txt";
  itk::WriteTransformToDisk<double>(affineTransform, affineName);

  using ImageType = itk::Image<signed short, 3>;
  ImageType::RegionType            region;
  ImageType::RegionType::SizeType  size;
  ImageType::RegionType::IndexType index;
  ImageType::SpacingType           spacing;
  ImageType::PointType             origin;

  size[0] = size[1] = size[2] = 10;
  region.SetSize(size);
  index[0] = index[1] = index[2] = 0;
  region.SetSize(size);
  region.SetIndex(index);
  spacing[0] = spacing[1] = spacing[2] = 2.0;
  origin[0] = -10; origin[1] = -10; origin[2] = 10;

  ImageType::Pointer testImage =
    itkUtil::AllocateImageFromRegionAndSpacing<ImageType>(region, spacing);

  vnl_random randgen;
  randgen.reseed( getpid() );
  for( itk::ImageRegionIterator<ImageType> it(testImage, testImage->GetLargestPossibleRegion() );
       !it.IsAtEnd(); ++it )
    {
    it.Set(static_cast<ImageType::PixelType>(randgen.lrand32(32767) ) );
    }

  std::string testImageName(argv[1]);
  testImageName += "/TransformConvertTestImage.nii.gz";
  itkUtil::WriteImage<ImageType>(testImage, testImageName);

  BSplineDeformableTransformType::Pointer bsplineTransform =
    CreateTransform<BSplineDeformableTransformType>();

  translation[0] = -1.0; translation[1] = 0.6; translation[2] = -0.5;
  affineTransform->Translate(translation);

  center[0] = 0.77; center[1] = -0.8; center[2] = 0.03;
  affineTransform->SetCenter(center);

  axis[0] = 0.0; axis[1] = -1.0; axis[2] = 0.0;
  affineTransform->Rotate3D(axis, 0.45);

  scale[0] = .8; scale[1] = .5; scale[2] = 0.3;
  affineTransform->Scale(scale);

  affineTransform->Shear(0, 1, 0.3);
  affineTransform->Shear(1, 0, 0.4);
  affineTransform->Shear(1, 2, 0.5);

  bsplineTransform->SetBulkTransform(affineTransform.GetPointer() );

  BSplineDeformableTransformType::PhysicalDimensionsType fixedPhysicalDimensions;
  BSplineDeformableTransformType::MeshSizeType           meshSize;
  for( unsigned int i = 0; i < 3; i++ )
    {
    fixedPhysicalDimensions[i] = spacing[i]
      * static_cast<double>(region.GetSize()[i] - 1);
    }
  meshSize.Fill( 5 - BSplineDeformableTransformType::SplineOrder );
  bsplineTransform->SetGridOrigin(origin);
  bsplineTransform->SetGridRegion(region);
  bsplineTransform->SetGridSpacing(spacing);

  BSplineDeformableTransformType::ParametersType parameters(bsplineTransform->GetNumberOfParameters() );
  parameters.Fill(0.0);
  bsplineTransform->SetParameters(parameters);

  std::string bsplineName(argv[1]);
  bsplineName += "/BSplineDeformableTransform.txt";
  itk::WriteTransformToDisk<double>(bsplineTransform, bsplineName);

  return EXIT_SUCCESS;
}
