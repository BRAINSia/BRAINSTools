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
#include "BRAINSCutUtilities.h"

//
// read/warp image
//
WorkingImagePointer
SmoothImage( const WorkingImagePointer& image, const float GaussianValue)
{
  if( GaussianValue < 0 + FLOAT_TOLERANCE )
    {
    std::cout << "Gaussian value is less than tolerance. "
              << "No smoothing occurs at this time"
              << std::endl;
    return image;
    }
  /*std::cout<<"Smooth Image with Gaussian value of :: "
           << GaussianValue
           <<std::endl;*/
  using SmoothingFilterType = itk::SmoothingRecursiveGaussianImageFilter<WorkingImageType, WorkingImageType>;
  SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();

  smoothingFilter->SetInput( image);
  smoothingFilter->SetSigma( GaussianValue );

  smoothingFilter->Update();

  return smoothingFilter->GetOutput();
}

WorkingImagePointer ReadImageByFilename( const std::string  & filename )
{
  std::cout << "********************************************************" << std::endl;
  std::cout << "ReadImageByFilename::: " << filename << std::endl;
  std::cout << "********************************************************" << std::endl;

  WorkingImagePointer readInImage;

  ReadInImagePointer inputImage = itkUtil::ReadImage<ReadInImageType>(filename.c_str() );

  readInImage = itkUtil::ScaleAndCast<ReadInImageType,
                                      WorkingImageType>(inputImage,
                                                        0.0F,
                                                        1.0F );
  return readInImage;
}

/* inline functions */

DisplacementFieldType::Pointer GetDeformationField( std::string filename)
{
  const bool useTransform( filename.find(".mat") == std::string::npos ||
                           filename.find(".h5") == std::string::npos );

  if( useTransform )
    {
    std::cout << "Return null deformation. Use transformation instead." << std::endl;
    return nullptr;
    }
  using DeformationReaderType = itk::ImageFileReader<DisplacementFieldType>;
  DeformationReaderType::Pointer deformationReader = DeformationReaderType::New();
  deformationReader->SetFileName( filename );
  deformationReader->Update();

  return deformationReader->GetOutput();
}

itk::Transform<double, 3, 3>::Pointer GetGenericTransform( std::string filename)
{
  const bool useDeformation( filename.find(".mat") != std::string::npos &&
                             filename.find(".h5") != std::string::npos &&
                             filename.find(".hdf5") != std::string::npos &&
                             filename.find(".txt") != std::string::npos );

  if( useDeformation )
    {
    std::cout << "Return null transformation. Use deformation instead." << std::endl;
    return nullptr;
    }
  return itk::ReadTransformFromDisk( filename );
}
