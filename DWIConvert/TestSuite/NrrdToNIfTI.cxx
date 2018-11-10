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
#include "DWIConvertUtils.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

const int DIMENSION(4);
using PixelType = unsigned short;
using Image4DType = itk::Image<PixelType, DIMENSION>;
using Image3DType = itk::Image<PixelType, DIMENSION - 1>;
using VectorImageType = itk::VectorImage<PixelType, 3>;

void
ConvertVectorImageTo4DImage(const VectorImageType * img,
                            Image4DType::Pointer & img4D)
{
  VectorImageType::SizeType      size3D(img->GetLargestPossibleRegion().GetSize() );
  VectorImageType::DirectionType direction3D(img->GetDirection() );
  VectorImageType::SpacingType   spacing3D(img->GetSpacing() );
  VectorImageType::PointType     origin3D(img->GetOrigin() );

  Image4DType::SizeType size4D;

  size4D[0] = size3D[0];
  size4D[1] = size3D[1];
  size4D[2] = size3D[2];
  size4D[3] = 3; // img->GetNumberOfComponentsPerPixel();

  Image4DType::DirectionType direction4D;
  Image4DType::SpacingType   spacing4D;
  Image4DType::PointType     origin4D;
  for( unsigned i = 0; i < 3; ++i )
    {
    for( unsigned j = 0; j < 3; ++j )
      {
      direction4D[i][j] = direction3D[i][j];
      }
    direction4D[3][i] = 0.0;
    direction4D[i][3] = 0.0;
    spacing4D[i] = spacing3D[i];
    origin4D[i] = origin3D[i];
    }
  direction4D[3][3] = 1.0;
  spacing4D[3] = 1.0;
  origin4D[3] = 0.0;

  img4D = Image4DType::New();
  img4D->SetRegions(size4D);
  img4D->SetDirection(direction4D);
  img4D->SetSpacing(spacing4D);
  img4D->SetOrigin(origin4D);
  img4D->Allocate();

  Image4DType::IndexType     index4D;
  VectorImageType::IndexType vectorIndex;
  for( vectorIndex[2] = 0; vectorIndex[2] < static_cast<VectorImageType::IndexType::IndexValueType>( size3D[2] );
       ++vectorIndex[2] )
    {
    index4D[2] = vectorIndex[2];
    for( vectorIndex[1] = 0; vectorIndex[1] < static_cast<VectorImageType::IndexType::IndexValueType>( size3D[1] );
         ++vectorIndex[1] )
      {
      index4D[1] = vectorIndex[1];
      for( vectorIndex[0] = 0; vectorIndex[0] < static_cast<VectorImageType::IndexType::IndexValueType>( size3D[0] );
           ++vectorIndex[0] )
        {
        index4D[0] = vectorIndex[0];
        const VectorImageType::PixelType pixel = img->GetPixel(vectorIndex);
        for( index4D[3] = 0; index4D[3] < static_cast<VectorImageType::IndexType::IndexValueType>( size4D[3] );
             ++index4D[3] )
          {
          img4D->SetPixel(index4D, pixel[index4D[3]]);
          }
        }
      }
    }
}

int main(int argc, char *argv[])
{
  if( argc != 4 )
    {
    std::cout << "ERROR: Must provide 3 arguments" << std::endl;
    std::cout << "USAGE: " << argv[0] << "input output outputName3D" << std::endl;
    return EXIT_FAILURE;
    }
  VectorImageType::Pointer img;
  std::string              inputName(argv[1]), outputName(argv[2]),
  outputName3D(argv[3]);

  ReadVectorVolume<VectorImageType>(img, inputName, false);
  std::cout << "input directions" << std::endl
            << img->GetDirection() << std::endl;

  Image4DType::Pointer img4D;
  Image3DType::Pointer img3D;

  ConvertVectorImageTo4DImage(img, img4D);
  WriteVolume<Image4DType>(img4D, outputName);

  ReadScalarVolume<Image3DType>(img3D, inputName, false);
  WriteVolume<Image3DType>(img3D, outputName3D);

  ReadScalarVolume<Image4DType>(img4D, outputName, false);
  std::cout << "output directions" << std::endl
            << img4D->GetDirection() << std::endl;

  ReadScalarVolume<Image3DType>(img3D, outputName3D,false);
  std::cout << "output directions" << std::endl
            << img3D->GetDirection() << std::endl;

  return 0;
}
