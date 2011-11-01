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

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"

#include "extractNrrdVectorIndexCLP.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);

  bool debug = true;
  if( debug )
    {
    std::cout << "Input Image: " <<  inputVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Vector Index: " <<  vectorIndex << std::endl;
    std::cout << "Set Image Orientation: " <<  setImageOrientation << std::endl << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short                   PixelType;
  typedef itk::VectorImage<PixelType, 3> NrrdImageType;

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( inputVolume );
  reader->Update();

  if( ( vectorIndex < 0 )
      || ( vectorIndex >= static_cast<int>( ( reader->GetOutput() )->GetVectorLength() ) ) )
    {
    const int maxIndex = reader->GetOutput()->GetVectorLength() - 1;
    std::cerr << "Invalid vector image index (" << vectorIndex << "), valid indexes are 0-" << maxIndex << std::endl;
    return 1;
    }

  typedef itk::Image<PixelType, 3>                                                IndexImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<NrrdImageType, IndexImageType> VectorSelectFilterType;
  typedef VectorSelectFilterType::Pointer                                         VectorSelectFilterPointer;

  VectorSelectFilterPointer SelectIndexImageFilter = VectorSelectFilterType::New();
  SelectIndexImageFilter->SetIndex( vectorIndex );
  SelectIndexImageFilter->SetInput( reader->GetOutput() );
  try
    {
    SelectIndexImageFilter->Update();
    }
  catch( itk::ExceptionObject e )
    {
    std::cout << e << std::endl;
    }

  /* Hack Required for Certain Output Image Types */
  itk::MetaDataDictionary       meta;
  IndexImageType::Pointer       indexImage = SelectIndexImageFilter->GetOutput();
  IndexImageType::DirectionType fixImageDir = indexImage->GetDirection();
#if ( ITK_VERSION_MAJOR < 4  )
#define EncapsulateMD(image, flag)                                       \
    {                                                                   \
    itk::EncapsulateMetaData                                            \
    <itk::SpatialOrientation::ValidCoordinateOrientationFlags>      \
      (meta, itk::ITK_CoordinateOrientation, flag);                   \
    image->SetMetaDataDictionary(meta);                                 \
    }
#else
#define EncapsulateMD(image, flag) {}
#endif
  if( setImageOrientation == "Axial"  ||  setImageOrientation == "AXIAL"  ||  setImageOrientation == "axial" )
    {
    fixImageDir.Fill(0);
    fixImageDir[0][0] = 1; fixImageDir[1][1] = 1; fixImageDir[2][2] = 1;
    indexImage->SetDirection( fixImageDir );
    EncapsulateMD(indexImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI);
    }
  else if( setImageOrientation == "Coronal"  ||  setImageOrientation == "CORONAL"  ||  setImageOrientation ==
           "coronal" )
    {
    fixImageDir.Fill(0);
    fixImageDir[0][0] = 1; fixImageDir[1][2] = 1; fixImageDir[2][1] = 1;
    indexImage->SetDirection( fixImageDir );
    EncapsulateMD(indexImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP);
    }
  else if( setImageOrientation == "Sagittal"  ||  setImageOrientation == "SAGITTAL"  ||  setImageOrientation ==
           "sagittal" )
    {
    fixImageDir.Fill(0);
    fixImageDir[0][2] = 1; fixImageDir[1][1] = 1; fixImageDir[2][0] = 1;
    indexImage->SetDirection( fixImageDir );
    EncapsulateMD(indexImage, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR);
    }
#undef EncapsulateMD // (image, flag)
  // else, leave it AsAcquired.

  typedef itk::ImageFileWriter<IndexImageType> WriterType;
  WriterType::Pointer imageWriter = WriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetInput( indexImage );
  imageWriter->SetFileName( outputVolume );

  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject e )
    {
    std::cout << e << std::endl;
    }

  // IndexImageType::Pointer indexImage = SelectIndexImageFilter->GetOutput();
  // std::cout << "Index Image: "<< indexImage  << std::endl;

  return EXIT_SUCCESS;
}
