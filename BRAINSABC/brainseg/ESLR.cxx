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
#include "ESLRCLP.h"
#include "BRAINSThreadControl.h"
#include "ExtractSingleLargestRegion.h"

#include <itksys/SystemTools.hxx>

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <BRAINSCommonLib.h>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );

  using ByteImageType = itk::Image< unsigned char, 3 >;
  using ReaderType = itk::ImageFileReader< ByteImageType >;

  ReaderType::Pointer myReader = ReaderType::New();
  myReader->SetFileName( inputVolume );
  myReader->Update();
  ByteImageType::Pointer myDirtyRegion = myReader->GetOutput();

  ByteImageType::Pointer myCleanRegion =
    ExtractSingleLargestRegion( low, high, openingSize, closingSize, safetySize, myDirtyRegion );

  if ( static_cast< int >( preserveOutside ) == true ) // For values outside the specified range,
                                                       // preserve those values.
  {
    std::cout << "PRESERVING OUTSIDE VALUES" << std::endl;
    itk::ImageRegionConstIterator< ByteImageType > dit( myDirtyRegion, myDirtyRegion->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< ByteImageType >      cit( myCleanRegion, myCleanRegion->GetLargestPossibleRegion() );
    dit.GoToBegin();
    cit.GoToBegin();
    while ( ( !cit.IsAtEnd() ) && ( !dit.IsAtEnd() ) )
    {
      if ( ( dit.Get() < low ) || ( dit.Get() > high ) ) // Outside of cleaning
                                                         // range, then
                                                         // preserve old value
      {
        cit.Set( dit.Get() );
      }
      ++cit;
      ++dit;
    }
  }

  using OutputWriterType = itk::ImageFileWriter< ByteImageType >;
  OutputWriterType::Pointer writer = OutputWriterType::New();

  writer->SetInput( myCleanRegion );
  writer->UseCompressionOn();
  const std::string fn = std::string( outputVolume );
  writer->SetFileName( fn.c_str() );
  writer->Update();

  return 0;
}
