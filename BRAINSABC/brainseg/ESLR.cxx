#include "ESLRCLP.h"
#include "BRAINSThreadControl.h"
#include "ExtractSingleLargestRegion.h"

#include <itksys/SystemTools.hxx>

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);

  typedef itk::Image<unsigned char, 3>        ByteImageType;
  typedef itk::ImageFileReader<ByteImageType> ReaderType;

  ReaderType::Pointer myReader = ReaderType::New();
  myReader->SetFileName(inputVolume);
  myReader->Update();
  ByteImageType::Pointer myDirtyRegion = myReader->GetOutput();

  ByteImageType::Pointer myCleanRegion
    = ExtractSingleLargestRegion(low, high, openingSize, closingSize, safetySize, myDirtyRegion);

  if( preserveOutside == true )  // For values outside the specified range,
                                 // preserve those values.
    {
    std::cout << "PRESERVING OUTSIDE VALUES" << std::endl;
    itk::ImageRegionConstIterator<ByteImageType> dit( myDirtyRegion, myDirtyRegion->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ByteImageType>      cit( myCleanRegion, myCleanRegion->GetLargestPossibleRegion() );
    dit.GoToBegin();
    cit.GoToBegin();
    while( ( !cit.IsAtEnd() ) && ( !dit.IsAtEnd() ) )
      {
      if( ( dit.Get() < low  ) || ( dit.Get() > high ) )  // Outside of cleaning
                                                          // range, then
                                                          // preserve old value
        {
        cit.Set( dit.Get() );
        }
      ++cit;
      ++dit;
      }
    }

  typedef itk::ImageFileWriter<ByteImageType> OutputWriterType;
  OutputWriterType::Pointer writer = OutputWriterType::New();

  writer->SetInput(  myCleanRegion );
  writer->UseCompressionOn();
  const std::string fn = std::string(outputVolume);
  writer->SetFileName( fn.c_str() );
  writer->Update();

  return 0;
}
