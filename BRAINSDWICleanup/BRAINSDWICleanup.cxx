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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "DWIMetaDataDictionaryValidator.h"
#include "BRAINSDWICleanupCLP.h"
#include <list>
#include <vector>
#include <list>

typedef signed short                        PixelType;
typedef itk::VectorImage<PixelType, 3>      NrrdImageType;
typedef itk::Image<PixelType, 3>            SingleComponentImageType;
typedef itk::ImageFileReader<NrrdImageType,
                             itk::DefaultConvertPixelTraits<PixelType> > ReaderType;

typedef itk::ImageFileWriter<NrrdImageType> WriterType;

typedef std::vector<std::string> GradStringVector;

template <typename TImage>
TImage *AllocVecImage(const itk::ImageBase<TImage::ImageDimension> *templateImage, unsigned long vecSize)
{
  typename TImage::Pointer p1 = TImage::New();
  p1->CopyInformation(templateImage);
  p1->SetRegions(templateImage->GetLargestPossibleRegion());
  p1->SetNumberOfComponentsPerPixel(vecSize);
  p1->Allocate();
  p1.GetPointer()->Register();
  return p1.GetPointer();
}

class inBadList
{
public:
  inBadList(const std::vector<int> &list) : valueList(list)
    {
    }
  bool operator() (const int &value)
    {
      for(unsigned int i = 0; i < valueList.size(); ++i)
        {
        if(valueList[i] == value)
          {
          return true;
          }
        }
      return false;
    }
private:
  const std::vector<int> &valueList;
};

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  if(badGradients.size() == 0)
    {
    std::cerr << "Missing bad gradient list" << std::endl;
    return 1;
    }
  if(inputVolume.empty())
    {
    std::cerr << "Missing input NRRD file name" << std::endl;
    return 1;
    }
  if(outputVolume.empty())
    {
    std::cerr << "Missing output NRRD file name" << std::endl;
    return 1;
    }
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( inputVolume );
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }
  std::cout << "Read Input Image..." << std::endl;
  NrrdImageType::Pointer inImage = imageReader->GetOutput();

  unsigned int numInputGradients = inImage->GetNumberOfComponentsPerPixel();
  unsigned int newGradientCount = numInputGradients - badGradients.size();

  //
  // make an index list containing the list of gradients/volumes to keep.
  std::list<int> keepIndices;
  for(unsigned int i = 0; i < numInputGradients; ++i)
    {
    keepIndices.push_back(i); // add all indicies
    }
  // remove the indices in badGradients
  inBadList pred(badGradients);
  keepIndices.remove_if(pred);

  if( keepIndices.size() != newGradientCount )
    {
    std::cerr << "ERROR: The number of gradients to be kept does not match the number of output components!" << std::endl
              << "*** Please check input list for bad gradients." << std::endl;
    return 1;
    }

  // create output volume
  NrrdImageType::Pointer outImage = AllocVecImage<NrrdImageType>(inImage,newGradientCount);

  // copy input vector pixels to output, skipping the bad gradients.
  itk::ImageRegionConstIterator<NrrdImageType> inIt(inImage,inImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<NrrdImageType> outIt(outImage,outImage->GetLargestPossibleRegion());
  for(; !inIt.IsAtEnd(); ++inIt, ++outIt)
    {
    NrrdImageType::PixelType inPix = inIt.Get();
    NrrdImageType::PixelType outpix(newGradientCount);
    std::list<int>::iterator keepIt = keepIndices.begin();
    for(unsigned int i = 0; i < newGradientCount; ++i, ++keepIt)
      {
      outpix[i] = inPix[*(keepIt)];
      }
    outIt.Set(outpix);
    }

  // deal with gradients in meta data
  DWIMetaDataDictionaryValidator    nrrdMetaDataValidator;
  nrrdMetaDataValidator.SetMetaDataDictionary( inImage->GetMetaDataDictionary() );

  // Get gradient table and update the gradient vectors based on keepIndices
  DWIMetaDataDictionaryValidator::GradientTableType inputGradTable = nrrdMetaDataValidator.GetGradientTable();

  // Now delete the gradient table to fill with new gradient values
  nrrdMetaDataValidator.DeleteGradientTable();

  DWIMetaDataDictionaryValidator::GradientTableType outputGradTable( newGradientCount );

  // add the good gradients to the outputGradTable
  std::list<int>::iterator keepIt = keepIndices.begin();
  for(unsigned int i = 0; i < keepIndices.size(); ++i, ++keepIt)
    {
    outputGradTable[i] = inputGradTable[*keepIt];
    }
  nrrdMetaDataValidator.SetGradientTable( outputGradTable );

  outImage->SetMetaDataDictionary( nrrdMetaDataValidator.GetMetaDataDictionary() );

  std::cout << "Write Output Image..." << std::endl;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( outImage );
  nrrdWriter->SetFileName( outputVolume );
  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    return 1;
    }
  return 0;
}
