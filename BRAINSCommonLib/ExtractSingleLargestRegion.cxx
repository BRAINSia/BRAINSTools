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
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"

#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include "itkMultiplyImageFilter.h"
#include "ExtractSingleLargestRegion.h"

itk::Image<unsigned char, 3>::Pointer
ExtractSingleLargestRegionFromMask(itk::Image<unsigned char, 3>::Pointer Mask,
                                   const int openingSize, const int closingSize, const int safetySize,
                                   itk::Image<unsigned char, 3>::Pointer inputLabelImage)
{
  typedef itk::Image<unsigned char, 3>                        ByteImageType;
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
  typedef
    itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
                                 StructElementType> DilateType;
  typedef itk::BinaryErodeImageFilter<ByteImageType, ByteImageType,
                                      StructElementType> ErodeType;
  typedef itk::ConnectedComponentImageFilter<ByteImageType, itk::Image<unsigned int, 3> > FilterType;

  ByteImageType::Pointer errodedImage = Mask;
  if( openingSize > 0 )
    {
    // Binary Erode
    StructElementType openStruct;
    openStruct.SetRadius(openingSize);
    openStruct.CreateStructuringElement();
    ErodeType::Pointer erode = ErodeType::New();
    erode->SetErodeValue(1);
    erode->SetKernel(openStruct);
    erode->SetInput(Mask);
    erode->Update();
    errodedImage = erode->GetOutput();
    }

  FilterType::Pointer labelConnectedComponentsFilter = FilterType::New();
  //  SimpleFilterWatcher watcher(labelConnectedComponentsFilter);
  //  watcher.QuietOn();
  labelConnectedComponentsFilter->SetInput( errodedImage );

  typedef itk::RelabelComponentImageFilter<itk::Image<unsigned int, 3>, ByteImageType> RelabelType;
  RelabelType::Pointer relabel = RelabelType::New();
  relabel->SetInput( labelConnectedComponentsFilter->GetOutput() );
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  typedef itk::BinaryThresholdImageFilter<ByteImageType, ByteImageType> BinaryThresholdFilter;
  BinaryThresholdFilter::Pointer labelThreshold = BinaryThresholdFilter::New();
  labelThreshold->SetInput( relabel->GetOutput() );
  labelThreshold->SetInsideValue(1);
  labelThreshold->SetOutsideValue(0);
  labelThreshold->SetLowerThreshold(1); // Only the largest label
  labelThreshold->SetUpperThreshold(1); // Only the largest label
  labelThreshold->Update();

  ByteImageType::Pointer largestLabel = labelThreshold->GetOutput();

  ByteImageType::Pointer dilateImage = largestLabel;
  if( closingSize > 0 )
    {
    // Dilate mask
    StructElementType closeStruct;
    closeStruct.SetRadius(openingSize + closingSize);
    closeStruct.CreateStructuringElement();
    DilateType::Pointer dil = DilateType::New();
    dil->SetDilateValue(1);
    dil->SetKernel(closeStruct);
    dil->SetInput( largestLabel );
    dil->Update();
    dilateImage = dil->GetOutput();
    }

  ByteImageType::Pointer safetyImage = dilateImage;
  const int              finalOpeningSize = ( closingSize - safetySize ) > 0 ? ( closingSize - safetySize ) : 0;
  if( finalOpeningSize > 0 )
    {
    // Binary Erode
    StructElementType remainderStruct;
    remainderStruct.SetRadius(finalOpeningSize);
    remainderStruct.CreateStructuringElement();
    ErodeType::Pointer finalErode = ErodeType::New();
    finalErode->SetErodeValue(1);
    finalErode->SetKernel(remainderStruct);
    finalErode->SetInput( dilateImage );
    finalErode->Update();
    safetyImage = finalErode->GetOutput();
    }

  typedef itk::MultiplyImageFilter<ByteImageType, ByteImageType, ByteImageType> multFilter;
  multFilter::Pointer myMult1 = multFilter::New();
  myMult1->SetInput1( safetyImage );
  myMult1->SetInput2( Mask );
  myMult1->Update();

  multFilter::Pointer myMult2 = multFilter::New();
  myMult2->SetInput1( myMult1->GetOutput() );
  myMult2->SetInput2(inputLabelImage);
  myMult2->Update();
  return myMult2->GetOutput();
}

itk::Image<unsigned char, 3>::Pointer
ExtractSingleLargestRegion(const unsigned char threshold_low, const unsigned char threshold_high,
                           const int openingSize, const int closingSize, const int safetySize,
                           itk::Image<unsigned char, 3>::Pointer inputLabelImage)
{
  typedef itk::Image<unsigned char, 3>                                  ByteImageType;
  typedef itk::BinaryThresholdImageFilter<ByteImageType, ByteImageType> BinaryThresholdFilter;
  BinaryThresholdFilter::Pointer threshold = BinaryThresholdFilter::New();
  threshold->SetInput(inputLabelImage);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->SetLowerThreshold(threshold_low);
  threshold->SetUpperThreshold( threshold_high );
  threshold->Update();

  return ExtractSingleLargestRegionFromMask(
    threshold->GetOutput(), openingSize, closingSize, safetySize, inputLabelImage);
}
