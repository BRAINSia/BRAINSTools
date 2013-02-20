/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    BinaryThresholdImageFilter.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

     Author Eunyoung Regina Kim
     This is modified from itk Binary Threshold Image Filter.
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkXorImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageDuplicator.h"

#include "CleanUpOverlapLabelsCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // check inputs
  if( inputBinaryVolumes.size() != outputBinaryVolumes.size() )
    {
    std::cout << "* Error! input and output has to be matched! "
              << std::endl;
    return EXIT_FAILURE;
    }
  typedef  unsigned char BinaryPixelType;

  const unsigned char Dim = 3;

  typedef itk::Image<BinaryPixelType,  Dim> InputBinaryImageType;

  typedef itk::ImageFileReader<InputBinaryImageType> InputBinaryVolumeReaderType;

  std::vector<InputBinaryImageType::Pointer> inputBinaryVolumeVector;
  std::vector<std::string>::iterator         inputBinaryVolumeStringIt;
  for( inputBinaryVolumeStringIt = inputBinaryVolumes.begin();
       inputBinaryVolumeStringIt < inputBinaryVolumes.end();
       ++inputBinaryVolumeStringIt )
    {
    std::cout << "* Read image: " << *inputBinaryVolumeStringIt << std::endl;
    InputBinaryVolumeReaderType::Pointer reader = InputBinaryVolumeReaderType::New();
    reader->SetFileName( *inputBinaryVolumeStringIt );
    reader->Update();

    inputBinaryVolumeVector.push_back( reader->GetOutput() );
    }

  // Add all labels

  typedef itk::ImageDuplicator<InputBinaryImageType> ImageDuplicatorType;

  ImageDuplicatorType::Pointer imageCopier = ImageDuplicatorType::New();
  imageCopier->SetInputImage( inputBinaryVolumeVector.front() );
  imageCopier->Update();

  InputBinaryImageType::Pointer sumVolume = imageCopier->GetOutput();

  typedef itk::AddImageFilter<InputBinaryImageType, InputBinaryImageType,
                              InputBinaryImageType> AddImageFilterType;

  AddImageFilterType::Pointer adder = AddImageFilterType::New();
  for( unsigned int i = 1; // should start from second image
       i < inputBinaryVolumeVector.size();
       ++i )
    {
    adder->SetInput1( sumVolume );
    adder->SetInput2( inputBinaryVolumeVector[i] );
    adder->Update();

    sumVolume = adder->GetOutput();
    }

  // threshold summed image

  typedef itk::BinaryThresholdImageFilter<InputBinaryImageType,
                                          InputBinaryImageType>  ThresholdFilterType;

  ThresholdFilterType::Pointer thresholder = ThresholdFilterType::New();
  thresholder->SetInput( sumVolume );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetLowerThreshold( 1 );

  thresholder->Update();
  sumVolume = thresholder->GetOutput();

  // process one label at a time

  typedef itk::AndImageFilter<InputBinaryImageType, InputBinaryImageType,
                              InputBinaryImageType> AndImageFilterType;

  typedef itk::XorImageFilter<InputBinaryImageType, InputBinaryImageType,
                              InputBinaryImageType> XORImageFilterType;

  typedef itk::ImageFileWriter<InputBinaryImageType> WriterType;
  for( unsigned int i = 0;
       i < inputBinaryVolumeVector.size();
       ++i )
    {
    // write out overlap between current mask and summed mask
    ThresholdFilterType::Pointer maskThresholder = ThresholdFilterType::New();
    maskThresholder->SetInput( inputBinaryVolumeVector[i] );
    maskThresholder->SetOutsideValue( 0 );
    maskThresholder->SetInsideValue( 1 );
    maskThresholder->SetLowerThreshold( 1 );

    AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
    andFilter->SetInput(0, sumVolume );
    andFilter->SetInput(1, maskThresholder->GetOutput() );
    andFilter->Update();

    std::cout << "Write binary volume : " << outputBinaryVolumes[i] << std::endl;

    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( andFilter->GetOutput() );
    writer->SetFileName( outputBinaryVolumes[i] );
    writer->Update();

    // subtract current mask from the summed mask
    XORImageFilterType::Pointer xorFilter = XORImageFilterType::New();
    xorFilter->SetInput1( sumVolume );
    xorFilter->SetInput2( maskThresholder->GetOutput() );

    xorFilter->Update();

    sumVolume = xorFilter->GetOutput();
    }

  return EXIT_SUCCESS;
}
