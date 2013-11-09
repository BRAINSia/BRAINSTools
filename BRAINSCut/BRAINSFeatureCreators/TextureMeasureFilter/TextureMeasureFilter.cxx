/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: CannySegmentationLevelSetImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 21:44:37 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

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

#include "itkScalarImageToTextureFeaturesFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "TextureMeasureFilterCLP.h"
#include "itkCastImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include <iostream>
#include <fstream>
#include <BRAINSCommonLib.h>

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if( inputVolume == "na"  || inputMaskVolume == "na" )
    {
    std::cout << "InputVolume Required!"
              << std::endl;
    }

  const     unsigned int Dimension = 3;

  typedef   double                              InputPixelType;
  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef itk::ImageFileReader<InputImageType>  InputImageReaderType;

  InputImageReaderType::Pointer inputImageReader = InputImageReaderType::New();

  inputImageReader->SetFileName( inputVolume );
  try
    {
    inputImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception Object Caught! " << std::endl;
    std::cerr << err << std::endl;
    throw;
    }

  /** Convert Type to unsigned char */
  typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType>
    RescaleFilterType;

  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

  rescaler->SetInput( inputImageReader->GetOutput() );
  rescaler->SetOutputMaximum(255);
  rescaler->SetOutputMinimum(0);

  /** read ROI mask */
  typedef unsigned char                            InternalPixelType;
  typedef itk::Image<InternalPixelType, Dimension> InternalImageType;
  typedef itk::CastImageFilter<InputImageType, InternalImageType>
    CasterType;

  CasterType::Pointer caster = CasterType::New();

  caster->SetInput( inputImageReader->GetOutput() );

  typedef itk::ImageFileReader<InputImageType> BinaryImageReaderType;

  BinaryImageReaderType::Pointer binaryImageReader = BinaryImageReaderType::New();

  binaryImageReader->SetFileName( inputMaskVolume);

  /** make sure binary mask */
  typedef itk::BinaryThresholdImageFilter<InputImageType, InternalImageType>
    BinaryFilterType;
  BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();

  binaryFilter->SetInput( binaryImageReader->GetOutput() );
  binaryFilter->SetLowerThreshold(1.0);
  binaryFilter->SetInsideValue(1);
  binaryFilter->SetOutsideValue(1);

  /** texture Computation */

  typedef itk::Statistics::ScalarImageToTextureFeaturesFilter<InternalImageType>
    TextureFilterType;

  TextureFilterType::FeatureNameVectorPointer requestedFeatures =
    TextureFilterType::FeatureNameVector::New();

  /** texture setting */
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::Energy );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::Entropy );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::Correlation );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::InverseDifferenceMoment );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::ClusterShade );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::ClusterProminence );
  requestedFeatures->push_back( TextureFilterType::TextureFeaturesFilterType::HaralickCorrelation );

  TextureFilterType::Pointer textureFilter = TextureFilterType::New();

  /** texture distance */
  TextureFilterType::OffsetVectorPointer requestedOffsets =
    TextureFilterType::OffsetVector::New();

  TextureFilterType::OffsetVector::ConstIterator offSetIt;
  std::cout << "TextureDistance,"
            << distance
            << std::endl;
  for( offSetIt  = textureFilter->GetOffsets()->Begin();
       offSetIt != textureFilter->GetOffsets()->End();
       ++offSetIt )
    {
    TextureFilterType::OffsetType tempOffset = offSetIt.Value();
    for( unsigned int t = 0; t < Dimension; t++ )
      {
      tempOffset[t] *= distance;
      }
    requestedOffsets->push_back( tempOffset );
    }

  textureFilter->SetInput( caster->GetOutput() );
  textureFilter->SetMaskImage( binaryFilter->GetOutput() );
  // textureFilter->SetInsidePixelValue( 1 );
  textureFilter->SetNumberOfBinsPerAxis(256);
  textureFilter->FastCalculationsOff();
  textureFilter->SetRequestedFeatures( requestedFeatures);
  textureFilter->SetOffsets( requestedOffsets );
  try
    {
    textureFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception Object Caught! " << std::endl;
    std::cerr << err << std::endl;
    throw;
    }

  /** Get output and display */

  TextureFilterType::FeatureValueVectorPointer means, stds;

  means = textureFilter->GetFeatureMeans();
  stds = textureFilter->GetFeatureStandardDeviations();

  const TextureFilterType::FeatureNameVector* featureNames = textureFilter->GetRequestedFeatures();

  TextureFilterType::FeatureValueVector::ConstIterator mIt;
  TextureFilterType::FeatureValueVector::ConstIterator sIt;
  TextureFilterType::FeatureNameVector::ConstIterator  nIt;

  /** output file ready */
  std::ofstream outputFileStream;
  try
    {
    outputFileStream.open( outputFilename.c_str() );
    }
  catch( std::exception& ex )
    {
    std::cout << "Opening file failed: " << outputFilename << std::endl;
    std::cout << ex.what() << std::endl;
    exit(EXIT_FAILURE);
    }

  int counter;
  for( counter = 0, mIt = means->Begin(), sIt = stds->Begin(), nIt = featureNames->Begin();
       mIt != means->End(); ++mIt, counter++, ++nIt, ++sIt )
    {
    int         name = nIt.Value();
    std::string outPrefixString("Feature, ");

    switch( name )
      {
      case 0:
        {
        outPrefixString +=  "Energy                  ,";
        }
        break;
      case 1:
        {
        outPrefixString +=  "Entropy                  ,";
        }
        break;
      case 2:
        {
        outPrefixString +=  "Correlation              ,";
        }
        break;
      case 3:
        {
        outPrefixString +=  "InverseDifferenceMoment  ,";
        }
        break;
      case 4:
        {
        outPrefixString +=  "Inertia                  ,";
        }
        break;
      case 5:
        {
        outPrefixString +=  "ClusterShade             ,";
        }
        break;
      case 6:
        {
        outPrefixString +=  "ClusterProminence        ,";
        }
        break;
      case 7:
        {
        outPrefixString +=  "HaralickCorrelation      ,";
        }
        break;
      default:
        outPrefixString +=  "Unknown Feature          ,";
      }
    outPrefixString +=  "Mean,";
    outputFileStream << outPrefixString;
    std::cout << outPrefixString;

    outputFileStream << mIt.Value();
    std::cout << mIt.Value();

    outputFileStream << ", std, ";
    std::cout << ", std, ";

    outputFileStream << sIt.Value();
    std::cout << sIt.Value();

    outputFileStream << "\n";
    std::cout << "\n";
    }

  outputFileStream.close();
  return 0;
}
