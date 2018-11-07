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

#include <iostream>
#include <fstream>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkMetaDataObject.h>
#include <itkIOCommon.h>
#include <itkPointSet.h>
#include <vtkPoints.h>
#include <itkDiffusionTensor3D.h>
#include <itkImageToImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkVector.h>
#include <vnl/vnl_vector.h>
#include <itkTextOutput.h>
#include <itkCommand.h>
#include <itkMath.h>

#include "itkDtiFastMarchingCostFilter.h"
#include "gtractCostFastMarchingCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>
namespace
{
// The following class is used to support callbacks
// on the filter in the pipeline that follows later
class ShowProgressObject
{
public:
  ShowProgressObject(itk::ProcessObject *o)
  {
    m_Process = o;
  }

  void ShowProgress()
  {
    std::cout << "Progress " << m_Process->GetProgress() << std::endl;
  }

  itk::ProcessObject::Pointer m_Process;
};
} // end namespace

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tensor Image: " <<  inputTensorVolume << std::endl;
    std::cout << "Anisotropy Image: " <<  inputAnisotropyVolume << std::endl;
    std::cout << "Starting Seeds LabelMap Image: " <<  inputStartingSeedsLabelMapVolume << std::endl;
    std::cout << "Starting Label: " <<  startingSeedsLabel << std::endl;
    std::cout << "Output Cost Image: " <<  outputCostVolume << std::endl;
    std::cout << "Output Speed Image: " <<  outputSpeedVolume << std::endl;
    std::cout << "Seed Threshold: " <<  seedThreshold << std::endl;
    std::cout << "Anisotropy Weight: " <<  anisotropyWeight << std::endl;
    std::cout << "Stopping Value: " <<  stoppingValue << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  typedef float                                                          PixelType;
  typedef itk::Image<PixelType, 3>                                       AnisotropyImageType;
  typedef itk::Image<PixelType, 3>                                       CostImageType;
  typedef itk::Image<PixelType, 3>                                       SpeedImageType;
  typedef double                                                         TensorElementType;
  typedef itk::DiffusionTensor3D<TensorElementType>                      TensorPixelType;
  typedef itk::Image<TensorPixelType, 3>                                 TensorImageType;
  typedef TensorImageType::IndexType                                     seedIndexType;
  typedef itk::DtiFastMarchingCostFilter<CostImageType, TensorImageType> FloatFMType;

  // Read Tensor Image to set principal eigenvector image
  typedef itk::ImageFileReader<TensorImageType> TensorImageReaderType;
  TensorImageReaderType::Pointer tensorImageReader = TensorImageReaderType::New();
  tensorImageReader->SetFileName( inputTensorVolume );

  try
    {
    tensorImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  TensorImageType::Pointer tensorImage = tensorImageReader->GetOutput();

  // Read Anisotropy Image
  typedef itk::ImageFileReader<AnisotropyImageType> AnisotropyImageReaderType;
  AnisotropyImageReaderType::Pointer anisotropyImageReader = AnisotropyImageReaderType::New();
  anisotropyImageReader->SetFileName(  inputAnisotropyVolume );

  try
    {
    anisotropyImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  AnisotropyImageType::Pointer anisotropyImage = anisotropyImageReader->GetOutput();

  typedef signed short                        MaskPixelType;
  typedef itk::Image<MaskPixelType, 3>        MaskImageType;
  typedef itk::ImageFileReader<MaskImageType> MaskImageReaderType;
  MaskImageReaderType::Pointer startingSeedImageReader = MaskImageReaderType::New();
  startingSeedImageReader->SetFileName( inputStartingSeedsLabelMapVolume );

  try
    {
    startingSeedImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  /* Threshold Starting Label Map */
  typedef itk::ThresholdImageFilter<MaskImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer startingThresholdFilter = ThresholdFilterType::New();
  startingThresholdFilter->SetInput( startingSeedImageReader->GetOutput() );
  startingThresholdFilter->SetLower( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->SetUpper( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->Update();

  MaskImageType::Pointer startingSeedMask = startingThresholdFilter->GetOutput();

  /* Now Generate the Seed points for the Fast Marching Algorithm */
  typedef itk::ImageRegionConstIterator<MaskImageType> ConstMaskIteratorType;
  ConstMaskIteratorType maskIt( startingSeedMask, startingSeedMask->GetLargestPossibleRegion() );
  typedef itk::ImageRegionConstIterator<AnisotropyImageType> ConstAnisotropyIteratorType;
  ConstAnisotropyIteratorType anisoIt( anisotropyImage, anisotropyImage->GetLargestPossibleRegion() );

  typedef FloatFMType::NodeType      NodeType; // NodeType is float type
  typedef FloatFMType::NodeContainer NodeContainer;
  typedef CostImageType::IndexType   seedIndexType;
  typedef std::list<seedIndexType>   SeedListType;
  typedef std::list<float>           SeedValueListType;
  SeedValueListType seedValueList;

  SeedListType seedList;

  unsigned int numSeeds = 0;
  anisoIt.GoToBegin();
  for( maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++anisoIt )
    {
    if( ( maskIt.Get() )  &&  ( anisoIt.Get() >= seedThreshold ) )
      {
      seedList.push_front( maskIt.GetIndex() );
      seedValueList.push_front( 0.0 );
      numSeeds++;
      }
    }

  // setup seed points
  NodeContainer::Pointer seedPoints = NodeContainer::New();
  NodeType               node;
  for( unsigned int j = 0; j < numSeeds; j++ )
    {
    seedIndexType seedIndex;
    seedIndex = seedList.back();
    seedList.pop_back();
    node.SetIndex( seedIndex );
    node.SetValue( 0.0 );
    seedPoints->InsertElement(j, node);
    }

  // Run the Fast Marching Algorithm if seeds are defined
  if( numSeeds > 0 )
    {
    /* Set Parameters for the Fast Marching Calculations */

    FloatFMType::Pointer marcher = FloatFMType::New();
    marcher->SetAlivePoints( seedPoints );
    marcher->SetOverrideOutputInformation( true );
#if  1
    // const CostImageType::SizeType size = tensorImage->GetLargestPossibleRegion().GetSize();
    marcher->SetOutputSize( tensorImage->GetLargestPossibleRegion().GetSize() );
    marcher->SetOutputRegion( tensorImage->GetLargestPossibleRegion() );
    marcher->SetOutputSpacing( tensorImage->GetSpacing() );
    marcher->SetOutputOrigin( tensorImage->GetOrigin() );
    marcher->SetOutputDirection( tensorImage->GetDirection() );
#else
    marcher->SetOutputParametersFromImage( tensorImage );
#endif
    marcher->SetInput( tensorImage );
    marcher->SetAnisotropyImage( anisotropyImage );
    marcher->SetAnisotropyWeight( anisotropyWeight );
    marcher->SetStoppingValue( stoppingValue );
    marcher->SetNormalizationFactor( 1 );
    std::cout << "NormalizationFactor: " << marcher->GetNormalizationFactor() << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    ShowProgressObject                                    progressWatch(marcher);
    itk::SimpleMemberCommand<ShowProgressObject>::Pointer command;
    command = itk::SimpleMemberCommand<ShowProgressObject>::New();
    command->SetCallbackFunction(&progressWatch,
                                 &ShowProgressObject::ShowProgress);
    marcher->AddObserver( itk::ProgressEvent(), command);

    itk::OutputWindow::SetInstance( itk::TextOutput::New().GetPointer() );

    // turn on debugging and update
    // marcher->DebugOn();
    marcher->Update();

    /* Save the Cost and Speed Images */
    typedef itk::ImageFileWriter<CostImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->UseCompressionOn();
    writer->SetInput( marcher->GetOutput() );
    writer->SetFileName( outputCostVolume );
    writer->Write();

    typedef itk::ImageFileWriter<SpeedImageType> SpeedWriterType;
    SpeedWriterType::Pointer writer2 = SpeedWriterType::New();
    writer2->UseCompressionOn();
    writer2->SetInput( marcher->GetOutputSpeedImage() );
    writer2->SetFileName( outputSpeedVolume );
    writer2->Write();
    }
  else
    {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << " !!!Warning: No Seeds Selected!!!" << std::endl;
    }
  return EXIT_SUCCESS;
}
