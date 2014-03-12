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

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/09 14:53:40 $
 Version:   $Revision: 1.0 $

 Copyright (c) University of Iowa Department of Radiology. All rights reserved.
 See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
 for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkNaryRelabelImageFilter.h"
#include "LabelMapsCLP.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Generate a labelmap out of label files:" << std::endl;
  for( unsigned int i = 0; i < labelFileList.size(); i++ )
    {
    std::cout << labelFileList[i] << std::endl;
    }
  std::cout << "---------------------------------------------------" << std::endl;

  unsigned char nLabels = labelFileList.size();

  // define image type
  const int dimension = 3;
  typedef unsigned char                    PixelType;
  typedef itk::Image<PixelType, dimension> ImageType;

  // read Inputimage
  typedef itk::ImageFileReader<ImageType> InputReaderType;

  typedef itk::NaryRelabelImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  for( int i = 0; i < nLabels; i++ )
    {
    InputReaderType::Pointer InputReader = InputReaderType::New();

    // input i
    InputReader->SetFileName( labelFileList[i] );
    InputReader->Update();

    filter->SetInput( i, InputReader->GetOutput() );
    }

  filter->SetBackgroundValue( 0.0 );
  filter->SetIgnoreCollision( 1 );   // has to be 1

  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter<ImageType> OutputWriterType;
  OutputWriterType::Pointer Writer = OutputWriterType::New();
  Writer->SetInput( filter->GetOutput() );
  Writer->SetFileName( outputImageFile.c_str() );
  Writer->Update();

  return EXIT_SUCCESS;
}
