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
#include <string>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiLabelSTAPLEImageFilter.h"
#include "itkTimeProbe.h"
#include "STAPLEAnalysisCLP.h"
#include <BRAINSCommonLib.h>


template <unsigned int ImageDimension>
int STAPLE(unsigned int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  itk::TimeProbe timer;
  timer.Start();

  using LabelImageType = itk::Image<unsigned int, ImageDimension>;

  using FilterType = itk::MultiLabelSTAPLEImageFilter<LabelImageType, LabelImageType>;

  typename FilterType::Pointer filter = FilterType::New();

  using StringVectorType = std::vector<std::string>;
  using StringVectorIteratorType = typename StringVectorType::const_iterator;

  StringVectorIteratorType currentLabel = inputLabelVolume.begin();
  for( unsigned int i = 0;
       i < inputLabelVolume.size();
       i++, ++currentLabel )
    {
    using ReaderType = itk::ImageFileReader<LabelImageType>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( *currentLabel );
    reader->Update();

    filter->SetInput( i, reader->GetOutput() );
    }
  filter->Update();
  for( unsigned int i = 0; i < inputLabelVolume.size(); i++ )
    {
//    std::cout << "Confusion matrix assessment: " << argv[i] << std::endl;
//    std::cout << filter->GetConfusionMatrix( i-3 ) << std::endl << std::endl;

    typename FilterType::ConfusionMatrixType conf =
      filter->GetConfusionMatrix( i );
    for( unsigned int j = 0; j < conf.rows(); j++ )
      {
      for( unsigned int k = 0; k < conf.cols(); k++ )
        {
        std::cout << conf( j, k );
        if( j == conf.rows() - 1 && k == conf.cols() - 1 )
          {
          std::cout << std::endl;
          }
        else
          {
          std::cout << ",";
          }
        }
      }
    }

  typename FilterType::PriorProbabilitiesType priors =
    filter->GetPriorProbabilities();

//  std::cout << "Prior probabilities" << std::endl;
//  for ( unsigned int i = 0; i < priors.GetSize(); ++i )
//    {
//    std::cout << "\t" << i << ": " << priors[i] << "\n";
//    }

  using WriterType = itk::ImageFileWriter<LabelImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  timer.Stop();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension outputImage "
              << "segmentationImage1 ... segmentationImageN" << std::endl;
    return EXIT_FAILURE;
    }

  switch( inputDimension )
    {
    case 2:
      STAPLE<2>( argc, argv);
      break;
    case 3:
      STAPLE<3>( argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
