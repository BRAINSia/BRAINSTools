#include "simpleEMCLP.h"

#include <itksys/SystemTools.hxx>

#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNumericTraits.h"
#include "itkMultiplyImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVersion.h"
#include "AtlasCropImageSource.h"

#include <vector>

// Use manually instantiated classes for the big program chunks
// #define MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"

// #undef MU_MANUAL_INSTANTIATION

#include "filterFloatImages.h"

#include <iostream>
#include <string>
#include <sstream>

#include <cstdlib>
// template <class inputPixelType,class outputPixelType,class priorPixelType>
template <class priorPixelType, class inputPixelType>
int simpleRunEMS( std::string t1Volume,
                  std::string t2Volume,
                  std::string pdVolume,
                  std::string templateVolume,
                  int maxIteration,
                  std::string OutputFileNamePrefix,
                  std::vector<std::string> &  priorsList,
                  std::vector<double> &  priorsWeightList,
                  bool warp,
                  int degreeOfBiasFieldCorrection,
                  double likelihoodTolerance)
{
  const int status = -1;
  // PARSE_ARGS;
  /** read parameters            **/
  /** - read target input images **/
  std::string T1FileName(t1Volume);
  std::string T2FileName(t2Volume);
  std::string PDFileName(pdVolume);

  typedef itk::Image<float, 3>         FloatImageType;
  typedef itk::Image<unsigned char, 3> ByteImageType;
  typedef itk::Image<short, 3>         ShortImageType;

  typedef FloatImageType::Pointer FloatImagePointer;
  typedef ByteImageType::Pointer  ByteImagePointer;
  typedef ShortImageType::Pointer ShortImagePointer;

  typedef itk::Image<inputPixelType, 3>                                    InputImageType;
  typedef itk::ImageFileReader<InputImageType>                             InputImageReaderType;
  typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType> RescaleImageFilterType;
  typedef typename InputImageType::Pointer                                 InputImagePointer;

  typedef itk::ImageFileReader<InputImageType> templateReaderType;
  typename templateReaderType::Pointer templateReader = templateReaderType::New();
  templateReader->SetFileName(templateVolume);
  templateReader->Update();
  std::vector<InputImagePointer> images;

  /* IMAGES MUST BE Non-Negative or the results will fail. */
  std::vector<std::string> inputImageFilenames;
  if( T1FileName == "" )
    {
    std::cerr << " T1 file is missing, other modality will be used as a reference image.\n  ";
    }
  else
    {
    typename InputImageReaderType::Pointer inputReaderT1 = InputImageReaderType::New();
    inputReaderT1->SetFileName(T1FileName);
    inputImageFilenames.push_back(T1FileName);
    inputReaderT1->Update();
    typename RescaleImageFilterType::Pointer rescalerT1 = RescaleImageFilterType::New();
    // Set the upper limit to 4096 in the case of floating point data.
    const inputPixelType outMin = 0;
    const inputPixelType outMax
      = ( vcl_numeric_limits<inputPixelType>::max() > 4096 ) ? 4096 : vcl_numeric_limits<inputPixelType>::max();
    rescalerT1->SetOutputMinimum(outMin);
    rescalerT1->SetOutputMaximum(outMax);
    rescalerT1->SetInput( inputReaderT1->GetOutput() );
    rescalerT1->Update();
    images.push_back( rescalerT1->GetOutput() );
    }

  if( T2FileName == "" )
    {
    std::cerr << " T2 file is missing.\n  ";
    }
  else
    {
    typename InputImageReaderType::Pointer inputReaderT2 = InputImageReaderType::New();
    inputReaderT2->SetFileName(T2FileName);
    inputImageFilenames.push_back(T2FileName);
    inputReaderT2->Update();
    typename RescaleImageFilterType::Pointer rescalerT2 = RescaleImageFilterType::New();
    // Set the upper limit to 4096 in the case of floating point data.
    const inputPixelType outMin = 0;
    const inputPixelType outMax
      = ( vcl_numeric_limits<inputPixelType>::max() > 4096 ) ? 4096 : vcl_numeric_limits<inputPixelType>::max();
    rescalerT2->SetOutputMinimum(outMin);
    rescalerT2->SetOutputMaximum(outMax);
    rescalerT2->SetInput( inputReaderT2->GetOutput() );
    rescalerT2->Update();
    images.push_back( rescalerT2->GetOutput() );
    }

  if( PDFileName == "" )
    {
    std::cerr << " PD file is missing.\n  ";
    }
  else
    {
    typename InputImageReaderType::Pointer inputReaderPD = InputImageReaderType::New();
    inputReaderPD->SetFileName(PDFileName);
    inputImageFilenames.push_back(PDFileName);
    inputReaderPD->Update();
    typename RescaleImageFilterType::Pointer rescalerPD = RescaleImageFilterType::New();
    // Set the upper limit to 4096 in the case of floating point data.
    const inputPixelType outMin = 0;
    const inputPixelType outMax
      = ( vcl_numeric_limits<inputPixelType>::max() > 4096 ) ? 4096 : vcl_numeric_limits<inputPixelType>::max();
    rescalerPD->SetOutputMinimum(outMin);
    rescalerPD->SetOutputMaximum(outMax);
    rescalerPD->SetInput( inputReaderPD->GetOutput() );
    rescalerPD->Update();
    images.push_back( rescalerPD->GetOutput() );
    }
  if( T1FileName == "" && T2FileName == "" && PDFileName == "" )
    {
    std::cerr << " We need at lease one or more target image.\n";
    return status;
    }
  // Compute the brain outline in order to cut images later after bias
  // correction
  const unsigned int closingSize = 7;
  const float        otsuPercentileThreshold = 0.01;
  //  typename InputImageType::Pointer HeadOutlineMaskImage =
  // FindLargestForgroundFilledMask<InputImageType>( images[0],
  // otsuPercentileThreshold, closingSize );
  typedef itk::LargestForegroundFilledMaskImageFilter<InputImageType> LFFMaskFilterType;
  typename LFFMaskFilterType::Pointer LFF = LFFMaskFilterType::New();
  LFF->SetInput(images[0]);
  LFF->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  LFF->SetClosingSize(closingSize);
  LFF->Update();
  typename InputImageType::Pointer HeadOutlineMaskImage = LFF->GetOutput();

  std::cerr << " Following image will be used as a reference target image \n";
  std::cerr << " Input Image 1:: " << images[0]->GetLargestPossibleRegion().GetSize() << std::endl;

  /** - read atlas images        **/
  /* atlas should have list of atlas for interest and background */
  /* and those atlas should including prior weight. */
  std::vector<std::string> PriorsList(priorsList);
  std::vector<double>      PriorWeights(priorsWeightList);

  typedef itk::Image<priorPixelType, 3>        PriorImageType;
  typedef itk::ImageFileReader<PriorImageType> PriorImageReaderType;
  typedef typename PriorImageType::Pointer     PriorImagePointer;

  std::vector<PriorImagePointer> priors;

  if( PriorsList.size() != PriorWeights.size() )
    {
    std::cerr << " Number of List of Prior Image should be same to the Prior Weights List\n";
    return status;
    }

  std::vector<std::string> priorImageFileNames;
  if( PriorsList.empty() )
    {
    std::cerr << " There is no prior image at all, we need at least one image\n";
    return status;
    }
  else
    {
    for( unsigned int i = 0; i < PriorsList.size(); i++ )
      {
      typename PriorImageReaderType::Pointer priorReader = PriorImageReaderType::New();
      std::string name = PriorsList.at(i);
      priorReader->SetFileName(  name  );
      priorImageFileNames.push_back(name);
      priorReader->Update();
      priors.push_back( priorReader->GetOutput() );
      std::cout << " Prior List :: " << i << "  :: " << name << "\n";
      }
    }

  typedef EMSegmentationFilter<InputImageType, PriorImageType> SegFilterType;
  typename SegFilterType::VectorType priorweights( PriorWeights.size() );
  for( unsigned int i = 0; i < PriorWeights.size(); i++ )
    {
    priorweights[i] = PriorWeights.at(i);
    }

  std::cerr <<  "Start segmentation...\n";

  typename SegFilterType::Pointer segfilter = SegFilterType::New();
  segfilter->DebugOn();

  // HACK:  Need for loop around this
  std::vector<InputImagePointer> templateImageList;
  templateImageList.clear();

  templateImageList.push_back( templateReader->GetOutput() );

  segfilter->SetTemplateImages( templateImageList );
  segfilter->SetInputImages(  images  );
  segfilter->SetPriors(  priors  );
  segfilter->SetMaximumIterations( maxIteration );
  segfilter->SetPriorWeights(priorweights);
  segfilter->SetMaxBiasDegree(degreeOfBiasFieldCorrection);
  if( likelihoodTolerance > 0.0 )   // override the default only if provided.
    {
    segfilter->SetLikelihoodTolerance(likelihoodTolerance);
    }
  if( warp )
    {
    segfilter->DoWarpingOn();
    }
  else
    {
    segfilter->DoWarpingOff();
    }
  segfilter->Update();

  // Write the labels
  std::cerr << "Writing labels...\n";
    {
    typedef itk::ImageFileWriter<ByteImageType> OutputWriterType;
    OutputWriterType::Pointer writer = OutputWriterType::New();

    writer->SetInput( segfilter->GetOutput() );
    writer->UseCompressionOn();
    std::string fn = std::string(OutputFileNamePrefix + "_labels.nii.gz");
    writer->SetFileName( fn.c_str() );
    writer->Update();
    }

  // Write the secondary outputs
  if( true )
    {
    std::cerr << "Writing filtered and bias corrected images...\n";
    std::vector<InputImagePointer> imgset = segfilter->GetCorrected();
    for( unsigned i = 0; i < imgset.size(); i++ )
      {
      typedef itk::CastImageFilter<InputImageType, ShortImageType> CasterType;
      typename CasterType::Pointer caster = CasterType::New();
      typename itk::MultiplyImageFilter<InputImageType,
                                        InputImageType>::Pointer multiplyFilter
        = itk::MultiplyImageFilter<InputImageType, InputImageType>::New();
      multiplyFilter->SetInput1(imgset[i]);
      multiplyFilter->SetInput2(HeadOutlineMaskImage);

      caster->SetInput( multiplyFilter->GetOutput() );
      caster->Update();
      std::string fn = std::string(inputImageFilenames[i] + "_corrected.nii.gz");

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      writer->SetInput( caster->GetOutput() );
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();
      writer->Update();
      }

    // Short posteriors
    std::cerr << "Writing posterior images...\n";
    std::vector<ShortImagePointer> probset = segfilter->GetShortPosteriors();
    for( unsigned int i = 0; i < ( probset.size() - 3 ); i++ )
      {
      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();
      writer->SetInput(    probset[i]    );
      writer->SetFileName(priorImageFileNames[i] + "_posterior.nii.gz");
      writer->UseCompressionOn();
      writer->Update();
      }
    }
  return status;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const int status
    = simpleRunEMS<float, float>(t1Volume,
                                 t2Volume,
                                 pdVolume,
                                 templateVolume,
                                 maxIteration,
                                 OutputFileNamePrefix,
                                 priorsList,
                                 priorsWeightList,
                                 warp,
                                 degreeOfBiasFieldCorrection,
                                 likelihoodTolerance );
  /*
      if(inputPixelType=="double") status =
      simpleRunEMS<double>(        t1Volume, t2Volume, pdVolume, templateVolume, maxIteration ,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="short")  status =
      simpleRunEMS<short>(         t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="ushort") status =
      simpleRunEMS<unsigned short>(t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="int")    status =
      simpleRunEMS<int>(           t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="uint")   status =
      simpleRunEMS<unsigned int>(  t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="long")   status =
      simpleRunEMS<long>(          t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="float")  status =
      simpleRunEMS<float>(         t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else if( inputPixelType=="uchar")  status =
      simpleRunEMS<unsigned char>( t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      else( inputPixelType=="char")   status =
      simpleRunEMS<char>(          t1Volume, t2Volume, pdVolume, templateVolume,
      OutputFileNamePrefix, priorsList, priorsWeightList, warp , degreeOfBiasFieldCorrection );
      */
  return status;
}
