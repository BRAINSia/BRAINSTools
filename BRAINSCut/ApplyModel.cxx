/** \file
 * \ingroup Main
 */

/**
 * \defgroup AM Apply Model
 * \ingroup Main
 */

#include "Utilities.h"
#include "itksys/SystemTools.hxx"
#include <NetConfiguration.h>
#include "NeuralParams.h"
#include "ANNParams.h"
#include "ApplyModel.h"
#include "NetConfigurationParser.h"
#include "itkIO.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_sample.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <itkHistogramMatchingImageFilter.h> // KEY_HE
/*Regina
  CleanUpOverlapArea is check overlapping area btw structures,
  usually btw adjacent structures, and
  choose the structure which has maximum ANN signal and
  surpress others below the given threshold value
 */
void         CleanUpOverlapArea(NetConfiguration & ANNConfiguration,
                                DataSet *curDataSet,
                                const std::string & fOutputDir,
                                const neural_scalar_type MaskThresh,
                                const int OutputVectorSize,
                                int verbose)
{
  if( verbose )
    {
    std::cerr << "This will print out process in more detail for debugging purpose\n"
              << std::endl;
    }
  // 1. Read every structure's ann cutout and store.
  std::cout << "Cleaning Mask Overlapping Area\n";
  ProbabilityMapList *probabilityMaps =    ANNConfiguration.Get<ProbabilityMapList>(
      "ProbabilityMapList");
  int structure = -1;

  std::string ANNCutOut_name[OutputVectorSize];
  std::string StructureID[OutputVectorSize];

  InternalImageType::Pointer ANNCutOut_Image[OutputVectorSize];
  const std::string          OutputDir = fOutputDir + "/";
  std::string                Type;
  std::string                ExtI;
  const std::string          ImageID( curDataSet->GetAttribute<StringValue>("Name") );

  typedef itk::AddImageFilter<InternalImageType, InternalImageType, InternalImageType> SumImagesType;
  SumImagesType::Pointer Summation_ANNCutOut = SumImagesType::New();
  for( ProbabilityMapList::iterator pmi = probabilityMaps->begin();
       pmi != probabilityMaps->end();
       ++pmi )
    {
    structure++;

    ProbabilityMapParser *probabilityMap    = dynamic_cast<ProbabilityMapParser *>( pmi->second );

    StructureID[structure] = probabilityMap->GetAttribute<StringValue>("StructureID");
    Type = "ANNContinuousPredictionMedianFiltered_";
    ExtI = ".nii.gz";

    ANNCutOut_name[structure] = OutputDir +  Type + StructureID[structure] + ImageID + ExtI;
    std::cout << " Read Image :: " << ANNCutOut_name[structure] << std::endl;
    try
      {
      ANNCutOut_Image[structure]     =
        itkUtil::ReadImage<InternalImageType>(ANNCutOut_name[structure]);
      }
    catch( ... )
      {
      std::cout << " Can't open file :: " << ANNCutOut_name[structure]
                << std::endl;
      }
    if( structure == 1 )    // Add First and Second Images
      {
      try
        {
        Summation_ANNCutOut->SetInput1(ANNCutOut_Image[0]);
        Summation_ANNCutOut->SetInput2(ANNCutOut_Image[1]);

        Summation_ANNCutOut->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Summation of ANNCutout:: exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        }
      }
    else if( structure > 1 )
      {
      std::cout << "Add " << structure + 1 << "th ANNCutout ... ... \n";
      Summation_ANNCutOut->SetInput1( Summation_ANNCutOut->GetOutput() );
      Summation_ANNCutOut->SetInput2(ANNCutOut_Image[structure]);
      try
        {
        Summation_ANNCutOut->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Summation of ANNCutout:: exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        }
      }
    }
  // 2. Find pixel which has more than one values above threshold.
  // iterate entire image and compare ANN output signal across the structures
  std::cout << " Iterator ... \n";
  typedef itk::ImageRegionIterator<InternalImageType> IteratorType;

  IteratorType Iterator_SumANNCutout;

  Iterator_SumANNCutout =
    IteratorType( Summation_ANNCutOut->GetOutput(),
                  Summation_ANNCutOut->GetOutput()->GetRequestedRegion() );
  Iterator_SumANNCutout.GoToBegin();
  /*for ( int strt = 0; strt < OutputVectorSize; strt++ )
  {
    strt_iter[strt] = IteratorType( ANNCutOut_Image[strt],
                                    ANNCutOut_Image[strt]->GetRequestedRegion() );
    strt_iter[strt].GoToBegin();
  }*/
  std::cout << " Find Maximum ... \n";
  while( !Iterator_SumANNCutout.IsAtEnd() )
    {
    if( Iterator_SumANNCutout.Get() > MaskThresh )  // Considering every voxel
      {                                             // summed up over threhsold
                                                    // value
                                                    // Fine maximum
      neural_scalar_type max = 0;
      int                max_structure = -1;
      // scan each ANNCutout....
      for( int strt = 0; strt < OutputVectorSize; strt++ )
        {
        if( ANNCutOut_Image[strt]->GetPixel( Iterator_SumANNCutout.GetIndex() ) > max )
          {
          max = ANNCutOut_Image[strt]->GetPixel( Iterator_SumANNCutout.GetIndex() );
          max_structure = strt;
          }
        }
      if( max < MaskThresh )
        {
        ANNCutOut_Image[max_structure]->SetPixel(Iterator_SumANNCutout.GetIndex(), MaskThresh + 0.01F);
        }
      for( int strt = 0; strt < OutputVectorSize; strt++ )
        {
        if( ANNCutOut_Image[strt]->GetPixel( Iterator_SumANNCutout.GetIndex() ) > MaskThresh && strt != max_structure )
          {
          ANNCutOut_Image[strt]->SetPixel(Iterator_SumANNCutout.GetIndex(), MaskThresh * 0.9F);
          }
        }
      }
    ++Iterator_SumANNCutout;
    }

  // Write ANN output and masks.
  std::cout << "Write ANN output and Mask after Cleaning... ...\n";
  for( int strt = 0; strt < OutputVectorSize; strt++ )
    {
    /*
     * Write Debugging Image
     */
    if( verbose > 5 )
      {
      Type = "ANNCutOut_median_cleanup";
      ExtI = ".nii.gz";

      std::string WriteID = OutputDir +  Type + StructureID[strt] + ImageID + ExtI;
      itkUtil::WriteImage<InternalImageType>(ANNCutOut_Image[strt], WriteID);
      }

    typedef itk::Image<unsigned char, 3> BinaryMaskImageType;
    BinaryMaskImageType::Pointer outputMask;
      {
      itk::ThresholdImageFilter<InternalImageType>::Pointer filterUpper =
        itk::ThresholdImageFilter<InternalImageType>::New();
      filterUpper->SetInput(ANNCutOut_Image[strt]);
      filterUpper->SetOutsideValue(255);    // Foreground
      // filterUpper-> SetInsideValue(0); //Background
      filterUpper->ThresholdAbove(MaskThresh);
      filterUpper->Update();
      outputMask = itkUtil::TypeCast<InternalImageType, BinaryMaskImageType>(
          filterUpper->GetOutput() );
      }
      {
      // Find the largest filled region with high probability
      // NOTE:  The most robust way to do this would be to find the largest
      // background labeled image,
      //       and then choose one of those locations as the seed.
      // For now just choose all the corners as seed points
      BinaryMaskImageType::SizeType ImageSize =
        outputMask->GetLargestPossibleRegion().GetSize();
      typedef itk::ConnectedThresholdImageFilter<BinaryMaskImageType,
                                                 BinaryMaskImageType> ConnectedThresholdFilterType;
      ConnectedThresholdFilterType::Pointer ConnectedThresholdFilter =
        ConnectedThresholdFilterType::New();
        {
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { 0, 0, 0}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { ImageSize[0] - 1, 0, 0}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { 0, ImageSize[1] - 1, 0}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { ImageSize[0] - 1, ImageSize[1] - 1, 0}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { 0, 0, ImageSize[2] - 1}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { ImageSize[0] - 1, 0, ImageSize[2] - 1}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation = { { 0, ImageSize[1] - 1, ImageSize[2] - 1}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
          {
          const BinaryMaskImageType::IndexType SeedLocation =
                { { ImageSize[0] - 1, ImageSize[1] - 1, ImageSize[2] - 1}};
          ConnectedThresholdFilter->SetSeed(SeedLocation);
          }
        }
      ConnectedThresholdFilter->SetReplaceValue(1);
      ConnectedThresholdFilter->SetUpper(0);
      ConnectedThresholdFilter->SetLower(0);
      ConnectedThresholdFilter->SetInput(outputMask);
      ConnectedThresholdFilter->Update();
      BinaryMaskImageType::Pointer ConnectedThresholdImage =
        ConnectedThresholdFilter->GetOutput();
      /*
       * write a debugging image
       */
      if( verbose > 5 )
        {
        const std::string ConnectedThresholdWriteID = OutputDir
          + Type
          + StructureID[strt]
          + ImageID
          +
          "_ConnectedThreshold.nii.gz";
        itkUtil::WriteImage<BinaryMaskImageType>(ConnectedThresholdImage,
                                                 ConnectedThresholdWriteID);
        }

      itk::ConnectedComponentImageFilter<BinaryMaskImageType, BinaryMaskImageType>
      ::Pointer myComponentsFilter =
        itk::ConnectedComponentImageFilter<BinaryMaskImageType,
                                           BinaryMaskImageType>::New();
      myComponentsFilter->SetInput(outputMask);
      myComponentsFilter->Update();

      BinaryMaskImageType::Pointer labeledImage = myComponentsFilter->GetOutput();
      /*
       * write a debugging image
       */
      if( verbose > 5 )
        {
        const std::string LabelWriteID = OutputDir + Type
          + StructureID[strt] + ImageID
          + "_labeled.nii.gz";
        itkUtil::WriteImage<BinaryMaskImageType>(labeledImage, LabelWriteID);
        }
      typedef itk::RelabelComponentImageFilter<BinaryMaskImageType,
                                               BinaryMaskImageType> RelabelType;
      RelabelType::Pointer relabel = RelabelType::New();
      relabel->SetInput(labeledImage);
      // RELAbel->SetMinimumObjectSize( HundredPercentRegionSize ); //Use the
      // size of the 100% region as a minimum.
      try
        {
        relabel->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Relabel: exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        }
      BinaryMaskImageType::Pointer MultiLabeledImage = relabel->GetOutput();

      /*
       * write a debugging image
       */
      if( verbose > 5 )
        {
        const std::string MultiLabelWriteID = OutputDir + Type
          + StructureID[strt] + ImageID
          + "_MultiLabeled.nii.gz";
        itkUtil::WriteImage<BinaryMaskImageType>(MultiLabeledImage,
                                                 MultiLabelWriteID);
        }
      // Thresholding Image so that the output Image have only one labeling map.
      typedef itk::BinaryThresholdImageFilter<BinaryMaskImageType,
                                              BinaryMaskImageType> ThresholdType;
      ThresholdType::Pointer thresholdLabel =  ThresholdType::New();
      thresholdLabel->SetInput( relabel->GetOutput() );
      thresholdLabel->SetLowerThreshold(1);
      thresholdLabel->SetUpperThreshold(1);
      thresholdLabel->Update();

      BinaryMaskImageType::Pointer OneLabeledMask =  thresholdLabel->GetOutput();

      // Openeing operation

      typedef itk::BinaryBallStructuringElement<BinaryMaskImageType::PixelType,
                                                3> KernelType;
      KernelType           ball;
      KernelType::SizeType ballSize;

      ballSize.Fill(20);
      ball.SetRadius(ballSize);
      ball.CreateStructuringElement();
      typedef itk::BinaryMorphologicalOpeningImageFilter<BinaryMaskImageType,
                                                         BinaryMaskImageType,
                                                         KernelType> ClosingFilterType;
      ClosingFilterType::Pointer closing = ClosingFilterType::New();
      closing->SetInput( OneLabeledMask );
      // DOING

      closing->SetKernel( ball );
      closing->SetForegroundValue(1);
      closing->Update();

      // TODO:: Give Proper Name if it is given in the xml file (net
      // configuration file)
      std::string OneLabelWriteID_Mask(
        curDataSet->GetMaskFilenameByType( StructureID[strt] ) );

      if( OneLabelWriteID_Mask == "" or OneLabelWriteID_Mask == "na" )
        {
        std::cout << " MaskName is " << OneLabelWriteID_Mask << std::endl;
        OneLabelWriteID_Mask = OutputDir
          + Type
          + StructureID[strt]
          + ImageID
          + "_OneLabeled.nii.gz";
        }
      std::cout << "Resulting ANN masks for test image " << ImageID.c_str()
                << " write to disk as " << OneLabelWriteID_Mask << std::endl;
      itkUtil::WriteImage<BinaryMaskImageType>(closing->GetOutput(),
                                               OneLabelWriteID_Mask);
      }
    }
}

// /Apply Model
/** detailed description of cutout
 * \ingroup AM
 */
int ApplyModel(NetConfiguration & ANNConfiguration,
               int verbose,
               bool /*TODO: DELETE  histogramEqualization */,
               bool /*TODO: DELETE multiStructureThreshold*/
               )
{
  /** generate registration */

  RegistrationConfigurationParser *regParams =
    ANNConfiguration.Get<RegistrationConfigurationParser>("RegistrationConfigurationParser");
  const std::string regID( regParams->GetAttribute<StringValue>("ID") );
  const std::string imageTypeToUse( regParams->GetAttribute<StringValue>(
                                      "ImageTypeToUse") );
  std::string ANNModelFilename;

  // Get the atlas dataset
  DataSet *         atlasDataSet = ANNConfiguration.GetAtlasDataSet();
  const std::string atlasImage( atlasDataSet->GetImageFilenameByType(
                                  imageTypeToUse.c_str() ) );

  NeuralParams *model = ANNConfiguration.Get<NeuralParams>("NeuralNetParams");

  if( model == 0 )
    {
    itkGenericExceptionMacro(<< "Model information (gaussian size, iris size, "
                             << "mask smoothing value, training vector filename) not given,"
                             << " so training vectors cannot be created.");
    }
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "** Using Model file: "
            << model->GetAttribute<StringValue>("TrainingModelFilename") <<  std::endl;
  std::cout << __FILE__ << " " << __LINE__ << std::endl;

  ApplyModelType *         applyModel = ANNConfiguration.Get<ApplyModelType>("ApplyModel");
  const neural_scalar_type MaskThresh = applyModel->GetAttribute<FloatValue>("MaskThresh");

  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << "Mask Thrshold Value :: " << MaskThresh << "\n";
  if( ( MaskThresh < MINIMUM_VALUE + PERCENT_MIN_MAX_TOLERANCE )
      || ( MaskThresh > ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) ) )
    {
    itkGenericExceptionMacro(<< "ERROR: Invalid mask threshold given as " << MaskThresh
                             << std::endl
                             << "Valid values are in range : ["
                             << MINIMUM_VALUE + PERCENT_MIN_MAX_TOLERANCE
                             << ","
                             << HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE
                             << "]");
    }

  /** Get Input Vector Size */
  ProbabilityMapList * ROIList = ANNConfiguration.Get<ProbabilityMapList>("ProbabilityMapList");

  int outputVectorSize = ROIList->size();

  unsigned int inputVectorSize = InputVectorSizeRequirement( outputVectorSize,
                                                             (ANNConfiguration.GetAtlasDataSet()->GetImageTypes() ).
                                                             size(),
                                                             model->GetAttribute<IntValue>("GradientProfileSize") );

  /** number of thread = 1 hard coded */
  const unsigned int numThreads = 1;

  /** generate registration from atlas to subject direction */
  GenerateRegistrations(ANNConfiguration, true, true,  numThreads);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Process
  //  -- Fixe the order of ROIs
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  std::map<int, std::string> MapOfROIOrder;
  int                        roiCounter = 0;
  for( ProbabilityMapList::iterator pmi = ROIList->begin();
       pmi != ROIList->end();
       ++pmi )
    {
    MapOfROIOrder.insert( std::pair<int, std::string>( roiCounter, pmi->first) );
    roiCounter++;
    }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // interate all the subjects
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  NetConfiguration::ApplyDataSetListType inputSubjectsList;

  inputSubjectsList = NetConfiguration::ApplyDataSetListType(
      ANNConfiguration.GetApplyDataSets() );
  for( NetConfiguration::ApplyDataSetListType::iterator subjectIt = inputSubjectsList.begin();
       subjectIt != inputSubjectsList.end();
       ++subjectIt )
    {
    AddSubjectInputVector(
      (*subjectIt),
      ANNConfiguration,
      regID,
      inputVectorSize,
      outputVectorSize,
      MapOfROIOrder,
      true);                      // last parameter: apply=false

    CleanUpOverlapArea(ANNConfiguration,
                       ( *subjectIt),
                       ( *subjectIt)->GetAttribute<StringValue>("OutputDir"),
                       MaskThresh,
                       outputVectorSize,
                       verbose);
    }

  // After Every Structure has been traced, Deal with Overlapping area.
  // Read ANNCutOut_median(ANN output) for each structure for one subject

  // Scan each voxel to find if there is more then one value greater than
  // threshold
  // surpress second top one

  return 0;
}

int ApplyModel(const std::string & XMLFile,
               int verbose,
               bool histogramEqualization,
               bool multiStructureThreshold
               )
{
  NetConfiguration *     ANNConfiguration;
  NetConfigurationParser ANNConfigurationParser = NetConfigurationParser( XMLFile);

  try
    {
    ANNConfiguration = ANNConfigurationParser.GetNetConfiguration();
    }
  catch( BRAINSCutExceptionStringHandler & ex )
    {
    throw;
    }
  if( ANNConfiguration->Verify() != true )
    {
    std::cerr << "XML file " << " is invalid." << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    ANNConfiguration->PrintSelf(std::cerr, 0);
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    return -1;
    }
  int rval;
  try
    {
    rval = ApplyModel(*ANNConfiguration,
                      verbose,
                      histogramEqualization,
                      multiStructureThreshold
                      );
    }
  catch( BRAINSCutExceptionStringHandler & ex )
    {
    std::cerr << ex.Error() << std::endl;
    rval = -1;
    }
  catch( ... )
    {
    std::cerr << "Unidentified exception" << std::endl;
    rval = -1;
    }
  return rval;
}
