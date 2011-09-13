/** \file
 * \ingroup Main
 */

/**
 * \defgroup AM Apply Model
 * \ingroup Main
 */

#include "Utilities.h"
#include "itksys/SystemTools.hxx"
#include <ProcessDescription.h>
#include "NeuralParams.h"
#include "ANNParams.h"
#include "SVMParams.h"
#include "ApplyModel.h"
#include "Parser.h"
#include "itkIO.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_sample.h"
#include "itkConnectedThresholdImageFilter.h"

#include <itkHistogramMatchingImageFilter.h> // KEY_HE
/*Regina
  CleanUpOverlapArea is check overlapping area btw structures,
  usually btw adjacent structures, and
  choose the structure which has maximum ANN signal and
  surpress others below the given threshold value
 */
void         CleanUpOverlapArea(ProcessDescription & ANNXMLObject,
                                DataSet *curDataSet,
                                const std::string & fOutputDir,
                                const float MaskThresh,
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
  ProbabilityMapList *probabilityMaps =    ANNXMLObject.Get<ProbabilityMapList>(
      "ProbabilityMapList");
  int structure = -1;

  std::string ANNCutOut_name[OutputVectorSize];
  std::string StructureID[OutputVectorSize];

  RealImageType::Pointer ANNCutOut_Image[OutputVectorSize];
  const std::string      OutputDir = fOutputDir + "/";
  const std::string      ImageID( curDataSet->GetAttribute<StringValue>("Name") );

  typedef itk::AddImageFilter<RealImageType, RealImageType, RealImageType> SumImagesType;
  SumImagesType::Pointer Summation_ANNCutOut = SumImagesType::New();
  for( ProbabilityMapList::iterator pmi = probabilityMaps->begin();
       pmi != probabilityMaps->end();
       ++pmi )
    {
    structure++;

    ProbabilityMap *probabilityMap    = dynamic_cast<ProbabilityMap *>( pmi->second );

    StructureID[structure] = probabilityMap->GetAttribute<StringValue>("StructureID");

    ANNCutOut_name[structure] = OutputDir +  "ANNCutOut_median" + StructureID[structure] + ImageID + ".nii.gz";
    std::cout << " Read Image :: " << ANNCutOut_name[structure] << std::endl;
    try
      {
      ANNCutOut_Image[structure]     =
        itkUtil::ReadImage<RealImageType>(ANNCutOut_name[structure]);
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
  typedef itk::ImageRegionIterator<RealImageType> IteratorType;

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
      float max = 0;
      int   max_structure = -1;
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
      std::string WriteID(OutputDir +  "ANNCutOut_median_cleanup" + StructureID[strt] + ImageID + ".nii.gz");
      itkUtil::WriteImage<RealImageType>(ANNCutOut_Image[strt], WriteID);
      }

    typedef itk::Image<unsigned char, 3> BinaryMaskImageType;
    BinaryMaskImageType::Pointer outputMask;
      {
      itk::ThresholdImageFilter<RealImageType>::Pointer filterUpper =
        itk::ThresholdImageFilter<RealImageType>::New();
      filterUpper->SetInput(ANNCutOut_Image[strt]);
      filterUpper->SetOutsideValue(255);    // Foreground
      // filterUpper-> SetInsideValue(0); //Background
      filterUpper->ThresholdAbove(MaskThresh);
      filterUpper->Update();
      outputMask = itkUtil::TypeCast<RealImageType, BinaryMaskImageType>(
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
        const std::string ConnectedThresholdWriteID(
          OutputDir
          + "ANNCutOut_median_cleanup"
          + StructureID[strt]
          + ImageID
          +
          "_ConnectedThreshold.nii.gz");
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
        const std::string LabelWriteID(OutputDir + "ANNCutOut_median_cleanup"
                                       + StructureID[strt] + ImageID
                                       + "_labeled.nii.gz");
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
        const std::string MultiLabelWriteID(OutputDir + "ANNCutOut_median_cleanup"
                                            + StructureID[strt] + ImageID
                                            + "_MultiLabeled.nii.gz");
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
      // TODO:: Give Proper Name if it is given in the xml file (net
      // configuration file)
      std::string OneLabelWriteID_Mask(
        curDataSet->GetMaskFilenameByType( StructureID[strt] ) );

      if( (OneLabelWriteID_Mask == "") || (OneLabelWriteID_Mask == "na") )
        {
        std::cout << " MaskName is " << OneLabelWriteID_Mask << std::endl;
        OneLabelWriteID_Mask = OutputDir
          + "ANNCutOut_median_cleanup"
          + StructureID[strt]
          + ImageID
          + "_OneLabeled.nii.gz";
        }
      std::cout << "Resulting ANN masks for test image " << ImageID.c_str()
                << " write to disk as " << OneLabelWriteID_Mask << std::endl;
      itkUtil::WriteImage<BinaryMaskImageType>(OneLabeledMask,
                                               OneLabelWriteID_Mask);
      }
    }
}

static void CreateOutputMasks(ProcessDescription & ANNXMLObject,
                              DataSet *curDataSet,
                              const std::string & ProbabilityMapName,
                              const std::string & RhoMapName,
                              const std::string & PhiMapName,
                              const std::string & ThetaMapName,
                              const std::string & StructureID,
                              const std::string & AtlasType,
                              // const std::string & AtlasName,
                              // const std::string & AtlasLandmark,
                              // const float & DefGaussian,
                              DataSet::TypeVector & ImageTypeList,
                              const std::string & OutputDir,
                              const int InputVectorSize,
                              const int GradientProfileSize,
                              // const int IrisSize,
                              // const float CutOutGaussianValue,
                              const float CutOutThresh,
                              const float MaskThresh,
                              bool doANN,
                              const std::string & ANNModelFilename,
                              bool doSVM,
                              const std::string & SVMModelFilename,
                              const int OutputVectorSize,
                              int verbose,
                              const int structure,
                              std::map<std::string, ProbabilityMapImageType::Pointer> MapOfSubjects,
                              std::map<std::string, ImageLinearInterpolatorType::Pointer> MapOfSubjectsInterpolators,
                              ProbabilityMapImageType::Pointer DeformedSpherical[],
                              ProbabilityMapImageType::Pointer DeformedProbabilityMap[]
                              )
{
  //  const std::string
  // LandmarkImage(curDataSet->GetAtlasFilenameByType(landmarkType));
  if( verbose > 0 )
    {
    std::cout << "CreateOutput Masks with Structure number : "
              << structure
              << " ...."
              << std::endl;
    }

  typedef std::list<float>                    FloatList;
  typedef std::list<RealImageType::IndexType> FloatIndexList;
  if( doSVM )
    {
    if( !itksys::SystemTools::FileExists( ProbabilityMapName.c_str() ) )
      {
      std::cout << "Missing Probability Map File  " << ProbabilityMapName
                << std::endl;
      exit(-1);
      }
    }
  // get the registration ID to pick the registration object
  // in the current dataset.
  RegistrationParams *regParams =
    ANNXMLObject.Get<RegistrationParams>("RegistrationParams");
  const std::string regID( regParams->GetAttribute<StringValue>("ID") );

  const std::string       ImageID( curDataSet->GetAttribute<StringValue>("Name") );
  const RegistrationType *transform = curDataSet->GetRegistrationWithID(regID);
  if( transform == 0 )
    {
    std::cout << "ERROR Can not find registration for image object: "
              << ImageID << std::endl;
    exit(-1);
    }
  //  const std::string
  // landmarkType(transform->GetAttribute<StringValue>("LandmarkType"));
  const std::string landmarkType("");

  // Reading Map Of Images
  // std::map<std::string, ProbabilityMapImageType::Pointer> MapOfSubjects;
  // std::map<std::string,
  // ImageLinearInterpolatorType::Pointer> MapOfSubjectsInterpolators;
#if 0
  if( doSVM )
    {
    for( DataSet::TypeVector::iterator ImageList = ImageTypeList.begin();
         ImageList != ImageTypeList.end();
         ImageList++ )
      {
      const std::string
        curFilename( curDataSet->GetImageFilenameByType(*ImageList) );

      // HACK:  radius is hardcoded to size 2, it should be taken from command
      // line arguments.
      ProbabilityMapImageType::SizeType radius;
      radius[0] = 2; radius[1] = 2; radius[2] = 2;
      MapOfSubjects[*ImageList] =
        ReadMedianFilteredImage<ProbabilityMapImageType>(curFilename, radius);
      ImageLinearInterpolatorType::Pointer CurrInterpolator =
        ImageLinearInterpolatorType::New();
      CurrInterpolator->SetInputImage(MapOfSubjects[*ImageList]);
      MapOfSubjectsInterpolators[*ImageList] = CurrInterpolator;
      }
    }
#endif
  // Use first image in list to get spacing and size
  ProbabilityMapImageType::Pointer ReferenceImage =
    MapOfSubjects[*( ImageTypeList.begin() )];
  const ProbabilityMapImageType::SpacingType ImageSpacing =
    ReferenceImage->GetSpacing();
  const ProbabilityMapImageType::PointType ImageOrigin =
    ReferenceImage->GetOrigin();
  const ProbabilityMapImageType::SizeType ImageSize =
    ReferenceImage->GetLargestPossibleRegion()
    .GetSize();
  itk::ContinuousIndex<float, 3> ImageOuterBoundaryLocationInMM;
  for( int i = 0; i < 3; i++ )
    {
    ImageOuterBoundaryLocationInMM[i] = ImageSize[i] * ImageSpacing[i]
      + ImageOrigin[i];
    }

  const std::string AtlasToSubjRegistrationFilename
    ( transform->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename") );
  std::cout << "Deforming "
            << ProbabilityMapName
            << " into coordinate space of "
            << ImageID
            << " with "
            << AtlasToSubjRegistrationFilename     << std::endl;
  // Registrations must exists at this point, so apply them
  const bool TransformFileExists = itksys::SystemTools::FileExists(
      AtlasToSubjRegistrationFilename.c_str() );
  // If the transform file does not exists, then use the image information to
  // automatically generate it
  if( !TransformFileExists )  // If transform needs estimating, then do
                              // preprocessing steps to input images
    {
    std::string errmsg(AtlasToSubjRegistrationFilename);
    errmsg += " does not exist";
    throw  ProcessObjectException(errmsg);
    }

  ProbabilityMapList *probabilityMaps =
    ANNXMLObject.Get<ProbabilityMapList>("ProbabilityMapList");
  int               NumberOfProbabilityMaps = probabilityMaps->size();
  const std::string Def_promap_Name = AtlasToSubjRegistrationFilename + "def_probmap.nii.gz";
  const bool        DeformationExists =
    itksys::SystemTools::FileExists( Def_promap_Name.c_str() );
  if( verbose > 0 && !DeformationExists )
    {
    std::cout << "write out def_probmap.nii.gz" << std::endl
              << "::: " << AtlasToSubjRegistrationFilename << std::endl;

    itkUtil::WriteImage<ProbabilityMapImageType>(DeformedProbabilityMap[structure],
                                                 Def_promap_Name);
    }
  if( verbose > 0 )
    {
    std::cout << "Gradient Image Filter Apply...." << std::endl;
    }

  itk::GradientImageFilter<ProbabilityMapImageType, float,  float>::Pointer GradientFilter =
    itk::GradientImageFilter<ProbabilityMapImageType, float, float>::New();
  GradientFilter->SetInput(DeformedProbabilityMap[structure]);
  GradientFilter->Update();

  if( verbose > 0 )
    {
    std::cout << "Get Spherical Deformed Map...." << std::endl;
    }
  itk::Image<itk::CovariantVector<float, 3>, 3>::Pointer ProbMapGradient = GradientFilter->GetOutput();
  RealImageType::Pointer                                 DeformedRhoMap,  DeformedPhiMap, DeformedThetaMap;
  if( doSVM )
    {
    DeformedRhoMap = ImageWarper<RealImageType>(AtlasToSubjRegistrationFilename,
                                                RhoMapName, ReferenceImage);
    DeformedPhiMap = ImageWarper<RealImageType>(AtlasToSubjRegistrationFilename,
                                                PhiMapName, ReferenceImage);
    DeformedThetaMap = ImageWarper<RealImageType>(AtlasToSubjRegistrationFilename,
                                                  ThetaMapName, ReferenceImage);
    }

  if( doANN )
    {
    DeformedRhoMap = DeformedSpherical[0];
    DeformedPhiMap = DeformedSpherical[1];
    DeformedThetaMap = DeformedSpherical[2];
    }
  /*REGINA:: removed #if*/
  if( verbose > 0 )
    {
    std::cout << "Get Min and Max...." << std::endl;
    }
  itk::Index<3> min, max;
  DefineBoundingBox(DeformedProbabilityMap[structure], min, max);
  if( min[0] > max[0] || min[1] > max[1] || min[2] > max[2] )
    {
    std::cout
      << "ERROR IN BOUNDING BOX!  Min can not be greater than max!!!"
      << __FILE__ << " " << __LINE__ << std::endl;
    itk::MinimumMaximumImageCalculator<ProbabilityMapImageType>::Pointer
      minmaxcalc =
      itk::MinimumMaximumImageCalculator<ProbabilityMapImageType>::New();
    minmaxcalc->SetImage(DeformedProbabilityMap[structure]);
    minmaxcalc->Compute();
    std::cout
      << "Minimum Value is: "
      << static_cast<unsigned int>( minmaxcalc->GetMinimum() )
      << " Maximum Value is: "
      << static_cast<unsigned int>( minmaxcalc->GetMaximum() )
      << std::endl;
    std::cout << "Its min index position is : "
              << minmaxcalc->GetIndexOfMinimum() << std::endl;
    std::cout << "Its max index position is : "
              << minmaxcalc->GetIndexOfMaximum() << std::endl;
    min[0] = max[0];
    min[1] = max[1];
    min[2] = max[2];
    }

  ProbabilityMapImageType::RegionType BBregion;
  itk::Size<3>                        BBRange;
  itk::Index<3>                       BBIndex;
  unsigned long                       range[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2]};
  long int                            index[3] = { min[0], min[1], min[2]};
  BBRange.SetSize(range);
  BBregion.SetSize(BBRange);
  BBIndex.SetIndex(index);
  BBregion.SetIndex(BBIndex);
  std::string m_OutputDir = OutputDir + "/";
    {
    // / DO ANN VERSION
    if( doANN == true )
      {
      FloatIndexList ANNIndexList;
      FloatList      ANNOutputList;
      std::cout << "Applying ANN Model to segment " << StructureID.c_str()
                << " in " << ImageID.c_str()
                << " using model " <<  ANNModelFilename
                << "\nwith cut out tolerance range ( "
        // <<  (MINIMUM_VALUE
        // + PERCENT_MIN_MAX_TOLERANCE )  << ","
                << (0.0F + PERCENT_MIN_MAX_TOLERANCE)
                << ( HUNDRED_PERCENT_VALUE
           - PERCENT_MIN_MAX_TOLERANCE ) << ")" << std::endl;

      neural_scalar_type *input = new neural_scalar_type[InputVectorSize];
      neural_net_type *   ANN = new neural_net_type();
      ANN->load( ANNModelFilename.c_str() );
      // print out parameters...?

      const unsigned int region_size = BBregion.GetNumberOfPixels();
      float              current_percentage = 0.0F;
      unsigned int       voxels_processed = 0;
      // TODO:  The following for loop should be parallelized. Each pass
      // throught the loop is independant of all other passes.
      for( itk::ImageRegionIterator<ProbabilityMapImageType> bbri(
             DeformedProbabilityMap[structure], BBregion);
           !bbri.IsAtEnd(); ++bbri )
        {
        if( voxels_processed > region_size * current_percentage )
          {
          std::cout << current_percentage * 100.0F << "% Completed"
                    << std::endl;
          current_percentage += 0.05F;
          }
        voxels_processed++;
        const ProbabilityMapImageType::IndexType CurrentIndex = bbri.GetIndex();
        const float                              DeformedProbabilityMapPixel =
          static_cast<float>( DeformedProbabilityMap[structure]->GetPixel(CurrentIndex) );
        double Output = MINIMUM_VALUE;
        if( verbose > 8 )
          {
          std::cout << "Current Pixel's Deformed Probability Map::: "
                    << DeformedProbabilityMap[structure]->GetPixel(CurrentIndex)
                    << std::endl;
          }
        if( DeformedProbabilityMapPixel > ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) )
          {
//          unsigned long                       HundredPercentRegionSize = 0;
          if( DeformedProbabilityMapPixel <
              ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
            {
            std::vector<float> inputvector(InputVectorSize);
            AddInputVector(inputvector,
                           //  DeformedProbabilityMap[structure],
                           ProbMapGradient,
                           DeformedRhoMap,
                           DeformedPhiMap,
                           DeformedThetaMap,
                           //  ImageTypeList,
                           MapOfSubjects,
                           MapOfSubjectsInterpolators,
                           DeformedProbabilityMap,
                           NumberOfProbabilityMaps,
                           CurrentIndex,
                           GradientProfileSize/*,
              IrisSize*/   );
              {
              std::vector<float>::const_iterator fi = inputvector.begin();
              for( int j = 0; j < InputVectorSize; j++, fi++ )
                {
                input[j] = *fi;
                }
              }
            if( verbose > 3 )
              {
              std::cout << "\n-- " << CurrentIndex << " -- I: ";
              for( int z = 0; z < InputVectorSize; z++ )
                {
                std::cout << input[z] << " ";
                }
              }
            // float CvANN_MLP::predict( const CvMat* _inputs, CvMat* _outputs )
            // output can be several from neural net
            neural_vector_type current_input = cvCreateMat( 1,
                                                            InputVectorSize,
                                                            CV_64FC1);
            cvInitMatHeader( current_input,
                             1,
                             InputVectorSize,
                             CV_64FC1,
                             input);
            neural_vector_type current_output = cvCreateMat( 1,
                                                             OutputVectorSize,
                                                             CV_64FC1);

            ANN->predict( current_input, current_output );
            // Output = HUNDRED_PERCENT_VALUE
            //               * static_cast< float >(
            // current_output->data.db[structure]);
            Output = HUNDRED_PERCENT_VALUE
              * static_cast<float>(
                CV_MAT_ELEM( *current_output, double, 0, structure) );
            }
          else
            {
//            HundredPercentRegionSize++;
            Output = HUNDRED_PERCENT_VALUE;
            }
          // Comments
          if( verbose > 3 )
            {
            std::cout << "O: " << Output << std::endl;
            }
          // ANNIndexList.push_back(CurrentIndex);
          // ANNOutputList.push_back(Output);
          // Dealing with boundary values
          if( Output > ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
            {
            ANNIndexList.push_back(CurrentIndex);
            ANNOutputList.push_back(HUNDRED_PERCENT_VALUE);
            }
          else if( Output > MINIMUM_VALUE + PERCENT_MIN_MAX_TOLERANCE )
            {
            ANNIndexList.push_back(CurrentIndex);
            ANNOutputList.push_back(Output);
            }
          else // if ( Output <= MINIMUM_VALUE + PERCENT_MIN_MAX_TOLERANCE )
            {
            ANNIndexList.push_back(CurrentIndex);
            ANNOutputList.push_back( MINIMUM_VALUE );
            }
          }
        else // if !{  DeformedProbabilityMapPixel > ( 0.0F +
             // PERCENT_MIN_MAX_TOLERANCE ) }
          {
          Output = MINIMUM_VALUE;
          ANNIndexList.push_back(CurrentIndex);
          ANNOutputList.push_back( Output );
          // Output = 0.0F;
          //          std::cout << '#';
          }
        }

      // std::string temp =ANNModelFilename + "Predict";
      // ANN->save( temp.c_str() );
      delete[] input;
        {
        RealImageType::Pointer ANNCutOut = RealImageType::New();
          {
          ANNCutOut->CopyInformation(MapOfSubjects[AtlasType.c_str()]);
          ANNCutOut->SetRegions( MapOfSubjects[AtlasType.c_str()]->GetLargestPossibleRegion() );
          ANNCutOut->Allocate();
          ANNCutOut->FillBuffer( MINIMUM_VALUE);
          }
        FloatList::const_iterator oiter = ANNOutputList.begin();
        for( FloatIndexList::const_iterator iiter = ANNIndexList.begin();
             iiter != ANNIndexList.end(); iiter++, oiter++ )
          {
          ANNCutOut->SetPixel(*iiter, *oiter);
          }
        std::string WriteID(m_OutputDir + "ANNCutOut" + StructureID + ImageID + ".nii.gz");
          {
          std::string destination_dir = itksys::SystemTools::GetFilenamePath(
              WriteID);
          itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
          }
        /*
         * Write out Deformed Probability Map only when Debugging
         */
        if( verbose > 5 )
          {
          itkUtil::WriteImage<RealImageType>(ANNCutOut, WriteID);
          std::string DefMapWriteID(m_OutputDir + "ANNCutOut" + StructureID
                                    + ImageID + "_defprobmap.nii.gz");
          itkUtil::WriteImage<RealImageType>(DeformedProbabilityMap[structure],
                                             DefMapWriteID);
          }
        // =======================================================
        // Median Filter, and then pull out the largest connected region to
        // write the final mask above the 50% level
        // =======================================================
        typedef itk::MedianImageFilter<RealImageType,
                                       RealImageType> MedianFilterType;
        MedianFilterType::Pointer medianFilter = MedianFilterType::New();
        RealImageType::SizeType   AllTheSame;
        AllTheSame.Fill(1);
        medianFilter->SetRadius(AllTheSame);
        medianFilter->SetInput(ANNCutOut);
        medianFilter->Update();

        RealImageType::Pointer medianImage = medianFilter->GetOutput();
        WriteID = m_OutputDir + "ANNCutOut_median" + StructureID + ImageID + ".nii.gz";
        itkUtil::WriteImage<RealImageType>(medianImage, WriteID);
        }
      }

      {
      // / DO SVM VERSION
      if( doSVM )
        {
        std::list<ProbabilityMapImageType::IndexType> SVMIndexList;
        std::list<ProbabilityMapImageType::PixelType> SVMOutputList;
        std::cout << "Applying SVM Model to segment " << StructureID.c_str()
                  << " with image " << ImageID.c_str()
                  << " using model " <<  SVMModelFilename
                  << " with cut out threshold of " << CutOutThresh <<  std::endl;
        struct svm_model *SVM = svm_load_model( SVMModelFilename.c_str() );
        if( SVM == NULL )
          {
          std::cout << "ERROR:  Can not read model file: " << SVMModelFilename
                    << std::endl;
          exit(-1);
          }
        double *prob_estimates = (double *)malloc( 2 * sizeof( double ) );
          {
          struct svm_node *x = (struct svm_node *)malloc(
              ( InputVectorSize + 1 ) * sizeof( struct svm_node ) );
          x[InputVectorSize].index = -1;
          for( itk::ImageRegionIterator<ProbabilityMapImageType> bbri(
                 DeformedProbabilityMap[structure], BBregion);
               !bbri.IsAtEnd(); ++bbri )
            {
            const ProbabilityMapImageType::IndexType CurrentIndex = bbri.GetIndex();
            const float                              DeformedProbabilityMapPixel =
              static_cast<float>( DeformedProbabilityMap[structure]->GetPixel(CurrentIndex) );
            if( DeformedProbabilityMapPixel > CutOutThresh )
              {
              std::cout << "InputvectorSize:: " << InputVectorSize << std::endl
                        << DeformedProbabilityMap[structure]->GetPixel(CurrentIndex) << std::endl
                        << DeformedRhoMap->GetPixel(CurrentIndex) << std::endl
                        << DeformedPhiMap->GetPixel(CurrentIndex) << std::endl;

              std::vector<float> inputvector(InputVectorSize);

              AddInputVector(inputvector,
                             //  DeformedProbabilityMap[structure],
                             ProbMapGradient,
                             DeformedRhoMap,
                             DeformedPhiMap,
                             DeformedThetaMap,
                             //  ImageTypeList,
                             MapOfSubjects,
                             MapOfSubjectsInterpolators,
                             DeformedProbabilityMap,
                             NumberOfProbabilityMaps,
                             CurrentIndex,
                             GradientProfileSize/*,
              IrisSize*/     );
                {
                std::vector<float>::const_iterator fi = inputvector.begin();
                for( int a = 0; a < InputVectorSize; a++, fi++ )
                  {
                  x[a].index = a + 1;
                  x[a].value = *fi;
                  }
                svm_predict_probability(SVM, x, prob_estimates);
                float Output = ( HUNDRED_PERCENT_VALUE * prob_estimates[1] );
                if( Output > ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
                  {
                  SVMIndexList.push_back(CurrentIndex);
                  SVMOutputList.push_back( static_cast<ProbabilityMapImageType::
                                                       PixelType>(
                                             HUNDRED_PERCENT_VALUE ) );
                  }
                else if( Output > 0.5F )
                  {
                  SVMIndexList.push_back(CurrentIndex);
                  SVMOutputList.push_back( static_cast<ProbabilityMapImageType::
                                                       PixelType>( Output ) );
                  }
                }
              }
            else if( DeformedProbabilityMapPixel >= ( 1.0F - CutOutThresh ) )
              {
              SVMIndexList.push_back(CurrentIndex);
              SVMOutputList.push_back( static_cast<ProbabilityMapImageType::
                                                   PixelType>(
                                         HUNDRED_PERCENT_VALUE ) );
              }
            }
          free(x);
          }
        svm_destroy_model(SVM);
          {
          ProbabilityMapImageType::Pointer SVMCutOut;
          CreateNewImageFromTemplate(SVMCutOut, MapOfSubjects[AtlasType.c_str()]);
            {
            std::list<ProbabilityMapImageType::PixelType>::const_iterator oiter =
              SVMOutputList.begin();
            for( std::list<ProbabilityMapImageType::IndexType>::const_iterator
                 iiter = SVMIndexList.begin();
                 iiter != SVMIndexList.end();
                 iiter++, oiter++ )
              {
              SVMCutOut->SetPixel(*iiter, *oiter);
              }
            }
          std::string WriteID(m_OutputDir + "SVMCutOut" + StructureID + ImageID + ".nii.gz");
            {
            std::string destination_dir = itksys::SystemTools::GetFilenamePath(
                WriteID);
            itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
            }
          itkUtil::WriteImage<ProbabilityMapImageType>( SVMCutOut, WriteID.c_str() );
          SVMCutOut = SmoothSingleImage<ProbabilityMapImageType>(
              SVMCutOut,
              0.5 );
          SVMCutOut->SetMetaDataDictionary(
            MapOfSubjects[AtlasType.c_str()]->GetMetaDataDictionary() );
          itk::ThresholdImageFilter<ProbabilityMapImageType>::Pointer filterLower =
            itk::ThresholdImageFilter<ProbabilityMapImageType>::New();
          filterLower->SetInput(SVMCutOut);
          filterLower->SetOutsideValue(0); // Background Value
          filterLower->ThresholdBelow(MaskThresh);
          filterLower->Update();
          SVMCutOut = filterLower->GetOutput();
          SVMCutOut->SetMetaDataDictionary(
            MapOfSubjects[AtlasType.c_str()]->GetMetaDataDictionary() );
          itk::ThresholdImageFilter<ProbabilityMapImageType>::Pointer filterUpper =
            itk::ThresholdImageFilter<ProbabilityMapImageType>::New();
          filterUpper->SetInput(SVMCutOut);
          filterUpper->SetOutsideValue(1); // Forground Value
          filterUpper->ThresholdAbove(MaskThresh);
          filterUpper->Update();
          WriteID = m_OutputDir + "SVMThresh" + StructureID + ImageID + ".mask";
            {
            std::string destination_dir = itksys::SystemTools::GetFilenamePath(
                WriteID);
            itksys::SystemTools::MakeDirectory( destination_dir.c_str() );
            }
          typedef itk::Image<unsigned char, 3> BinaryMaskImageType;
          BinaryMaskImageType::Pointer outputMask =
            itkUtil::TypeCast<ProbabilityMapImageType, BinaryMaskImageType>(
              filterUpper->GetOutput() );
          itkUtil::WriteImage<BinaryMaskImageType>(outputMask, WriteID);
          std::cout << "Resulting SVM masks for test image " << ImageID
                    << " written to disk as " << WriteID << std::endl;
          }
        }
      }
    }
}

// /Apply Model
/** detailed description of cutout
 * \ingroup AM
 */
int ApplyModel(ProcessDescription & ANNXMLObject,
               int verbose,
               bool histogramEqualization,
               bool multiStructureThreshold,
               bool doTest)
{
  // find out the registration parameters

  RegistrationParams *regParams =
    ANNXMLObject.Get<RegistrationParams>("RegistrationParams");
  const std::string regID( regParams->GetAttribute<StringValue>("ID") );
  const std::string imageTypeToUse( regParams->GetAttribute<StringValue>(
                                      "ImageTypeToUse") );
  std::string ANNModelFilename, SVMModelFilename;
  // Get the atlas dataset
  DataSet *         atlasDataSet = ANNXMLObject.GetAtlasDataSet();
  const std::string atlasImage( atlasDataSet->GetImageFilenameByType(
                                  imageTypeToUse.c_str() ) );

  NeuralParams *model = ANNXMLObject.Get<NeuralParams>("NeuralNetParams");

  if( model == 0 )
    {
    std::cout << "Model information (gaussian size, iris size, "
      "mask smoothing value, training vector filename) not given,"
      " so training vectors cannot be created." << std::endl;
    exit(-1);
    }
  std::cout << "** Using Model file: "
            << model->GetAttribute<StringValue>("TrainingModelFilename") <<  std::endl;

  bool            doANN = 1;
  int             Iterations = 0;
  ApplyModelType *applyModel =
    ANNXMLObject.Get<ApplyModelType>("ApplyModel");
  if( applyModel == 0 )
    {
    std::cout << "Cut Out Information (Gaussian, threshold, "
              << " write option, and write directory) not given,"
              << " so output cannot be computed."
              << std::endl;
    exit(-1);
    }
  ANNParams *annParams = ANNXMLObject.Get<ANNParams>("ANNParams");
  if( annParams == 0 )
    {
    doANN = 0;
    }
  else
    {
    Iterations      = annParams->GetAttribute<IntValue>("Iterations");
    // XXXXXX

    if(  applyModel->GetAttributeIfExist<IntValue>("IterationNumber") != 0
         && applyModel->GetAttributeIfExist<IntValue>("IterationNumber") < Iterations + 1 )
      {
      Iterations = applyModel->GetAttributeIfExist<IntValue>("IterationNumber");
      }
    // Iteration number should be optimal one if exist
    ANNModelFilename = model->GetAttribute<StringValue>("TrainingModelFilename");
    ANNModelFilename += "ANN";
    if( doTest )
      {
      std::string ANNTestVectorFilename =
        model->GetAttribute<StringValue>("TestVectorFileName");

      ANNTestVectorFilename += "ANN.hdr";
      std::fstream ANNHeaderTestVectorStream;
      ANNHeaderTestVectorStream.open( ANNTestVectorFilename.c_str(),
                                      std::ios::in
                                      | std::ios::binary );

      if( !ANNHeaderTestVectorStream.is_open() )
        {
        std::cout
          << "Error: Writing Header File: \
                      Could not open ANN Test vector file name of \""
          << ANNTestVectorFilename.c_str()
          << "\""
          << std::endl;
        exit(-1); // break?
        }
      char currentline[MAX_LINE_SIZE];
      while( ANNHeaderTestVectorStream.getline(currentline,
                                               MAX_LINE_SIZE - 1) )
        {
        std::string        temp;
        std::istringstream iss(currentline, std::istringstream::in);
        iss >> temp;
        if( temp == ANNModelFilename )
          {
          iss >> Iterations;
          std::cout << "================================================="
                    << std::endl
                    << "** Using optimally trained point "
                    << Iterations << std::endl;
          std::cout << "================================================="
                    << std::endl;
          }
        }

      ANNHeaderTestVectorStream.close();
      }
    // TODO Insert the GetAttribute of Iteration from ApplyModel
    }
  bool       doSVM = 1;
  SVMParams *svModel = ANNXMLObject.Get<SVMParams>("SVMParams");
  if( svModel == 0 )
    {
    doSVM = 0;
    }

  if( ( !doANN ) && ( !doSVM ) )
    {
    std::cout << "No models (ANN, SVM) chosen" << std::cout;
    return -1;
    }

  const float MaskThresh = applyModel->GetAttribute<FloatValue>("MaskThresh");
  std::cout << "Mask Thrshold Value :: " << MaskThresh << "\n";
  if( ( MaskThresh < MINIMUM_VALUE + PERCENT_MIN_MAX_TOLERANCE )
      || ( MaskThresh > ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) ) )
    {
    std::cout << "ERROR: Invalid mask threshold given as " << MaskThresh
              << std::endl;
    std::cout << "Valid values are in range : [" << MINIMUM_VALUE
      + PERCENT_MIN_MAX_TOLERANCE << "," << HUNDRED_PERCENT_VALUE
      - PERCENT_MIN_MAX_TOLERANCE << "]" << std::endl;
    exit(-1);
    }

  DataSet::TypeVector ImageTypeList =
    ANNXMLObject.GetAtlasDataSet()->ImageTypes();
  int NumberOfImageTypes = ImageTypeList.size();

  if( !NumberOfImageTypes )
    {
    std::cout << "No images types found (ApplyModel). Cannot compute neural net output."
              << std::endl;
    exit(-1);
    }
  const std::string landmarkType("");
  const std::string SubjectLandmark("");
  const std::string AtlasLandmark("");
  /*=======================================================================*/
  /* Scan Each Probability Map and Apply Trained Model to the Target*/
  if( doANN )
    {
    ProbabilityMapList *probabilityMaps =
      ANNXMLObject.Get<ProbabilityMapList>("ProbabilityMapList");
    int          OutputVectorSize = probabilityMaps->size();
    unsigned int inputVectorSizeRequirement =
      InputVectorSizeRequirement(OutputVectorSize,
                                 NumberOfImageTypes,
                                 model->GetAttribute<IntValue>("GradientProfileSize")/*,
                                                   model->GetAttribute<IntValue>("IrisSize")*/);
    if( inputVectorSizeRequirement == 0 )
      {
      std::cout << "Input vector size must not be zero!" << std::endl;
      }
    /* Scan Each Target and Apply Trained Model to the Target*/
    ProcessDescription::TrainDataSetListType
                                                       subjectDataSets( ANNXMLObject.GetApplyDataSets() );
    ProcessDescription::TrainDataSetListType::iterator dsIt;

    if( subjectDataSets.size() == 0 )
      {
      std::cerr << "No \"Apply\" datasets found" << std::endl;
      exit(-1);
      }
    /*=============================  Scan each subject to apply*/
    char temp[10];
    sprintf(temp, "%09d", Iterations);
    ANNModelFilename += temp;
    if( !itksys::SystemTools::FileExists( ANNModelFilename.c_str() ) )
      {
      std::cout << "ANN Model File " << ANNModelFilename
                << " does not exist, so output vectors cannot be computed."
                << std::endl;
      exit(-1);
      }
    std::cout << "================================================="
              << std::endl
              << "** Using Specified trained point "
              << Iterations << std::endl;
    std::cout << "================================================="
              << std::endl;
    for( dsIt = subjectDataSets.begin();
         dsIt != subjectDataSets.end();
         ++dsIt )
      {
      if( !multiStructureThreshold )
        {
        const std::string       ImageID( ( *dsIt )->GetAttribute<StringValue>("Name") );
        const RegistrationType *transform = ( *dsIt )->GetRegistrationWithID(regID);
        if( transform == 0 )
          {
          std::cout << "ERROR Can not find registration for image object: "
                    << ImageID << std::endl;
          exit(-1);
          }
        std::cout << "Starting Registrations." << std::endl;

        // TODO:  numThreads should be exposed on the command line.
        const unsigned int numThreads = 1;
        // itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
        GenerateRegistrations(ANNXMLObject, true, true, numThreads);
        // Reading Subject Images as a Map
        std::map<std::string, ProbabilityMapImageType::Pointer> MapOfSubjects;
        std::map<std::string, ProbabilityMapImageType::Pointer> MapOfAtlas;
        std::map<std::string, ImageLinearInterpolatorType::Pointer>
        MapOfSubjectsInterpolators;
        for( DataSet::TypeVector::iterator ImageList = ImageTypeList.begin();
             ImageList != ImageTypeList.end();
             ++ImageList )
          {
          const std::string
            curFilename( ( *dsIt )->GetImageFilenameByType(*ImageList) );
          if( verbose > 3 )
            {
            std::cout << "current file name: " << curFilename << std::endl;
            }
          ProbabilityMapImageType::SizeType radius;
          radius[0] = 0; radius[1] = 0; radius[2] = 0;
          MapOfSubjects[*ImageList] =
            ReadMedianFilteredImage<ProbabilityMapImageType>(curFilename, radius);

          ImageLinearInterpolatorType::Pointer CurrInterpolator =
            ImageLinearInterpolatorType::New();
          CurrInterpolator->SetInputImage(MapOfSubjects[*ImageList]);
          MapOfSubjectsInterpolators[*ImageList] = CurrInterpolator;
          // Histogram Equalization Option. KEY_HE
          if( histogramEqualization )
            {
            std::cout << "- Histogram Equalization Performs on the Images " << *ImageList
                      << std::endl;
            std::string atlasVolumeFileName = atlasDataSet->
              GetImageFilenameByType( *ImageList);
            // Check if altlas image is given for the type
            std::transform( atlasVolumeFileName.begin(),
                            atlasVolumeFileName.end(),
                            atlasVolumeFileName.begin(), ::tolower);

            if( atlasVolumeFileName == "na" )
              {
              std::cout << "  !!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!! " << std::endl
                        << "  -  If atlas data set of the " << *ImageList << std::endl
                        << "     is not given, Histogram Equalization Process will be skipped "
                        << std::endl
                        << "  !!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!! " << std::endl;
              break;
              }
            MapOfAtlas[*ImageList] =  ReadMedianFilteredImage<ProbabilityMapImageType>(
                atlasDataSet->GetImageFilenameByType( *ImageList), radius );
            typedef itk::HistogramMatchingImageFilter<ProbabilityMapImageType,
                                                      ProbabilityMapImageType> HEFilterType;
            HEFilterType::Pointer HEFilter = HEFilterType::New();
            HEFilter->SetReferenceImage( MapOfAtlas[*ImageList] );
            HEFilter->SetInput( MapOfSubjects[*ImageList] );
            HEFilter->SetNumberOfHistogramLevels( 100);
            HEFilter->SetNumberOfMatchPoints( 15);
            HEFilter->ThresholdAtMeanIntensityOn();   // TOOD:: Check if we want
                                                      // to do this.
            HEFilter->Update();
            // overwrite the input image to the Histogram Equalized one.
            MapOfSubjects[*ImageList] = HEFilter->GetOutput();
            }
          }
        ProbabilityMapImageType::Pointer DeformedSpherical[3];
        ProbabilityMapImageType::Pointer DeformedProbabilityMap[OutputVectorSize];

        ProbabilityMapImageType::Pointer ReferenceImage =
          MapOfSubjects[*( ImageTypeList.begin() )];

        const ProbabilityMapImageType::SpacingType ImageSpacing   =
          ReferenceImage->GetSpacing();
        const ProbabilityMapImageType::PointType ImageOrigin      =
          ReferenceImage->GetOrigin();
        const ProbabilityMapImageType::SizeType ImageSize =
          ReferenceImage->GetLargestPossibleRegion().GetSize();
        const std::string AtlasToSubjRegistrationFilename
          ( transform->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename") );
        // Registrations must exists at this point, so apply them
        const bool TransformFileExists =
          itksys::SystemTools::FileExists( AtlasToSubjRegistrationFilename.c_str() );
        // If the transform file does not exists, then use the image information
        // to
        // automatically generate it

        if( !TransformFileExists )  // If transform needs estimating, then do
                                    // preprocessing steps to input images
          {
          std::string errmsg(AtlasToSubjRegistrationFilename);
          errmsg += " does not exist";
          throw  ProcessObjectException(errmsg);
          }
        int tempIndex = 0;
        for(  ProbabilityMapList::iterator pmi = probabilityMaps->begin();
              pmi != probabilityMaps->end();
              ++pmi,  ++tempIndex )
          {
          ProbabilityMap *probabilityMap      = dynamic_cast<ProbabilityMap *>
            ( pmi->second );
          try
            {
            DeformedProbabilityMap[tempIndex] = ImageWarper<ProbabilityMapImageType>
                (AtlasToSubjRegistrationFilename,
                probabilityMap->
                GetAttribute<StringValue>("Filename"),
                ReferenceImage);
            }
          catch( ProcessObjectException & ex )
            {
            std::string errmsg(probabilityMap->GetAttribute<StringValue>("Filename") );
            errmsg += " does not exist";
            std::cout << ex.Error();
            // throw ProcessObjectException(errmsg);
            }
          catch( ... )
            {
            std::cout << " Caught a unknown exception "
                      << probabilityMap->GetAttribute<StringValue>("Filename")
                      << " has a problem ( Check the validity of the file please) "
                      << std::endl;

            exit(-1);
            }

          if( verbose > 0 )
            {
            std::cout << "*** Read Deformed Prob Map..." << tempIndex << std::endl
                      << "  => " << probabilityMap->GetAttribute<StringValue>("Filename")
                      << std::endl;
            }
          if( DeformedProbabilityMap[tempIndex].IsNull() )
            {
            std::cout << "No transformation " << AtlasToSubjRegistrationFilename
                      << " filename found. Cannot compute registration." << std::endl;
            exit(-1);
            }
          if( pmi ==  probabilityMaps->begin() )
            {
            DeformedSpherical[0] = ImageWarper<ProbabilityMapImageType>
                (AtlasToSubjRegistrationFilename,
                probabilityMap->GetAttribute<StringValue>(
                  "rho"),
                ReferenceImage);
            DeformedSpherical[1] = ImageWarper<ProbabilityMapImageType>
                (AtlasToSubjRegistrationFilename,
                probabilityMap->GetAttribute<StringValue>(
                  "phi"),
                ReferenceImage);
            DeformedSpherical[2] = ImageWarper<ProbabilityMapImageType>
                (AtlasToSubjRegistrationFilename,
                probabilityMap->GetAttribute<StringValue>(
                  "theta"),
                ReferenceImage);
            }
          }

        int structure = 0;
        for( ProbabilityMapList::iterator pmi = probabilityMaps->begin();
             pmi != probabilityMaps->end();
             ++pmi, ++structure )
          {
          ProbabilityMap *probabilityMap
            = dynamic_cast<ProbabilityMap *>( pmi->second );

          CreateOutputMasks
            (ANNXMLObject,
            ( *dsIt ),
            probabilityMap->GetAttribute<StringValue>("Filename"),    //
                                                                      // ProbabiltyMapNameInput
            probabilityMap->GetAttribute<StringValue>("rho"),         // rho
                                                                      // file
                                                                      // name
            probabilityMap->GetAttribute<StringValue>("phi"),         // phi
                                                                      // file
                                                                      // name
            probabilityMap->GetAttribute<StringValue>("theta"),       // theta
                                                                      // file
                                                                      // name
            probabilityMap->GetAttribute<StringValue>("StructureID"), //
                                                                      // StructureID
            imageTypeToUse,                                           //
                                                                      // AtlasType
                                                                      // atlasImage,
                                                                      //      //
                                                                      // AtlasName
                                                                      // AtlasLandmark,
                                                                      // AtlasLandmark
                                                                      // probabilityMap->GetAttribute<FloatValue>("Gaussian"),
                                                                      //     //
                                                                      // DefGaussian
            ImageTypeList,                                            //
                                                                      // ImageTypeList
            ( *dsIt )->GetAttribute<StringValue>("OutputDir"),        //
                                                                      // OutputDir
            inputVectorSizeRequirement,
            model->GetAttribute<IntValue>("GradientProfileSize"),
            applyModel->GetAttribute<FloatValue>("CutOutThresh"), //
                                                                  // CutOutThresh
            MaskThresh,                                           //
                                                                  // MaskThresh
            doANN,
            ANNModelFilename,
            doSVM,
            SVMModelFilename,
            OutputVectorSize,
            verbose,
            structure,
            MapOfSubjects,
            MapOfSubjectsInterpolators,
            DeformedSpherical,
            DeformedProbabilityMap
            );
          }
        }
      CleanUpOverlapArea(ANNXMLObject,
                         ( *dsIt ),
                         ( *dsIt )->GetAttribute<StringValue>("OutputDir"),
                         MaskThresh,
                         OutputVectorSize,
                         verbose);
      } // End of Subject Scanning
    }   // end of doANN

#if 0
  if( doSVM )
    {
    SVMModelFilename = model->GetAttribute<StringValue>("Filename");
    SVMModelFilename += "SVM";
    SVMModelFilename += probabilityMap->GetAttribute<StringValue>(
        "StructureID");
    if( !itksys::SystemTools::FileExists( SVMModelFilename.c_str() ) )
      {
      std::cout
        <<
        "SVM Model File does not exist, so output vectors cannot be computed."
        << std::endl;
      exit(-1);
      }

    CreateOutputMasks
      (ANNXMLObject,
      ( *dsIt ),
      probabilityMap->GetAttribute<StringValue>("Filename"),
      probabilityMap->GetAttribute<StringValue>("rho"),
      probabilityMap->GetAttribute<StringValue>("phi"),
      probabilityMap->GetAttribute<StringValue>("theta"),
      probabilityMap->GetAttribute<StringValue>("StructureID"),
      imageTypeToUse,
      atlasImage,
      AtlasLandmark,
      probabilityMap->GetAttribute<FloatValue>("Gaussian"),
      std::map<std::string, ProbabilityMapImageType::Pointer> MapOfSubjects,
      applyModel->GetAttribute<StringValue>("OutputDir"),
      inputVectorSizeRequirement,
      model->GetAttribute<IntValue>("GradientProfileSize"),
      /* model->GetAttribute<IntValue>("IrisSize"), */
      applyModel->GetAttribute<FloatValue>("CutOutGaussian"),
      applyModel->GetAttribute<FloatValue>("CutOutThresh"),
      MaskThresh,
      doANN,
      ANNModelFilename,
      doSVM,
      SVMModelFilename,
      OutputVectorSize,
      verbose);
    }
#endif

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
               bool multiStructureThreshold,
               bool doTest)
{
  ProcessDescription ANNXMLObject;

  try
    {
    ReadXML(XMLFile.c_str(), ANNXMLObject);
    }
  catch( ProcessObjectException & ex )
    {
    std::cerr << " this error possibly caused by no present of file\n";
    std::cerr << ex.Error() << std::endl;
    exit(-1);
    }
  if( ANNXMLObject.Verify() != true )
    {
    std::cerr << "XML file " << " is invalid." << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    ANNXMLObject.PrintSelf(std::cerr, 0);
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    std::cerr << "FULL DUMP ===============================" << std::endl;
    return -1;
    }
  int rval;
  try
    {
    rval = ApplyModel(ANNXMLObject,
                      verbose,
                      histogramEqualization,
                      multiStructureThreshold,
                      doTest);
    }
  catch( ProcessObjectException & ex )
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
