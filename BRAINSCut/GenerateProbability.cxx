#include <NetConfiguration.h>
#include <itksys/SystemTools.hxx>
#include <list>
#include <vector>
#include <iostream>
#include <itkIO.h>
#include "Utilities.h"
#include "NetConfigurationParser.h"

static InternalImageType::Pointer
GaussianSmooth(const InternalImageType::Pointer image,
               const double Sigma)
{
  InternalImageType::SpacingType imageSpacing = image->GetSpacing();

  for( int currDim = 0; currDim < InternalImageType::ImageDimension; currDim++ )
    {
    if( Sigma < imageSpacing[currDim] )
      {
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:  Smoothing with radius (" << Sigma
                << ") less than pixel size " <<  imageSpacing << std::endl;
      std::cout << "WARNING:  This is likely not what you wanted." << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      std::cout << "WARNING:" << std::endl;
      }
    }

  g_GaussianFilterType::Pointer g_gaussianFilterX = g_GaussianFilterType::New();
  g_gaussianFilterX->SetNormalizeAcrossScale(true);
  g_gaussianFilterX->SetSigma(Sigma);
  g_gaussianFilterX->SetDirection(0);
  g_gaussianFilterX->SetInput(image);
  g_gaussianFilterX->ReleaseDataFlagOn();
  g_GaussianFilterType::Pointer g_gaussianFilterY = g_GaussianFilterType::New();
  g_gaussianFilterY->SetNormalizeAcrossScale(true);
  g_gaussianFilterY->SetSigma(Sigma);
  g_gaussianFilterY->SetDirection(1);
  g_gaussianFilterY->SetInput( g_gaussianFilterX->GetOutput() );
  g_gaussianFilterY->ReleaseDataFlagOn();
  g_GaussianFilterType::Pointer g_gaussianFilterZ = g_GaussianFilterType::New();
  g_gaussianFilterZ->SetNormalizeAcrossScale(true);
  g_gaussianFilterZ->SetSigma(Sigma);
  g_gaussianFilterZ->SetDirection(2);
  g_gaussianFilterZ->SetInput( g_gaussianFilterY->GetOutput() );
  g_gaussianFilterZ->ReleaseDataFlagOn();
  g_gaussianFilterZ->Update();
  InternalImageType::Pointer gaussianImage = g_gaussianFilterZ->GetOutput();
  // NOTE:  Sometimes the "SetNormalizeAcrossScale" will make the zeros slightly
  // above zero.
  InternalImageType::Pointer returnSmoothed =
    itkUtil::ScaleAndCast<InternalImageType, InternalImageType>(
      gaussianImage,
      0,
      HUNDRED_PERCENT_VALUE);
  return returnSmoothed;
}

//
// this is to remove warnings about assigning a float to a long int.
// I could have just added a cast, but it would make hairy code even uglier,
inline InternalImageType::IndexType::IndexValueType
TruncatedHalf(const InternalImageType::SizeType::SizeValueType & v)
{
  return static_cast<InternalImageType::IndexType::IndexValueType>
         ( static_cast<double>( v ) * 0.5 );
}

// Generate PolarImages about the center of the Atlas Image.
static void GeneratePolarImages(DataSet *atlasDataSet, /*const
                                                                 * std::string
                                                                 *
                                                                 *&ANNConfigurationFilename,*/
                                const InternalImageType::Pointer AtlasImage)
{
  InternalImageType::IndexType CenterIndex;

    {
    const InternalImageType::SizeType
      size( AtlasImage->GetLargestPossibleRegion().GetSize() );
    CenterIndex[0] = TruncatedHalf(size[0]);
    CenterIndex[1] = TruncatedHalf(size[1]);
    CenterIndex[2] = TruncatedHalf(size[2]);
    }
  itk::Point<float, 3> ImageCenterLocationInMM;

  AtlasImage->TransformIndexToPhysicalPoint(CenterIndex,
                                            ImageCenterLocationInMM);
  // std::cout << AtlasImage << std::endl;
  std::cout << "CENTER LOC " << ImageCenterLocationInMM << std::endl;

  InternalImageType::Pointer rhoImage, phiImage, thetaImage;
  CreateNewFloatImageFromTemplate(rhoImage, AtlasImage);
  CreateNewFloatImageFromTemplate(phiImage, AtlasImage);
  CreateNewFloatImageFromTemplate(thetaImage, AtlasImage);

  itk::ImageRegionIterator<InternalImageType> it(
    AtlasImage, AtlasImage->GetLargestPossibleRegion() );
  it.GoToBegin();
  itk::ImageRegionIterator<InternalImageType> rhoit(
    rhoImage, rhoImage->GetLargestPossibleRegion() );
  rhoit.GoToBegin();
  itk::ImageRegionIterator<InternalImageType> phiit(
    phiImage, phiImage->GetLargestPossibleRegion() );
  phiit.GoToBegin();
  itk::ImageRegionIterator<InternalImageType> thetait(
    thetaImage, thetaImage->GetLargestPossibleRegion() );
  thetait.GoToBegin();
  itk::Point<float, 3> CurrentLocationInMM;
  itk::Point<float, 3> LocationWithRespectToCenterOfImageInMM;
  while( !it.IsAtEnd() )
    {
    const InternalImageType::IndexType CurrentIndex = it.GetIndex();
    AtlasImage->TransformIndexToPhysicalPoint(CurrentIndex, CurrentLocationInMM);
    for( unsigned i = 0; i < 3; i++ )
      {
      LocationWithRespectToCenterOfImageInMM[i] =
        CurrentLocationInMM[i] - ImageCenterLocationInMM[i];
      }
    if( CurrentIndex[0] == ( CenterIndex[0] ) && CurrentIndex[1] ==
        ( CenterIndex[1] ) && CurrentIndex[2] == ( CenterIndex[2] ) )
      {
      std::cout << "CENTER_MATH AT (" << CurrentIndex << "): "
                << LocationWithRespectToCenterOfImageInMM << " = "
                << CurrentLocationInMM << " - " << ImageCenterLocationInMM
                << std::endl;
      }
    float rho, phi, theta;
    XYZToSpherical(LocationWithRespectToCenterOfImageInMM,
                   rho, phi, theta);

    rhoit.Set(rho);
    phiit.Set(phi);
    thetait.Set(theta);
    ++it;
    ++rhoit;
    ++phiit;
    ++thetait;
    }

  /* Implementation 1.
    * std::string ANNConfigurationFilenamePath =
    * itksys::SystemTools::GetFilenamePath(ANNConfigurationFilename);
    * if ( ANNConfigurationFilenamePath == "" )
    *  {
    *  ANNConfigurationFilenamePath = "." ;
    *  }
    * std::string RhoMapName(ANNConfigurationFilenamePath); RhoMapName +=
    * "/rho.nii.gz";
    * std::string PhiMapName(ANNConfigurationFilenamePath); PhiMapName +=
    * "/phi.nii.gz";
    * std::string ThetaMapName(ANNConfigurationFilenamePath); ThetaMapName +=
    * "/theta.nii.gz";
    */
  // Implementation 2.
  // REGINA:: Read Rho, phi, theta file name from xml
  std::cout << " Write rho,phi and theta" << std::endl;
  std::string RhoMapName = atlasDataSet->GetSpatialLocationFilenameByType("rho");
  std::string PhiMapName = atlasDataSet->GetSpatialLocationFilenameByType("phi");
  std::string ThetaMapName = atlasDataSet->GetSpatialLocationFilenameByType("theta");
  // Check if rho,phi and theta file exists.
  itkUtil::WriteImage<InternalImageType>(rhoImage, RhoMapName);
  itkUtil::WriteImage<InternalImageType>(phiImage, PhiMapName);
  itkUtil::WriteImage<InternalImageType>(thetaImage, ThetaMapName);
}

static int GenerateProbabilityMaps(NetConfiguration & ANNConfiguration,
                                   bool verbose)
{
  // Get the atlas dataset
  if( verbose )
    {
    std::cout << "* GenerateProbabilityMaps" << std::endl;
    }

  // find out the registration parameters
  RegistrationConfigurationParser *regParams =
    ANNConfiguration.Get<RegistrationConfigurationParser>("RegistrationConfigurationParser");
  const std::string imageTypeToUse
  (
    regParams->GetAttribute<StringValue>(
      "ImageTypeToUse") );
  const std::string regID
  (
    regParams->GetAttribute<StringValue>("ID") );

  // generate the registrations needed all in one go, using
  // multiple processors if possible
  const unsigned int numThreads = 1; //

  // itk::MultiThreader::GetGlobalDefaultNumberOfThreads();

  GenerateRegistrations(ANNConfiguration, false, false,  numThreads);

  // Get Atlas Image Name
  DataSet *   atlasDataSet = ANNConfiguration.GetAtlasDataSet();
  std::string AtlasImageFilename( atlasDataSet->GetImageFilenameByType(
                                    imageTypeToUse) );
  // HACK:  radius is hardcoded to size 2, it should be taken from command line
  // arguments.
  InternalImageType::SizeType radius;
  radius[0] = 0; radius[1] = 0; radius[2] = 0;
  InternalImageType::Pointer AtlasImage =
    ReadMedianFilteredImage<InternalImageType>(AtlasImageFilename,
                                               radius);

  GeneratePolarImages(atlasDataSet, /*probFilename,*/ AtlasImage);

  ProbabilityMapList *probMaps = ANNConfiguration.Get<ProbabilityMapList>(
      "ProbabilityMapList");
  const DataSet::StringVectorType ProbMapByMaskStringVectorType =
    probMaps->CollectAttValues<ProbabilityMapParser>("StructureID");
  const unsigned int probMapCount = ProbMapByMaskStringVectorType.size();
  assert(probMapCount > 0);

  std::vector<unsigned int>               accumulatorcounter(probMapCount);
  std::vector<InternalImageType::Pointer> accumulators(probMapCount);
  for( unsigned int accumindex = 0;
       accumindex < accumulators.size();
       accumindex++ )
    {
    CreateNewFloatImageFromTemplate(accumulators[accumindex], AtlasImage);
    accumulatorcounter[accumindex] = 0;
    }
  // For each dataset, create deformation fields if they don't already exist.
  std::list<DataSet *> dataSets = ANNConfiguration.GetTrainDataSets();
  // For each mask type, you need an accumulator image, so find all mask types,
  // and initialize the images to zeros with size and resolution of the
  // AtlasImage.
  DataSet::StringVectorType allMaskTypes( ( *( dataSets.begin() ) )->GetMaskTypes() );
  // Traverse all datasets looking for mask types.
  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end();
       ++it )
    {
    DataSet::StringVectorType currDatasetMaskTypes( ( *it )->GetMaskTypes() );
    for( DataSet::StringVectorType::iterator maskit = currDatasetMaskTypes.begin();
         maskit != currDatasetMaskTypes.end();
         ++maskit )
      {
      DataSet::StringVectorType::iterator currmaskiterator = std::find(
          allMaskTypes.begin(),
          allMaskTypes.end(),
          *maskit);
      if( currmaskiterator == allMaskTypes.end() )
        {
        allMaskTypes.push_back(*maskit);
        }
      }
    }
  /*REGIAN:: Get rid if "#if0" */
  std::cout << "===== BEGIN probMaps ====" << std::endl;
  for( unsigned int pmindex = 0; pmindex < probMaps->size(); pmindex++ )
    {
    std::cout << ProbMapByMaskStringVectorType[pmindex] << std::endl;
    DataSet::StringVectorType::iterator currmaskiterator = std::find(
        allMaskTypes.begin(),
        allMaskTypes.end(),
        ProbMapByMaskStringVectorType[pmindex]);
    if( currmaskiterator == allMaskTypes.end() )
      {
      itkGenericExceptionMacro(<< "ERROR:  No training data sets found with structure: "
                               << ProbMapByMaskStringVectorType[pmindex]);
      }
    else
      {
      std::cout << "Found Match For: " << *currmaskiterator << std::endl;
      }
    }
  std::cout << "===== END probMaps ====" << std::endl;
  // #####################Start quick processing to ensure that files exists
  // ########################
    {
    bool AllGood = true;
    // For each of the datasets do:
    for( std::list<DataSet *>::iterator it = dataSets.begin();
         it != dataSets.end(); ++it )
      {
      // get subject image
      std::string SubjectImage( ( *it )->GetImageFilenameByType(imageTypeToUse) );
      // Get SubjtoAtlasRegistrationFilename
      // RegistrationType *reg = (*it)->Get<RegistrationType>("Registration");
      const RegistrationType *reg = ( *it )->GetRegistrationWithID(regID);
      std::string             SubjToAtlasRegistrationFilename
      (
        reg->GetAttribute<StringValue>(
          "SubjToAtlasRegistrationFilename") );
      std::string landmarkType("");

      if( !itksys::SystemTools::FileExists( SubjToAtlasRegistrationFilename.
                                            c_str() ) )
        {
        std::string errmsg(SubjToAtlasRegistrationFilename);
        errmsg += " does not exist";
        throw  BRAINSCutExceptionStringHandler(errmsg);
        }
      // for each mask in this dataset, register it to atlas space
      // and then add it to the accumulator.
      for( unsigned i = 0; i < probMapCount; i++ )
        {
        const std::string MaskName( ( *it )->GetMaskFilenameByType(
                                      ProbMapByMaskStringVectorType[i]) );
        if( MaskName == "" || !itksys::SystemTools::FileExists( MaskName.c_str() ) )
          {
          continue;
          }
        if( !itksys::SystemTools::FileExists( MaskName.c_str() ) )
          {
          std::cerr << "Required mask is missing " << MaskName << std::endl;
          AllGood = false;
          }
        }
      }
    if( !AllGood )
      {
      itkGenericExceptionMacro(<< "Can not continue until there are no missing files.");
      }
    std::cout << "File existance verification Done." << std::endl;
    }
  // #####################Stop quick processing to ensure that files exists
  // ########################
  // For each of the datasets do:
  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end(); ++it )
    {
    // get subject image
    std::string SubjectImage( ( *it )->GetImageFilenameByType(imageTypeToUse) );
    // Get SubjtoAtlasRegistrationFilename
    // RegistrationType *reg = (*it)->Get<RegistrationType>("Registration");
    const RegistrationType *reg = ( *it )->GetRegistrationWithID(regID);
    std::string             SubjToAtlasRegistrationFilename
    (
      reg->GetAttribute<StringValue>(
        "SubjToAtlasRegistrationFilename") );
    /*REGINA:: Remove #if 0 */
    std::string landmarkType("");
    // for each mask in this dataset, register it to atlas space
    // and then add it to the accumulator.
    for( unsigned i = 0; i < probMapCount; i++ )
      {
      const std::string MaskName( ( *it )->GetMaskFilenameByType(
                                    ProbMapByMaskStringVectorType[i]) );
      if( MaskName == "" || !itksys::SystemTools::FileExists( MaskName.c_str() ) )
        {
        std::cout << "---Mask missing for " << SubjectImage << " "
                  << SubjToAtlasRegistrationFilename << " " << MaskName.c_str()
                  << std::endl;
        continue;
        }
      if( !itksys::SystemTools::FileExists( MaskName.c_str() ) )
        {
        itkGenericExceptionMacro(<< "Required mask is missing " << MaskName);
        }
      std::cout << "Deforming " << MaskName << " with "
                << SubjToAtlasRegistrationFilename
                << "to atlas " << AtlasImageFilename << std::endl;

      // TODO
      // - Check if this Image Warper deals properly with Binary deformation
      InternalImageType::Pointer DeformedMask =
        ImageWarper<InternalImageType>(SubjToAtlasRegistrationFilename, MaskName, AtlasImage);

      // It has to be make sure that read-in mask image has binary values with
      // one and zero. Otherwise, summing will mess up the average calculation.
      typedef itk::BinaryThresholdImageFilter<InternalImageType, InternalImageType>
        UnifyMaskLabelFilterType;
      UnifyMaskLabelFilterType::Pointer unifyMaskLabelFilter =
        UnifyMaskLabelFilterType::New();

      unifyMaskLabelFilter->SetInput( DeformedMask );
      unifyMaskLabelFilter->SetInsideValue( 1.0F );
      unifyMaskLabelFilter->SetOutsideValue( 0.0F );
      unifyMaskLabelFilter->SetLowerThreshold( 0.5F );
      unifyMaskLabelFilter->Update();

      if( verbose > 1 )
        {
        std::cout << "=======================Deformed Mask Information \n"
                  << " * Origin : " << DeformedMask->GetOrigin()
                  << std::endl;
        }
      VerifyNonZeroImage<InternalImageType>(DeformedMask, MaskName);
      try
        {
        typedef itk::AddImageFilter<InternalImageType, InternalImageType,
                                    InternalImageType> AddType;
        AddType::Pointer add = AddType::New();
        add->SetInput1( unifyMaskLabelFilter->GetOutput() );
        add->SetInput2(accumulators[i]);
        add->Update();
        accumulators[i] = add->GetOutput();
        accumulatorcounter[i] = accumulatorcounter[i] + 1;
        }
      catch( itk::ExceptionObject & e )
        {
        throw;
        }
      }
    }
  // now, for each accumulator image, divide by the #of datasets excluding
  // the atlas dataset.
  for( unsigned int i = 0; i < probMapCount; i++ )
    {
    std::cout
      << " =====================================================================";
    std::cout << " Number of Probability Map :: " << i << " / "
              << probMapCount << std::endl;
    std::cout
      << " =====================================================================";
    // find the corresponding ProbabilityMapParser
    ProbabilityMapParser *probMapObject =
      probMaps->GetMatching<ProbabilityMapParser>( "StructureID",
                                                   ProbMapByMaskStringVectorType[i].c_str() );
    if( probMapObject == 0 )
      {
      itkGenericExceptionMacro(<< "Missing Probability map for map structure: "
                               << ProbMapByMaskStringVectorType[i]);
      }
    const std::string probFilename
    (
      probMapObject->GetAttribute<StringValue>(
        "Filename") );
    VerifyNonZeroImage<InternalImageType>(accumulators[i],
                                          probFilename + "AccumulatorPreDiv");
    const double AveragerValue =
      ( 1.0 / static_cast<double>( accumulatorcounter[i] ) );
    std::cout << "#######AveragerValue " << AveragerValue << " for "
              << ProbMapByMaskStringVectorType[i] <<  " " << i << std::endl;
    accumulators[i] = ImageMultiplyConstant<InternalImageType>(accumulators[i],
                                                               AveragerValue);
    VerifyNonZeroImage<InternalImageType>(accumulators[i],
                                          probFilename + "AccumulatorPostDiv");
    const double GaussianValue =
      probMapObject->GetAttribute<FloatValue>("Gaussian");
    InternalImageType::Pointer currentImage;
    if( GaussianValue > 1e-2 )  /*Use a very small sigma value as decision
                                  * rule*/
      {
      currentImage = GaussianSmooth(accumulators[i], GaussianValue);
      }
    else
      {
      currentImage =
        itkUtil::ScaleAndCast<InternalImageType, InternalImageType>(
          accumulators[i],
          0,
          HUNDRED_PERCENT_VALUE);
      }
    itk::ImageRegionIterator<InternalImageType> it(
      currentImage, currentImage->GetLargestPossibleRegion() );
    it.GoToBegin();
    while( !it.IsAtEnd() )
      {
      if( it.Value() > ( HUNDRED_PERCENT_VALUE - PERCENT_MIN_MAX_TOLERANCE ) )
        {
        it.Value() = HUNDRED_PERCENT_VALUE;
        }
      if( it.Value() < ( 0.0F + PERCENT_MIN_MAX_TOLERANCE ) )
        {
        it.Value() = 0.0F;
        }
      ++it;
      }

    // HACK-- Should just make all the images float,
    // small roundoff errors are causing huge problems with 254 instead of
    // HUNDRED_PERCENT_VALUE for values. down stream
    // depending on gaussian values.

    VerifyNonZeroImage<InternalImageType>(currentImage, probFilename);
    try
      {
      // std::cout << currentImage << "\n  " << probFilename << std::endl;
      itkUtil::WriteImage<InternalImageType>( currentImage,
                                              probFilename.c_str() );
      }
    catch( ... )
      {
      std::cerr << "Error writing " << probFilename << std::endl;
      return -1;
      }
    } // end of probability map loop
  return 0;
}

int GenerateProbability(const std::string & XMLFile,
                        int verbose,
                        bool validate)
{
  if( verbose )
    {
    std::cerr << "This will print out process in more detail for debugging purpose\n"
              << std::endl;
    }
  NetConfiguration *     ANNConfiguration;
  NetConfigurationParser ANNConfigurationParser = NetConfigurationParser( XMLFile );

  ANNConfiguration = ANNConfigurationParser.GetNetConfiguration();

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
  if( validate )
    {
    std::cout << " *~*~*~*~*~*~*~*~*Start Validation.... *~*~*~*~*~*~*~*~*"  << std::endl;
    ANNConfigurationParser.ValidateDataSets();
    std::cout << " *~*~*~*~*~*~*~*~*End of Validation.... *~*~*~*~*~*~*~*~*"  << std::endl;
    }
  int rval = -1;
  try
    {
    rval = GenerateProbabilityMaps(*ANNConfiguration,  verbose);
    }
  catch( BRAINSCutExceptionStringHandler & ex )
    {
    std::cerr << ex.Error() << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Unidentified exception" << std::endl;
    }
  return rval;
}
