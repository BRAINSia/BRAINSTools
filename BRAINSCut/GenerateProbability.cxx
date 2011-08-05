#include <ProcessDescription.h>
#include <itksys/SystemTools.hxx>
#include <list>
#include <vector>
#include <iostream>
#include <itkIO.h>
#include "Utilities.h"
#include "Parser.h"

// ValidateDataSets function will read throught he XML files,
// and for each data set, it will ensure that all
// the image files are in the same physical space.
bool ValidateDataSets(ProcessDescription & ANNXMLObject)
{
  // HACK:  Needed to speed up testing.
  std::list<DataSet *> dataSets = ANNXMLObject.GetTrainDataSets();

  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end(); ++it )
    {
    DataSet::TypeVector allImageTypes( ( *it )->ImageTypes() );
    const std::string   FirstImageName(
      ( *it )->GetImageFilenameByType(allImageTypes[0]) );
    ProbabilityMapImageType::Pointer FirstImage =
      itkUtil::ReadImage<ProbabilityMapImageType>(FirstImageName);
    FirstImage =
      itkUtil::ScaleAndCast<ProbabilityMapImageType, ProbabilityMapImageType>(
        FirstImage,
        0,
        HUNDRED_PERCENT_VALUE);
    for( unsigned int imindex = 1; imindex < allImageTypes.size(); imindex++ )
      {
      const std::string CurrentImageName( ( *it )->GetImageFilenameByType(
                                            allImageTypes[imindex]) );
      if( !itksys::SystemTools::FileExists( CurrentImageName.c_str() ) )
        {
        std::string errmsg(CurrentImageName);
        errmsg += " does not exist";
        throw  ProcessObjectException(errmsg);
        }
      ProbabilityMapImageType::Pointer CurrentImage =
        itkUtil::ReadImage<ProbabilityMapImageType>(CurrentImageName);
      CurrentImage =
        itkUtil::ScaleAndCast<ProbabilityMapImageType,
                              ProbabilityMapImageType>(
          CurrentImage,
          0,
          HUNDRED_PERCENT_VALUE);
      if( !itkUtil::ImagePhysicalDimensionsAreIdentical<
            ProbabilityMapImageType, ProbabilityMapImageType>(FirstImage,
                                                              CurrentImage) )
        {
        std::string errmsg(CurrentImageName);
        errmsg += " and ";
        errmsg += FirstImageName;
        errmsg += " differ";
        std::cout << "============" << FirstImageName << "===============\n"
                  << FirstImage << std::endl;
        std::cout << "============" << CurrentImageName
                  << "===============\n" << CurrentImage << std::endl;
        throw  ProcessObjectException(errmsg);
        }
      } // Each Image

    DataSet::TypeVector allMaskTypes( ( *it )->MaskTypes() );
    for( unsigned maskindex = 0; maskindex < allMaskTypes.size(); maskindex++ )
      {
      const std::string CurrentMaskName( ( *it )->GetMaskFilenameByType(
                                           allMaskTypes[maskindex]) );
      if( !itksys::SystemTools::FileExists( CurrentMaskName.c_str() ) )
        {
        std::string errmsg(CurrentMaskName);
        errmsg += " does not exist";
        throw  ProcessObjectException(errmsg);
        }
      ProbabilityMapImageType::Pointer CurrentMask =
        itkUtil::ReadImage<ProbabilityMapImageType>(CurrentMaskName);
      CurrentMask =
        itkUtil::ScaleAndCast<ProbabilityMapImageType,
                              ProbabilityMapImageType>(
          CurrentMask,
          0,
          HUNDRED_PERCENT_VALUE);
      if( !itkUtil::ImagePhysicalDimensionsAreIdentical<
            ProbabilityMapImageType, ProbabilityMapImageType>(FirstImage,
                                                              CurrentMask) )
        {
        std::string errmsg(CurrentMaskName);
        errmsg += " and ";
        errmsg += FirstImageName;
        errmsg += " differ";
        std::cout << "============" << FirstImageName << "===============\n"
                  << FirstImage << std::endl;
        std::cout << "============" << CurrentMaskName
                  << "===============\n" << CurrentMask << std::endl;
        throw  ProcessObjectException(errmsg);
        }
      } // Each Mask
    }   // Each DataSet
  return true;
}

static ProbabilityMapImageType::Pointer
GaussianSmooth(const RealImageType::Pointer image,
               const double Sigma)
{
  RealImageType::SpacingType imageSpacing = image->GetSpacing();

  for( int currDim = 0; currDim < RealImageType::ImageDimension; currDim++ )
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
  RealImageType::Pointer gaussianImage = g_gaussianFilterZ->GetOutput();
  // NOTE:  Sometimes the "SetNormalizeAcrossScale" will make the zeros slightly
  // above zero.
  ProbabilityMapImageType::Pointer returnSmoothed =
    itkUtil::ScaleAndCast<RealImageType, ProbabilityMapImageType>(
      gaussianImage,
      0,
      HUNDRED_PERCENT_VALUE);
  return returnSmoothed;
}

//
// this is to remove warnings about assigning a float to a long int.
// I could have just added a cast, but it would make hairy code even uglier,
inline ProbabilityMapImageType::IndexType::IndexValueType
TruncatedHalf(const ProbabilityMapImageType::SizeType::SizeValueType & v)
{
  return static_cast<ProbabilityMapImageType::IndexType::IndexValueType>
         ( static_cast<double>( v ) * 0.5 );
}

// Generate PolarImages about the center of the Atlas Image.
static void GeneratePolarImages(ProbabilityMap *probMapObject, /*const
                                                                 * std::string
                                                                 *
                                                                 *&ANNXMLObjectFilename,*/
                                const ProbabilityMapImageType::Pointer AtlasImage)
{
  ProbabilityMapImageType::IndexType CenterIndex;

    {
    const ProbabilityMapImageType::SizeType
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

  RealImageType::Pointer rhoImage, phiImage, thetaImage;
  CreateNewFloatImageFromTemplate(rhoImage, AtlasImage);
  CreateNewFloatImageFromTemplate(phiImage, AtlasImage);
  CreateNewFloatImageFromTemplate(thetaImage, AtlasImage);

  itk::ImageRegionIterator<ProbabilityMapImageType> it(
    AtlasImage, AtlasImage->GetLargestPossibleRegion() );
  it.GoToBegin();
  itk::ImageRegionIterator<RealImageType> rhoit(
    rhoImage, rhoImage->GetLargestPossibleRegion() );
  rhoit.GoToBegin();
  itk::ImageRegionIterator<RealImageType> phiit(
    phiImage, phiImage->GetLargestPossibleRegion() );
  phiit.GoToBegin();
  itk::ImageRegionIterator<RealImageType> thetait(
    thetaImage, thetaImage->GetLargestPossibleRegion() );
  thetait.GoToBegin();
  itk::Point<float, 3> CurrentLocationInMM;
  itk::Point<float, 3> LocationWithRespectToCenterOfImageInMM;
  while( !it.IsAtEnd() )
    {
    const ProbabilityMapImageType::IndexType CurrentIndex = it.GetIndex();
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
    * std::string ANNXMLObjectFilenamePath =
    * itksys::SystemTools::GetFilenamePath(ANNXMLObjectFilename);
    * if ( ANNXMLObjectFilenamePath == "" )
    *  {
    *  ANNXMLObjectFilenamePath = "." ;
    *  }
    * std::string RhoMapName(ANNXMLObjectFilenamePath); RhoMapName +=
    * "/rho.nii.gz";
    * std::string PhiMapName(ANNXMLObjectFilenamePath); PhiMapName +=
    * "/phi.nii.gz";
    * std::string ThetaMapName(ANNXMLObjectFilenamePath); ThetaMapName +=
    * "/theta.nii.gz";
    */
  // Implementation 2.
  // REGINA:: Read Rho, phi, theta file name from xml
  std::cout << " Write rho,phi and theta" << std::endl;
  std::string RhoMapName = probMapObject->GetAttribute<StringValue>("rho");
  std::string PhiMapName = probMapObject->GetAttribute<StringValue>("phi");
  std::string ThetaMapName = probMapObject->GetAttribute<StringValue>("theta");
  // Check if rho,phi and theta file exists.
  itkUtil::WriteImage<RealImageType>(rhoImage, RhoMapName);
  itkUtil::WriteImage<RealImageType>(phiImage, PhiMapName);
  itkUtil::WriteImage<RealImageType>(thetaImage, ThetaMapName);
}

static int GenerateProbabilityMaps(ProcessDescription & ANNXMLObject,
                                   bool verbose)
{
  // Get the atlas dataset
  if( verbose )
    {
    std::cout << "* GenerateProbabilityMaps" << std::endl;
    }
  DataSet *atlasDataSet = ANNXMLObject.GetAtlasDataSet();

  // find out the registration parameters
  RegistrationParams *regParams =
    ANNXMLObject.Get<RegistrationParams>("RegistrationParams");
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

  GenerateRegistrations(ANNXMLObject, false, false, numThreads);

  // Get Atlas Image Name
  std::string AtlasImageFilename( atlasDataSet->GetImageFilenameByType(
                                    imageTypeToUse) );
  // HACK:  radius is hardcoded to size 2, it should be taken from command line
  // arguments.
  ProbabilityMapImageType::SizeType radius;
  radius[0] = 0; radius[1] = 0; radius[2] = 0;
  ProbabilityMapImageType::Pointer AtlasImage =
    ReadMedianFilteredImage<ProbabilityMapImageType>(AtlasImageFilename,
                                                     radius);

  ProbabilityMapList *probMaps = ANNXMLObject.Get<ProbabilityMapList>(
      "ProbabilityMapList");
  const DataSet::TypeVector ProbMapByMaskTypeVector =
    probMaps->CollectAttValues<ProbabilityMap>("StructureID");
  const unsigned int probMapCount = ProbMapByMaskTypeVector.size();
  assert(probMapCount > 0);

  std::vector<unsigned int>           accumulatorcounter(probMapCount);
  std::vector<RealImageType::Pointer> accumulators(probMapCount);
  for( unsigned int accumindex = 0;
       accumindex < accumulators.size();
       accumindex++ )
    {
    CreateNewFloatImageFromTemplate(accumulators[accumindex], AtlasImage);
    accumulatorcounter[accumindex] = 0;
    }
  // For each dataset, create deformation fields if they don't already exist.
  std::list<DataSet *> dataSets = ANNXMLObject.GetTrainDataSets();
  // For each mask type, you need an accumulator image, so find all mask types,
  // and initialize the images to zeros with size and resolution of the
  // AtlasImage.
  DataSet::TypeVector allMaskTypes( ( *( dataSets.begin() ) )->MaskTypes() );
  // Traverse all datasets looking for mask types.
  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end();
       ++it )
    {
    DataSet::TypeVector currDatasetMaskTypes( ( *it )->MaskTypes() );
    for( DataSet::TypeVector::iterator maskit = currDatasetMaskTypes.begin();
         maskit != currDatasetMaskTypes.end();
         ++maskit )
      {
      DataSet::TypeVector::iterator currmaskiterator = std::find(
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
    std::cout << ProbMapByMaskTypeVector[pmindex] << std::endl;
    DataSet::TypeVector::iterator currmaskiterator = std::find(
        allMaskTypes.begin(),
        allMaskTypes.end(),
        ProbMapByMaskTypeVector[pmindex]);
    if( currmaskiterator == allMaskTypes.end() )
      {
      std::cout << "ERROR:  No training data sets found with structure: "
                << ProbMapByMaskTypeVector[pmindex] <<  std::endl;
      exit(-1);
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
        throw  ProcessObjectException(errmsg);
        }
      // for each mask in this dataset, register it to atlas space
      // and then add it to the accumulator.
      for( unsigned i = 0; i < probMapCount; i++ )
        {
        const std::string MaskName( ( *it )->GetMaskFilenameByType(
                                      ProbMapByMaskTypeVector[i]) );
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
      std::cerr << "Can not continue until there are no missing files."
                << std::endl;
      exit(-1);
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
                                    ProbMapByMaskTypeVector[i]) );
      if( MaskName == "" || !itksys::SystemTools::FileExists( MaskName.c_str() ) )
        {
        std::cout << "---Mask missing for " << SubjectImage << " "
                  << SubjToAtlasRegistrationFilename << " " << MaskName.c_str()
                  << std::endl;
        continue;
        }
      if( !itksys::SystemTools::FileExists( MaskName.c_str() ) )
        {
        std::cerr << "Required mask is missing " << MaskName << std::endl;
        exit(-1);
        }
      std::cout << "Deforming " << MaskName << " with "
                << SubjToAtlasRegistrationFilename
                << "to atlas " << AtlasImageFilename << std::endl;

      // TODO
      // - Check if this Image Warper deals properly with Binary deformation
      RealImageType::Pointer DeformedMask =
        ImageWarper<RealImageType>(SubjToAtlasRegistrationFilename, MaskName, AtlasImage);

      // It has to be make sure that read-in mask image has binary values with
      // one and zero. Otherwise, summing will mess up the average calculation.
      typedef itk::BinaryThresholdImageFilter<RealImageType, RealImageType>
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
      VerifyNonZeroImage<RealImageType>(DeformedMask, MaskName);
      try
        {
        typedef itk::AddImageFilter<RealImageType, RealImageType,
                                    RealImageType> AddType;
        AddType::Pointer add = AddType::New();
        add->SetInput1( unifyMaskLabelFilter->GetOutput() );
        add->SetInput2(accumulators[i]);
        add->Update();
        accumulators[i] = add->GetOutput();
        accumulatorcounter[i] = accumulatorcounter[i] + 1;
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception in Image Accumulation." << std::endl;
        std::cerr << e.GetDescription() << std::endl;
        std::cerr << e.GetLocation() << std::endl;
        exit(-1);
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
    // find the corresponding ProbabilityMap
    ProbabilityMap *probMapObject =
      probMaps->GetMatching<ProbabilityMap>( "StructureID",
                                             ProbMapByMaskTypeVector[i].c_str() );
    if( probMapObject == 0 )
      {
      std::cout << "Missing Probability map for map structure: "
                << ProbMapByMaskTypeVector[i] << std::endl;
      exit(-1);
      }
    const std::string probFilename
    (
      probMapObject->GetAttribute<StringValue>(
        "Filename") );
    VerifyNonZeroImage<RealImageType>(accumulators[i],
                                      probFilename + "AccumulatorPreDiv");
    const double AveragerValue =
      ( 1.0 / static_cast<double>( accumulatorcounter[i] ) );
    std::cout << "#######AveragerValue " << AveragerValue << " for "
              << ProbMapByMaskTypeVector[i] <<  " " << i << std::endl;
    accumulators[i] = ImageMultiplyConstant<RealImageType>(accumulators[i],
                                                           AveragerValue);
    VerifyNonZeroImage<RealImageType>(accumulators[i],
                                      probFilename + "AccumulatorPostDiv");
    const double GaussianValue =
      probMapObject->GetAttribute<FloatValue>("Gaussian");
    ProbabilityMapImageType::Pointer currentImage;
    if( GaussianValue > 1e-2 )  /*Use a very small sigma value as decision
                                  * rule*/
      {
      currentImage = GaussianSmooth(accumulators[i], GaussianValue);
      }
    else
      {
      currentImage =
        itkUtil::ScaleAndCast<RealImageType, ProbabilityMapImageType>(
          accumulators[i],
          0,
          HUNDRED_PERCENT_VALUE);
      }
    itk::ImageRegionIterator<ProbabilityMapImageType> it(
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

    VerifyNonZeroImage<ProbabilityMapImageType>(currentImage, probFilename);
    try
      {
      // std::cout << currentImage << "\n  " << probFilename << std::endl;
      itkUtil::WriteImage<ProbabilityMapImageType>( currentImage,
                                                    probFilename.c_str() );
      }
    catch( ... )
      {
      std::cerr << "Error writing " << probFilename << std::endl;
      return -1;
      }
    // generate Rho/Phi/Theta images
    //    if ( i == 0 )                  // we need every polor image
    //      {
    GeneratePolarImages(probMapObject, /*probFilename,*/ AtlasImage);
    //      }
    } // end of probability map loop
  return 0;
}

int GenerateProbability(const std::string & XMLFile, int verbose, bool validate)
{
  if( verbose )
    {
    std::cerr << "This will print out process in more detail for debugging purpose\n"
              << std::endl;
    }
  ProcessDescription ANNXMLObject;
  try
    {
    ReadXML(XMLFile.c_str(), ANNXMLObject);
    }
  catch( ProcessObjectException & ex )
    {
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
  if( validate )
    {
    std::cout << " *~*~*~*~*~*~*~*~*Start Validation.... *~*~*~*~*~*~*~*~*"  << std::endl;
    ValidateDataSets(ANNXMLObject);
    std::cout << " *~*~*~*~*~*~*~*~*End of Validation.... *~*~*~*~*~*~*~*~*"  << std::endl;
    }
  int rval = -1;
  try
    {
    rval = GenerateProbabilityMaps(ANNXMLObject, verbose);
    }
  catch( ProcessObjectException & ex )
    {
    std::cerr << ex.Error() << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Unidentified exception" << std::endl;
    }
  return rval;
}
