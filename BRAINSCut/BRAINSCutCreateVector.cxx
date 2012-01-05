#include "BRAINSCutCreateVector.h"

BRAINSCutCreateVector
::BRAINSCutCreateVector(std::string netConfigurationFilename)
  : BRAINSCutPrimary( netConfigurationFilename)
{
  // GenerateRegistrations(BRAINSCutNetConfiguration, true, false, 1);

  SetRegistrationParametersFromNetConfiguration();

  SetRegionsOfInterestFromNetConfiguration();

  SetAtlasDataSet();
  SetRhoPhiThetaFromNetConfiguration();

  SetANNModelConfiguration();
  SetGradientSizeFromNetConfiguration();

  std::cout << "Get Normalization From NetConfiguration ";
  normalization = GetNormalizationFromNetConfiguration();
  std::cout << "(" << normalization << ")" << std::endl;
}

void
BRAINSCutCreateVector
::CreateVectors()
{
  typedef NetConfiguration::TrainDataSetListType::iterator
    TrainSubjectIteratorType;

  int numberOfInputVector = 0;

  /* open up the output stream */
  SetTrainingVectorFilenameFromNetConfiguration();
  std::ofstream outputVectorStream;

  std::string vectorFileDirectory = itksys::SystemTools::GetFilenamePath( vectorFilename.c_str() );

  if( !itksys::SystemTools::FileExists( vectorFileDirectory.c_str(), false ) )
    {
    std::cout << " OutputVectorDirectory does not exist. Create as following:: "
              << vectorFileDirectory.c_str()
              << std::endl;
    itksys::SystemTools::MakeDirectory( vectorFileDirectory.c_str() );
    }

  outputVectorStream.open( vectorFilename.c_str(),
                           std::ios::out | std::ios::binary | std::ios::trunc);

  if( !outputVectorStream.good() )
    {
    itkGenericExceptionMacro(<< "Error: Could not open ANN vector file: "
                             << vectorFilename)
    }
  for( TrainSubjectIteratorType subjectIt = trainDataSetList.begin();
       subjectIt != trainDataSetList.end();
       ++subjectIt )
    {
    numberOfInputVector += CreateSubjectVectors( *(*subjectIt), outputVectorStream);
    std::cout << "Number Of InputVector : " << numberOfInputVector << std::endl;
    }
  outputVectorStream.close();

  WriteHeaderFile( this->inputVectorSize, outputVectorSize, numberOfInputVector);
}

void
BRAINSCutCreateVector
::SetTrainingDataSetFromNetConfiguration()
{
  try
    {
    trainDataSetList = BRAINSCutNetConfiguration.GetTrainDataSets();
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    throw e;
    }
}

int
BRAINSCutCreateVector
::CreateSubjectVectors( DataSet& subject, std::ofstream& outputStream )
{
  std::map<std::string, WorkingImagePointer> deformedSpatialLocationImageList;

  GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  GetImagesOfSubjectInOrder(imagesOfInterest, subject);

  std::map<std::string, WorkingImagePointer> deformedROIs;
  GetDeformedROIs(deformedROIs, subject);

  FeatureInputVector inputVectorGenerator;

  inputVectorGenerator.SetGradientSize( gradientSize );
  inputVectorGenerator.SetImagesOfInterestInOrder( imagesOfInterest );
  inputVectorGenerator.SetImagesOfSpatialLocation( deformedSpatialLocationImageList );
  inputVectorGenerator.SetCandidateROIs( deformedROIs);
  inputVectorGenerator.SetROIInOrder( roiIDsInOrder);
  inputVectorGenerator.SetInputVectorSize();
  inputVectorGenerator.SetNormalization( normalization );

  inputVectorSize = inputVectorGenerator.GetInputVectorSize(); // TODO
  outputVectorSize = roiIDsInOrder.size();

  /* now iterate through the roi */
  unsigned int roiIDsOrderNumber = 0;

  int numberOfVectors = 0;
  for( DataSet::StringVectorType::iterator roiTyIt = roiIDsInOrder.begin();
       roiTyIt != roiIDsInOrder.end();
       ++roiTyIt ) // roiTyIt = Region of Interest Type Iterator
    {
    ProbabilityMapParser* roiDataSet =
      roiDataList->GetMatching<ProbabilityMapParser>( "StructureID", (*roiTyIt).c_str() );

    if( roiDataSet->GetAttribute<StringValue>("GenerateVector") == "true" )
      {
      /* get input vector */
      InputVectorMapType roiInputVector = inputVectorGenerator.GetFeatureInputOfROI( *roiTyIt );

      /*
       * get paired output vector
       * * subjectROIBinaryName = given answer = ground truth of segmentation = manual
       */
      std::string subjectROIBinaryName = GetROIBinaryFilename( subject, *roiTyIt );
      std::cout << __LINE__ << "::" << __FILE__ << "::" << subjectROIBinaryName << std::endl;

      OutputVectorMapType roiOutputVector = GetPairedOutput( deformedROIs, *roiTyIt,
                                                             subjectROIBinaryName, roiIDsOrderNumber );
      WriteCurrentVectors( roiInputVector, roiOutputVector, outputStream );
      numberOfVectors += roiInputVector.size();
      }

    roiIDsOrderNumber++;
    }

  return numberOfVectors;
}

inline std::string
BRAINSCutCreateVector
::GetROIBinaryFilename( DataSet& subject, std::string roiName)
{
  return std::string(subject.GetMaskFilenameByType( roiName ) );
}

OutputVectorMapType
BRAINSCutCreateVector
::GetPairedOutput( std::map<std::string, WorkingImagePointer>& deformedROIs,
                   std::string roiName,
                   std::string subjectROIBinaryFilename,
                   int roiNumber)
{
  std::cout << __LINE__ << "::" << __FILE__ << "::" << roiName << std::endl;
  std::cout << __LINE__ << "::" << __FILE__ << "::" << subjectROIBinaryFilename << std::endl;
  WorkingImagePointer subjectROIBinaryImage = ReadImageByFilename( subjectROIBinaryFilename );

  itk::ImageRegionIterator<WorkingImageType> it( deformedROIs.find( roiName )->second,
                                                 deformedROIs.find( roiName )->second->GetLargestPossibleRegion() );

  it.GoToBegin();

  OutputVectorMapType currentOutputVectorMap;
  while( !it.IsAtEnd() )
    {
    if( it.Value() > 0.0F )
      {
      const WorkingImageType::IndexType itIndex        = it.GetIndex();
      const int                         itKeyFromIndex = FeatureInputVector::HashKeyFromIndex( itIndex );

      /* fill the output vector with zeros */
      OutputVectorType oneRowOutputVector( deformedROIs.size(), 0.0F ); // working pixel type =float

      /* get the answer */
      oneRowOutputVector[roiNumber] = GetBinaryValue( subjectROIBinaryImage->GetPixel(itIndex) );

      // currentOutputVectorMap.insert( std::pair<CompatableWorkingIndexType,OutputVectorType>( itIndex,
      // oneRowOutputVector));

      currentOutputVectorMap[itKeyFromIndex] = oneRowOutputVector;
      OutputVectorMapType::const_iterator myIterator = currentOutputVectorMap.find(  itKeyFromIndex );
      }
    ++it;
    }

  return OutputVectorMapType( currentOutputVectorMap );
}

void
BRAINSCutCreateVector
::SetTrainingVectorFilenameFromNetConfiguration()
{
  vectorFilename = annModelConfiguration->GetAttribute<StringValue>("TrainingVectorFilename");
  vectorFilename += "ANN";
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "Write vector file at " << vectorFilename << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}

void
BRAINSCutCreateVector
::WriteCurrentVectors( InputVectorMapType& pairedInput,
                       OutputVectorMapType& pairedOutput,
                       std::ofstream& outputStream )
{
  int bufferSize       = (inputVectorSize + outputVectorSize + 1);

  for( InputVectorMapType::iterator it = pairedInput.begin();
       it != pairedInput.end();
       it++ )
    {
    scalarType *bufferToWrite = new scalarType[bufferSize];
    bufferToWrite[bufferSize - 1] = LineGuard;
    for( int i = 0; i < outputVectorSize; i++ )
      {
      if( pairedOutput.find( it->first) == pairedOutput.end() )
        {
        std::cout << "No output compute for this "
                  << it->first << " at " << FeatureInputVector::HashIndexFromKey( it->first )
                  << std::endl;
        }
      InputVectorType inputVector = (pairedOutput.find( it->first) )->second;

      bufferToWrite[i] = inputVector[i];
      }
    for( int j = 0; j < inputVectorSize; j++ )
      {
      bufferToWrite[outputVectorSize + j] = (it->second)[j];
      }

    outputStream.write( (const char *) bufferToWrite, bufferSize * sizeof(scalarType) );
    delete[] bufferToWrite;
    }
}

void
BRAINSCutCreateVector
::WriteHeaderFile( int LocalinputVectorSize, int LocaloutputVectorSize, int numberOfInputVector)
{
  const std::string headerFilename = vectorFilename + ".hdr";
  std::ofstream     outputStream;

  outputStream.open( headerFilename.c_str(),
                     std::ios::out | std::ios::binary );

  outputStream <<  "IVS " << LocalinputVectorSize     << std::endl
               <<  "OVS " << LocaloutputVectorSize    << std::endl
               <<  "TVC " << numberOfInputVector << std::endl;

  outputStream.close();
}

inline scalarType
BRAINSCutCreateVector
::GetBinaryValue( scalarType value)
{
  if( value > 0.5F )
    {
    return 1.0F;
    }
  else
    {
    return 0.0F;
    }
}
