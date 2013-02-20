#include "BRAINSCutCreateVector.h"

BRAINSCutCreateVector
::BRAINSCutCreateVector( BRAINSCutDataHandler dataHandler ) :
  m_inputVectorSize(0),
  m_outputVectorSize(0),
  m_myDataHandler(dataHandler)
{
  m_myDataHandler.SetRegistrationParameters();
  m_myDataHandler.SetAtlasDataSet();
  m_myDataHandler.SetRegionsOfInterest();
  m_myDataHandler.SetRhoPhiTheta();
  m_myDataHandler.SetTrainingVectorConfiguration();
  m_myDataHandler.SetGradientSize();
  m_myDataHandler.SetNormalization();
}

void
BRAINSCutCreateVector
::CreateVectors()
{
  typedef BRAINSCutConfiguration::TrainDataSetListType::iterator
    TrainSubjectIteratorType;

  int numberOfInputVector = 0;

  /* open up the output stream */
  m_myDataHandler.SetTrainVectorFilename();

  std::ofstream outputVectorStream;

  const std::string vectorFilename = m_myDataHandler.GetTrainVectorFilename();
  const std::string vectorFileDirectory
    = itksys::SystemTools::GetFilenamePath( vectorFilename.c_str() );

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
  for( TrainSubjectIteratorType subjectIt = m_trainDataSetList.begin();
       subjectIt != m_trainDataSetList.end();
       ++subjectIt )
    {
    numberOfInputVector += CreateSubjectVectors( *(*subjectIt), outputVectorStream);
    std::cout << "Number Of InputVector : " << numberOfInputVector << std::endl;
    }
  outputVectorStream.close();

  WriteHeaderFile( vectorFilename, this->m_inputVectorSize, m_outputVectorSize, numberOfInputVector);
}

void
BRAINSCutCreateVector
::SetTrainingDataSet()
{
  try
    {
    m_trainDataSetList = m_myDataHandler.GetTrainDataSet();
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

  m_myDataHandler.GetDeformedSpatialLocationImages( deformedSpatialLocationImageList, subject );

  WorkingImageVectorType imagesOfInterest;
  m_myDataHandler.ReadImagesOfSubjectInOrder(imagesOfInterest, subject);

  std::map<std::string, WorkingImagePointer> deformedROIs;
  m_myDataHandler.GetDeformedROIs(deformedROIs, subject);

  FeatureInputVector inputVectorGenerator;

  inputVectorGenerator.SetGradientSize( m_myDataHandler.GetGradientSize() );
  inputVectorGenerator.SetImagesOfInterestInOrder( imagesOfInterest );
  inputVectorGenerator.SetImagesOfSpatialLocation( deformedSpatialLocationImageList );
  inputVectorGenerator.SetCandidateROIs( deformedROIs);
  inputVectorGenerator.SetROIInOrder( m_myDataHandler.GetROIIDsInOrder() );
  inputVectorGenerator.SetInputVectorSize();
  inputVectorGenerator.SetNormalizationMethod( m_myDataHandler.GetNormalizationMethod() );

  m_inputVectorSize = inputVectorGenerator.GetInputVectorSize(); // TODO
  m_outputVectorSize = m_myDataHandler.GetROIIDsInOrder().size();

  /* now iterate through the roi */
  unsigned int roiIDsOrderNumber = 0;

  int numberOfVectors = 0;
  // for( DataSet::StringVectorType::iterator roiTyIt = m_myDataHandler.GetROIIDsInOrder().begin();
  //     roiTyIt != m_myDataHandler.GetROIIDsInOrder().end();
  //     ++roiTyIt ) // roiTyIt = Region of Interest Type Iterator
  while( roiIDsOrderNumber <  m_myDataHandler.GetROIIDsInOrder().size() )
    {
    std::string currentROI( m_myDataHandler.GetROIIDsInOrder()[roiIDsOrderNumber] );

    ProbabilityMapParser* roiDataSet =
      m_myDataHandler.GetROIDataList()->GetMatching<ProbabilityMapParser>( "StructureID", currentROI.c_str() );

    if( roiDataSet->GetAttribute<StringValue>("GenerateVector") == "true" )
      {
      /* get input vector */
      InputVectorMapType roiInputVector = inputVectorGenerator.ComputeAndGetFeatureInputOfROI( currentROI );

      /*
       * get paired output vector
       * * subjectROIBinaryName = given answer = ground truth of segmentation = manual
       */
      std::string subjectROIBinaryName = GetROIBinaryFilename( subject, currentROI );

      OutputVectorMapType roiOutputVector = GetPairedOutput( deformedROIs, currentROI,
                                                             subjectROIBinaryName, roiIDsOrderNumber );
      WriteCurrentVectors( roiInputVector, roiOutputVector, outputStream );
      numberOfVectors += roiInputVector.size();
      }

    ++roiIDsOrderNumber;
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
  if( !itksys::SystemTools::FileExists( subjectROIBinaryFilename.c_str(), false ) )
    {
    std::cout << " Subject binary file of "
              << roiName
              << " named as " << subjectROIBinaryFilename
              << " does not exist. "
              << std::endl;
    exit(EXIT_FAILURE);
    }
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
      OutputVectorType oneRowOutputVector( deformedROIs.size(), ZeroPercentValue ); // working pixel type =float

      /* get the answer */
      oneRowOutputVector[roiNumber] = GetBinaryValue( subjectROIBinaryImage->GetPixel(itIndex) );

      // currentOutputVectorMap.insert( std::pair<CompatableWorkingIndexType,OutputVectorType>( itIndex,
      // oneRowOutputVector));

      currentOutputVectorMap[itKeyFromIndex] = oneRowOutputVector;
      }
    ++it;
    }

  return OutputVectorMapType( currentOutputVectorMap );
}

void
BRAINSCutCreateVector
::WriteCurrentVectors( InputVectorMapType& pairedInput,
                       OutputVectorMapType& pairedOutput,
                       std::ofstream& outputStream )
{
  int bufferSize       = (m_inputVectorSize + m_outputVectorSize + 1);

  for( InputVectorMapType::iterator it = pairedInput.begin();
       it != pairedInput.end();
       ++it )
    {
    scalarType *bufferToWrite = new scalarType[bufferSize];
    bufferToWrite[bufferSize - 1] = LineGuard;
    for( int i = 0; i < m_outputVectorSize; ++i )
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
    for( int j = 0; j < m_inputVectorSize; ++j )
      {
      bufferToWrite[m_outputVectorSize + j] = (it->second)[j];
      }

    outputStream.write( (const char *) bufferToWrite, bufferSize * sizeof(scalarType) );
    delete[] bufferToWrite;
    }
}

void
BRAINSCutCreateVector
::WriteHeaderFile( std::string vectorFilename,
                   int LocalinputVectorSize, int LocaloutputVectorSize, int numberOfInputVector)
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
    return HundredPercentValue;
    }
  else
    {
    return ZeroPercentValue;
    }
}
