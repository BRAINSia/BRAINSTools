#include "BRAINSCutVectorTrainingSet.h"
#include "BRAINSCutExceptionStringHandler.h"
#include "ShuffleVectors.h"

#include <fstream>
// #include <vnl/vnl_random.h>

// #include <fcntl.h>

// ---------------------------//
BRAINSCutVectorTrainingSet
::BRAINSCutVectorTrainingSet( const std::string vectorFilename)
  : trainingVectorFilename( vectorFilename),
  trainingHeaderFilename( vectorFilename + ".hdr" ),
  totalVectorSize(0),
  inputVectorSize(0),
  outputVectorSize(0),
  shuffled(false),
  recordSize(0),
  bufferRecordSize(0),
  numberOfSubSet(1),
  currentTrainingSubSet(NULL)

{
  // trainingVectorFilename = vectorFilename;
  // trainingHeaderFilename += ".hdr";
}

BRAINSCutVectorTrainingSet
::~BRAINSCutVectorTrainingSet()
{
  cvReleaseMat( &(currentTrainingSubSet->pairedInput) );
  cvReleaseMat( &(currentTrainingSubSet->pairedOutput) );
  cvReleaseMat( &(currentTrainingSubSet->pairedOutputRF) );
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::ReadHeaderFileInformation()
{
  std::ifstream headerFileStream;

  headerFileStream.open( trainingHeaderFilename.c_str(),
                         std::ios::in );

  if( !headerFileStream.is_open() )
    {
    std::string error = "Cannot Open the file of ";
    error += trainingHeaderFilename;
    error += " ";
    error += __FILE__;
    throw BRAINSCutExceptionStringHandler( error );
    }
  for( int tags = 0; tags < 4; tags++ )  // Read header file for input/output
                                         // vectors and total vectors.
    {
    std::string temp;
    char        currentline[MAXIMUMCHAR];
    headerFileStream.getline(currentline, MAXIMUMCHAR - 1);
    std::istringstream iss(currentline, std::istringstream::in);
    iss >> temp;
    if( temp == "IVS" )
      {
      iss >> inputVectorSize;
      }
    else if( temp == "OVS" )
      {
      iss >> outputVectorSize;
      }
    else if( temp == "TVC" )
      {
      iss >> totalVectorSize;
      }
    else if( temp == "SHUFFLED" )
      {
      iss >> shuffled;
      }
    }
  headerFileStream.close();
}

// ---------------------------//
int
BRAINSCutVectorTrainingSet
::GetInputVectorSize()
{
  return inputVectorSize;
}

// ---------------------------//
int
BRAINSCutVectorTrainingSet
::GetOutputVectorSize()
{
  return outputVectorSize;
}

// ---------------------------//
int
BRAINSCutVectorTrainingSet
::GetTotalVectorSize()
{
  return totalVectorSize;
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::SetBufferRecordSize()
{
  bufferRecordSize = ( inputVectorSize + outputVectorSize + LineGuardSize );
}

void
BRAINSCutVectorTrainingSet
::SetShuffled(bool shuffledTrue)
{
  shuffled = shuffledTrue;
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::SetRecordSize()
{
  recordSize = ( inputVectorSize + outputVectorSize + LineGuardSize ) * sizeof( scalarType );
}

// ---------------------------//
unsigned int
BRAINSCutVectorTrainingSet
::GetRecordSize()
{
  return recordSize;
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::PrintDebuggingMessage(std::string msg )
{
  std::cout << " ***** note                                                   **** " << std::endl
            << " *" << msg << std::endl
            << " ******************************************************************" << std::endl;

  return;
}

// ---------------------------//
inline
int
GetFileStreamToWrite( std::string filename)
{
  int fileStreamToWrite = open( filename.c_str(),
                                O_WRONLY | O_CREAT | O_TRUNC,
                                S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH );

  if( fileStreamToWrite == -1 )
    {
    std::string msg( "Cannot Open FileStream of " );
    msg += filename;
    throw BRAINSCutExceptionStringHandler( msg );
    }
  return fileStreamToWrite;
}

// ---------------------------//
inline
void
GetFileStreamToRead( std::string filename, std::ifstream& fileStreamToRead)
{
  if( !itksys::SystemTools::FileExists( filename.c_str() ) )
    {
    std::string msg( "Vector File has not been created. " + filename );
    throw BRAINSCutExceptionStringHandler( msg);
    }
  try
    {
    fileStreamToRead.open( filename.c_str(),
                           std::ios::in | std::ios::binary );
    }
  catch( std::ifstream::failure e )
    {
    std::cout << "Exception opening file::"
              << filename << std::endl
              << e.what() << std::endl;
    }
  catch( BRAINSCutExceptionStringHandler& e )
    {
    std::cout << e.Error();
    exit(EXIT_FAILURE);
    }
  if( !fileStreamToRead.is_open() )
    {
    std::string msg( "Cannot Open FileStream of " );
    msg += filename;
    throw BRAINSCutExceptionStringHandler( msg );
    }
  return;
}

// ---------------------------//
scalarType *
BRAINSCutVectorTrainingSet
::ReadBufferFromFileStream( std::ifstream& fileStream )
{
  scalarType * buffer = new scalarType[bufferRecordSize];

  fileStream.read( (char *)buffer, recordSize );
  /*
  for( unsigned int i = 0; i < bufferRecordSize; i++ )
    {
    std::cout << buffer[i] << " ";
    }
  std::cout << std::endl;
  */

  if( buffer[bufferRecordSize - 1] != LineGuard )
    {
    throw ( BRAINSCutExceptionStringHandler( "Record not properly terminated by sentinel value") );
    }
  return buffer;
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::RandomizeTrainingVector()
{
  PrintDebuggingMessage( "this shuffling process will override current version of vector *" );

  std::string temporaryResultFilename = trainingVectorFilename;
  temporaryResultFilename += "Shuffled";

  const int        samplingProportion = 1; // shuffle order only (without up/down sampling)
  ShuffleVectors * my_ShuffleVector = new ShuffleVectors(  trainingVectorFilename,
                                                           temporaryResultFilename,
                                                           samplingProportion );

  my_ShuffleVector->ReadHeader();
  my_ShuffleVector->Shuffling();
  my_ShuffleVector->WriteHeader();

  if( rename( temporaryResultFilename.c_str(), trainingVectorFilename.c_str() ) != 0 )
    {
    std::string msg = "Fail to rename from " + temporaryResultFilename + " to " + trainingVectorFilename;
    throw ( BRAINSCutExceptionStringHandler( msg) );
    }
  else
    {
    std::string msg = "The " + temporaryResultFilename + " successfully renamed to " + trainingVectorFilename;
    PrintDebuggingMessage( msg );
    }
}

// ---------------------------//
unsigned int
BRAINSCutVectorTrainingSet
::GetNumberOfSubSet()
{
  return numberOfSubSet;
}

void
BRAINSCutVectorTrainingSet
::SetNumberOfSubSet( const unsigned int count )
{
  numberOfSubSet = count;
}

// ---------------------------//
pairedTrainingSetType *
BRAINSCutVectorTrainingSet
::GetTrainingSubSet( unsigned int count )
{
  if( !(count < numberOfSubSet) )
    {
    throw (  BRAINSCutExceptionStringHandler( "Specified SubSet is not valid") );
    }
  else if( count == currentSubSetID && currentTrainingSubSet != NULL )
    {
    return currentTrainingSubSet;
    }
  else
    {
    std::cout << " Read in New Sub Set" << std::endl;
    SetTrainingSubSet( count );
    return currentTrainingSubSet;
    }
}

// ---------------------------//
void
BRAINSCutVectorTrainingSet
::SetTrainingSubSet( unsigned int count )
{
  currentSubSetID = count;
  unsigned int subSetSize = totalVectorSize / numberOfSubSet;
  std::cout << totalVectorSize << "/" << numberOfSubSet << " = " << subSetSize << std::endl;

  std::ifstream readInFile;
  GetFileStreamToRead( trainingVectorFilename, readInFile );

  /* set to the location of this subset */
  std::ios::off_type seekval = static_cast<std::ios::off_type>( recordSize * subSetSize * count );
  readInFile.seekg( seekval, std::ios::beg);

  scalarType * currentBuffer = new scalarType[bufferRecordSize];

  scalarType * pairedInputBuffer = new scalarType[subSetSize * inputVectorSize];
  scalarType * pairedOutputBuffer = new scalarType[subSetSize * outputVectorSize];
  scalarType * pairedOutputBufferRF = new scalarType[subSetSize];  // RandomForest
  for( unsigned int i = 0; i < subSetSize  && !readInFile.eof(); i++ )
    {
    currentBuffer = ReadBufferFromFileStream( readInFile );
    /* move this to one line buffer for open cv matrix type */
    scalarType tempOutput = 0;
    for( int j = 0; j < outputVectorSize; j++ )
      {
      pairedOutputBuffer[i * outputVectorSize + j] = currentBuffer[j];
      if( currentBuffer[j] > 0.5F && tempOutput == 0 )
        {
        tempOutput = j + 1;
        }
      else if(  currentBuffer[j] > 0.5F && tempOutput != 0 )
        {
        std::cout << "A voxel belongs to more than a structure" << std::endl;
        exit(EXIT_FAILURE);
        }
      }
    pairedOutputBufferRF[i] = tempOutput;
    for( int j = 0; j < inputVectorSize; j++ )
      {
      pairedInputBuffer[i * inputVectorSize + j] = currentBuffer[j + outputVectorSize];
      }
    }
  delete[] currentBuffer;

  currentTrainingSubSet = new pairedTrainingSetType;
  currentTrainingSubSet->pairedInput = cvCreateMat( subSetSize, inputVectorSize, CV_32FC1 );
  cvInitMatHeader( currentTrainingSubSet->pairedInput,
                   subSetSize,
                   inputVectorSize,
                   CV_32FC1,
                   pairedInputBuffer);
  currentTrainingSubSet->pairedOutput = cvCreateMat( subSetSize, outputVectorSize, CV_32FC1);
  cvInitMatHeader( currentTrainingSubSet->pairedOutput,
                   subSetSize,
                   outputVectorSize,
                   CV_32FC1,
                   pairedOutputBuffer);

  currentTrainingSubSet->pairedOutputRF = cvCreateMat( subSetSize, 1, CV_32FC1);
  cvInitMatHeader( currentTrainingSubSet->pairedOutputRF,
                   subSetSize,
                   1,
                   CV_32FC1,
                   pairedOutputBufferRF);

  if( currentTrainingSubSet->pairedInput->rows == 0 )
    {
    throw (  BRAINSCutExceptionStringHandler( "Vector File is Empty!! ") );
    }
}
