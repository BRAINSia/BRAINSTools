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
  currentTrainingSubSet(nullptr)

{
  // trainingVectorFilename = vectorFilename;
  // trainingHeaderFilename += ".hdr";
}

BRAINSCutVectorTrainingSet
::~BRAINSCutVectorTrainingSet()
{
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
  catch( std::ifstream::failure & e )
    {
    std::cout << "Exception opening file::"
              << filename << std::endl
              << e.what() << std::endl;
    }
  catch( BRAINSCutExceptionStringHandler & e )
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
    delete [] buffer;
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

  constexpr int samplingProportion = 1; // shuffle order only (without up/down sampling)
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
  std::cout<<currentTrainingSubSet<<std::endl;
  if( !(count < numberOfSubSet) )
    {
    throw (  BRAINSCutExceptionStringHandler( "Specified SubSet is not valid") );
    }
  else if( count == currentSubSetID && currentTrainingSubSet != nullptr )
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
    std::cout<<std::endl;
    }
  delete[] currentBuffer;

  currentTrainingSubSet = new pairedTrainingSetType;
  currentTrainingSubSet->pairedInput = cv::Mat( subSetSize,
                   inputVectorSize,
                   CV_32F,
                   pairedInputBuffer);
  currentTrainingSubSet->pairedOutput = cv::Mat( subSetSize,
                   outputVectorSize,
                   CV_32F,
                   pairedOutputBuffer);

  currentTrainingSubSet->pairedOutputRF = cv::Mat( subSetSize,
                   1,
                   CV_32F,
                   pairedOutputBufferRF);

  std::cout<< " n row input= " << (currentTrainingSubSet->pairedInput).rows
           << " n row output = " << (currentTrainingSubSet->pairedOutput).rows
           << std::endl;
  if( (currentTrainingSubSet->pairedInput).rows != (currentTrainingSubSet->pairedOutput).rows)
    {
    throw (  BRAINSCutExceptionStringHandler( "Input and Output row do not match!! ") );
    }
  if( (currentTrainingSubSet->pairedInput).rows == 0)
    {
    throw (  BRAINSCutExceptionStringHandler( "Vector File is Empty!! ") );
    }
}
