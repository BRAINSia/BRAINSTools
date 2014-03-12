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
/** Authors: Marius Staring and Stefan Klein **/
#include "castconverthelpers.h"

int DicomFileConverterScalarA( const std::string & inputPixelComponentType,
                               const std::string & outputPixelComponentType, const std::string & inputDirectoryName,
                               const std::string & outputFileName, int inputDimension )
{
  /** Support for 3D images. */
  if( inputDimension == 3 )
    {
    /** From int to something else. */
    callCorrectReadDicomWriterMacro( int, unsigned char );
    callCorrectReadDicomWriterMacro( int, char );
    callCorrectReadDicomWriterMacro( int, unsigned short );
    callCorrectReadDicomWriterMacro( int, short );
    callCorrectReadDicomWriterMacro( int, unsigned int );
    callCorrectReadDicomWriterMacro( int, int );
    callCorrectReadDicomWriterMacro( int, unsigned long );
    callCorrectReadDicomWriterMacro( int, long );
    callCorrectReadDicomWriterMacro( int, float );
    callCorrectReadDicomWriterMacro( int, double );

    /** From float to something else. */
    callCorrectReadDicomWriterMacro( float, unsigned char );
    callCorrectReadDicomWriterMacro( float, char );
    callCorrectReadDicomWriterMacro( float, unsigned short );
    callCorrectReadDicomWriterMacro( float, short );
    callCorrectReadDicomWriterMacro( float, unsigned int );
    callCorrectReadDicomWriterMacro( float, int );
    callCorrectReadDicomWriterMacro( float, unsigned long );
    callCorrectReadDicomWriterMacro( float, long );
    callCorrectReadDicomWriterMacro( float, float );
    callCorrectReadDicomWriterMacro( float, double );

    /** From double to something else. */
    callCorrectReadDicomWriterMacro( double, unsigned char );
    callCorrectReadDicomWriterMacro( double, char );
    callCorrectReadDicomWriterMacro( double, unsigned short );
    callCorrectReadDicomWriterMacro( double, short );
    callCorrectReadDicomWriterMacro( double, unsigned int );
    callCorrectReadDicomWriterMacro( double, int );
    callCorrectReadDicomWriterMacro( double, unsigned long );
    callCorrectReadDicomWriterMacro( double, long );
    callCorrectReadDicomWriterMacro( double, float );
    callCorrectReadDicomWriterMacro( double, double );
    }
  else
    {
    std::cerr << "Dimension equals " << inputDimension << ", which is not supported." << std::endl;
    std::cerr << "Only 3D images are supported." << std::endl;
    return 1;
    } // end if over inputDimension

  /** Return a value. */
  return 0;
} // end support for SCALAR pixel type
