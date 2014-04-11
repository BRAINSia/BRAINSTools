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
#include "castconverthelpers.h"

int FileConverterScalar2DA( const std::string & inputPixelComponentType,
                            const std::string & outputPixelComponentType, const std::string & inputFileName,
                            const std::string & outputFileName, int inputDimension )
{
  enum { ImageDims = 2 };

  if( inputDimension == ImageDims )
    {
    /** From int to something else. */
    callCorrectReadWriterMacro( int, unsigned char, ImageDims );
    callCorrectReadWriterMacro( int, char, ImageDims );
    callCorrectReadWriterMacro( int, unsigned short, ImageDims );
    callCorrectReadWriterMacro( int, short, ImageDims );
    callCorrectReadWriterMacro( int, unsigned int, ImageDims );
    callCorrectReadWriterMacro( int, int, ImageDims );
    callCorrectReadWriterMacro( int, unsigned long, ImageDims );
    callCorrectReadWriterMacro( int, long, ImageDims );
    callCorrectReadWriterMacro( int, float, ImageDims );
    callCorrectReadWriterMacro( int, double, ImageDims );

    /** From unsigned long to something else. */
    callCorrectReadWriterMacro( unsigned long, unsigned char, ImageDims );
    callCorrectReadWriterMacro( unsigned long, char, ImageDims );
    callCorrectReadWriterMacro( unsigned long, unsigned short, ImageDims );
    callCorrectReadWriterMacro( unsigned long, short, ImageDims );
    callCorrectReadWriterMacro( unsigned long, unsigned int, ImageDims );
    callCorrectReadWriterMacro( unsigned long, int, ImageDims );
    callCorrectReadWriterMacro( unsigned long, unsigned long, ImageDims );
    callCorrectReadWriterMacro( unsigned long, long, ImageDims );
    callCorrectReadWriterMacro( unsigned long, float, ImageDims );
    callCorrectReadWriterMacro( unsigned long, double, ImageDims );

    /** From long to something else. */
    callCorrectReadWriterMacro( long, unsigned char, ImageDims );
    callCorrectReadWriterMacro( long, char, ImageDims );
    callCorrectReadWriterMacro( long, unsigned short, ImageDims );
    callCorrectReadWriterMacro( long, short, ImageDims );
    callCorrectReadWriterMacro( long, unsigned int, ImageDims );
    callCorrectReadWriterMacro( long, int, ImageDims );
    callCorrectReadWriterMacro( long, unsigned long, ImageDims );
    callCorrectReadWriterMacro( long, long, ImageDims );
    callCorrectReadWriterMacro( long, float, ImageDims );
    callCorrectReadWriterMacro( long, double, ImageDims );

    /** From float to something else. */
    callCorrectReadWriterMacro( float, unsigned char, ImageDims );
    callCorrectReadWriterMacro( float, char, ImageDims );
    callCorrectReadWriterMacro( float, unsigned short, ImageDims );
    callCorrectReadWriterMacro( float, short, ImageDims );
    callCorrectReadWriterMacro( float, unsigned int, ImageDims );
    callCorrectReadWriterMacro( float, int, ImageDims );
    callCorrectReadWriterMacro( float, unsigned long, ImageDims );
    callCorrectReadWriterMacro( float, long, ImageDims );
    callCorrectReadWriterMacro( float, float, ImageDims );
    callCorrectReadWriterMacro( float, double, ImageDims );

    /** From double to something else. */
    callCorrectReadWriterMacro( double, unsigned char, ImageDims );
    callCorrectReadWriterMacro( double, char, ImageDims );
    callCorrectReadWriterMacro( double, unsigned short, ImageDims );
    callCorrectReadWriterMacro( double, short, ImageDims );
    callCorrectReadWriterMacro( double, unsigned int, ImageDims );
    callCorrectReadWriterMacro( double, int, ImageDims );
    callCorrectReadWriterMacro( double, unsigned long, ImageDims );
    callCorrectReadWriterMacro( double, long, ImageDims );
    callCorrectReadWriterMacro( double, float, ImageDims );
    callCorrectReadWriterMacro( double, double, ImageDims );
    }
  else
    {
    std::cerr << "Dimension equals " << inputDimension << ", which is not supported." << std::endl;
    std::cerr << "Only " << ImageDims << "D images are supported." << std::endl;
    return 1;
    } // end if over inputDimension

  /** Return a value. */
  return 0;
} // end support for SCALAR pixel type
