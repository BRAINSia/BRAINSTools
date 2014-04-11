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

int FileConverterScalar2D( const std::string & inputPixelComponentType,
                           const std::string & outputPixelComponentType, const std::string & inputFileName,
                           const std::string & outputFileName, int inputDimension )
{
  enum { ImageDims = 2 };

  if( inputDimension == ImageDims )
    {
    /** From unsigned char to something else. */
    callCorrectReadWriterMacro( unsigned char, unsigned char, ImageDims );
    callCorrectReadWriterMacro( unsigned char, char, ImageDims );
    callCorrectReadWriterMacro( unsigned char, unsigned short, ImageDims );
    callCorrectReadWriterMacro( unsigned char, short, ImageDims );
    callCorrectReadWriterMacro( unsigned char, unsigned int, ImageDims );
    callCorrectReadWriterMacro( unsigned char, int, ImageDims );
    callCorrectReadWriterMacro( unsigned char, unsigned long, ImageDims );
    callCorrectReadWriterMacro( unsigned char, long, ImageDims );
    callCorrectReadWriterMacro( unsigned char, float, ImageDims );
    callCorrectReadWriterMacro( unsigned char, double, ImageDims );

    /** From char to something else. */
    callCorrectReadWriterMacro( char, unsigned char, ImageDims );
    callCorrectReadWriterMacro( char, char, ImageDims );
    callCorrectReadWriterMacro( char, unsigned short, ImageDims );
    callCorrectReadWriterMacro( char, short, ImageDims );
    callCorrectReadWriterMacro( char, unsigned int, ImageDims );
    callCorrectReadWriterMacro( char, int, ImageDims );
    callCorrectReadWriterMacro( char, unsigned long, ImageDims );
    callCorrectReadWriterMacro( char, long, ImageDims );
    callCorrectReadWriterMacro( char, float, ImageDims );
    callCorrectReadWriterMacro( char, double, ImageDims );

    /** From unsigned short to something else. */
    callCorrectReadWriterMacro( unsigned short, unsigned char, ImageDims );
    callCorrectReadWriterMacro( unsigned short, char, ImageDims );
    callCorrectReadWriterMacro( unsigned short, unsigned short, ImageDims );
    callCorrectReadWriterMacro( unsigned short, short, ImageDims );
    callCorrectReadWriterMacro( unsigned short, unsigned int, ImageDims );
    callCorrectReadWriterMacro( unsigned short, int, ImageDims );
    callCorrectReadWriterMacro( unsigned short, unsigned long, ImageDims );
    callCorrectReadWriterMacro( unsigned short, long, ImageDims );
    callCorrectReadWriterMacro( unsigned short, float, ImageDims );
    callCorrectReadWriterMacro( unsigned short, double, ImageDims );

    /** From short to something else. */
    callCorrectReadWriterMacro( short, unsigned char, ImageDims );
    callCorrectReadWriterMacro( short, char, ImageDims );
    callCorrectReadWriterMacro( short, unsigned short, ImageDims );
    callCorrectReadWriterMacro( short, short, ImageDims );
    callCorrectReadWriterMacro( short, unsigned int, ImageDims );
    callCorrectReadWriterMacro( short, int, ImageDims );
    callCorrectReadWriterMacro( short, unsigned long, ImageDims );
    callCorrectReadWriterMacro( short, long, ImageDims );
    callCorrectReadWriterMacro( short, float, ImageDims );
    callCorrectReadWriterMacro( short, double, ImageDims );

    /** From unsigned int to something else. */
    callCorrectReadWriterMacro( unsigned int, unsigned char, ImageDims );
    callCorrectReadWriterMacro( unsigned int, char, ImageDims );
    callCorrectReadWriterMacro( unsigned int, unsigned short, ImageDims );
    callCorrectReadWriterMacro( unsigned int, short, ImageDims );
    callCorrectReadWriterMacro( unsigned int, unsigned int, ImageDims );
    callCorrectReadWriterMacro( unsigned int, int, ImageDims );
    callCorrectReadWriterMacro( unsigned int, unsigned long, ImageDims );
    callCorrectReadWriterMacro( unsigned int, long, ImageDims );
    callCorrectReadWriterMacro( unsigned int, float, ImageDims );
    callCorrectReadWriterMacro( unsigned int, double, ImageDims );
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
