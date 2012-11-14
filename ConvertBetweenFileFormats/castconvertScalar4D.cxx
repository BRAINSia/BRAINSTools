/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    castconvertScalar4D.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "castconverthelpers.h"

int FileConverterScalar4D( const std::string & inputPixelComponentType,
                           const std::string & outputPixelComponentType, const std::string & inputFileName,
                           const std::string & outputFileName, int inputDimension )
{
  enum { ImageDims = 4 };

  if( inputDimension == ImageDims )
    {
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
