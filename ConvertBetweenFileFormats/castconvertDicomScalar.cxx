/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    castconvertDicomScalar.cxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/** Authors: Marius Staring and Stefan Klein **/
#include "castconverthelpers.h"

int DicomFileConverterScalar( const std::string & inputPixelComponentType,
                              const std::string & outputPixelComponentType, const std::string & inputDirectoryName,
                              const std::string & outputFileName, int inputDimension )
{
  /** Support for 3D images. */
  if( inputDimension == 3 )
    {
    /** From unsigned char to something else. */
    callCorrectReadDicomWriterMacro( unsigned char, unsigned char );
    callCorrectReadDicomWriterMacro( unsigned char, char );
    callCorrectReadDicomWriterMacro( unsigned char, unsigned short );
    callCorrectReadDicomWriterMacro( unsigned char, short );
    callCorrectReadDicomWriterMacro( unsigned char, unsigned int );
    callCorrectReadDicomWriterMacro( unsigned char, int );
    callCorrectReadDicomWriterMacro( unsigned char, unsigned long );
    callCorrectReadDicomWriterMacro( unsigned char, long );
    callCorrectReadDicomWriterMacro( unsigned char, float );
    callCorrectReadDicomWriterMacro( unsigned char, double );

    /** From char to something else. */
    callCorrectReadDicomWriterMacro( char, unsigned char );
    callCorrectReadDicomWriterMacro( char, char );
    callCorrectReadDicomWriterMacro( char, unsigned short );
    callCorrectReadDicomWriterMacro( char, short );
    callCorrectReadDicomWriterMacro( char, unsigned int );
    callCorrectReadDicomWriterMacro( char, int );
    callCorrectReadDicomWriterMacro( char, unsigned long );
    callCorrectReadDicomWriterMacro( char, long );
    callCorrectReadDicomWriterMacro( char, float );
    callCorrectReadDicomWriterMacro( char, double );

    /** From unsigned short to something else. */
    callCorrectReadDicomWriterMacro( unsigned short, unsigned char );
    callCorrectReadDicomWriterMacro( unsigned short, char );
    callCorrectReadDicomWriterMacro( unsigned short, unsigned short );
    callCorrectReadDicomWriterMacro( unsigned short, short );
    callCorrectReadDicomWriterMacro( unsigned short, unsigned int );
    callCorrectReadDicomWriterMacro( unsigned short, int );
    callCorrectReadDicomWriterMacro( unsigned short, unsigned long );
    callCorrectReadDicomWriterMacro( unsigned short, long );
    callCorrectReadDicomWriterMacro( unsigned short, float );
    callCorrectReadDicomWriterMacro( unsigned short, double );

    /** From short to something else. */
    callCorrectReadDicomWriterMacro( short, unsigned char );
    callCorrectReadDicomWriterMacro( short, char );
    callCorrectReadDicomWriterMacro( short, unsigned short );
    callCorrectReadDicomWriterMacro( short, short );
    callCorrectReadDicomWriterMacro( short, unsigned int );
    callCorrectReadDicomWriterMacro( short, int );
    callCorrectReadDicomWriterMacro( short, unsigned long );
    callCorrectReadDicomWriterMacro( short, long );
    callCorrectReadDicomWriterMacro( short, float );
    callCorrectReadDicomWriterMacro( short, double );

    /** From unsigned int to something else. */
    callCorrectReadDicomWriterMacro( unsigned int, unsigned char );
    callCorrectReadDicomWriterMacro( unsigned int, char );
    callCorrectReadDicomWriterMacro( unsigned int, unsigned short );
    callCorrectReadDicomWriterMacro( unsigned int, short );
    callCorrectReadDicomWriterMacro( unsigned int, unsigned int );
    callCorrectReadDicomWriterMacro( unsigned int, int );
    callCorrectReadDicomWriterMacro( unsigned int, unsigned long );
    callCorrectReadDicomWriterMacro( unsigned int, long );
    callCorrectReadDicomWriterMacro( unsigned int, float );
    callCorrectReadDicomWriterMacro( unsigned int, double );
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
