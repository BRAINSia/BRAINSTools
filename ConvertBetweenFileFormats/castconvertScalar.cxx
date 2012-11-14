/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    castconvertScalar.cxx
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

extern int FileConverterScalar2D( const std::string & inputPixelComponentType,
                                  const std::string & outputPixelComponentType, const std::string & inputFileName,
                                  const std::string & outputFileName, int inputDimension );

extern int FileConverterScalar3D( const std::string & inputPixelComponentType,
                                  const std::string & outputPixelComponentType, const std::string & inputFileName,
                                  const std::string & outputFileName, int inputDimension );

extern int FileConverterScalar2DA( const std::string & inputPixelComponentType,
                                   const std::string & outputPixelComponentType, const std::string & inputFileName,
                                   const std::string & outputFileName, int inputDimension );

extern int FileConverterScalar3DA( const std::string & inputPixelComponentType,
                                   const std::string & outputPixelComponentType, const std::string & inputFileName,
                                   const std::string & outputFileName, int inputDimension );

extern int FileConverterScalar4D( const std::string & inputPixelComponentType,
                                  const std::string & outputPixelComponentType, const std::string & inputFileName,
                                  const std::string & outputFileName, int inputDimension );

extern int FileConverterScalar4DA( const std::string & inputPixelComponentType,
                                   const std::string & outputPixelComponentType, const std::string & inputFileName,
                                   const std::string & outputFileName, int inputDimension );

int FileConverterScalar( const std::string & inputPixelComponentType,
                         const std::string & outputPixelComponentType, const std::string & inputFileName,
                         const std::string & outputFileName, int inputDimension )
{
  /** Support for 2D images. */
  if( inputDimension == 2 )
    {
    const int ret_value = FileConverterScalar2D(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension )
      ||                  FileConverterScalar2DA(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension );
    if( ret_value != 0 )
      {
      return ret_value;
      }
    }
  else if( inputDimension == 3 )
    {
    const int ret_value = FileConverterScalar3D(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension )
      ||                  FileConverterScalar3DA(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension );
    if( ret_value != 0 )
      {
      return ret_value;
      }
    } // end support for 3D images
  else if( inputDimension == 4 )
    {
    const int ret_value = FileConverterScalar4D(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension )
      ||                  FileConverterScalar4DA(
        inputPixelComponentType, outputPixelComponentType,
        inputFileName, outputFileName, inputDimension );
    if( ret_value != 0 )
      {
      return ret_value;
      }
    } // end support for 4D images
  else
    {
    std::cerr << "Dimension equals " << inputDimension << ", which is not supported." << std::endl;
    std::cerr << "Only 2D, 3D, and 4D images are supported." << std::endl;
    return 1;
    } // end if over inputDimension

  /** Return a succes value. */
  return 0;
} // end support for SCALAR pixel type
