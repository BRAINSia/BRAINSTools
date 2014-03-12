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
