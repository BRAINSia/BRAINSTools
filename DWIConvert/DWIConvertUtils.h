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
#ifndef DWIConvertUtils_h
#define DWIConvertUtils_h
#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorRescaleIntensityImageFilter.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkVariableLengthVector.h"
#include <vcl_compiler.h>
#include <iostream>
#include "cmath"
#include <iostream>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <string>
#include "itkNumberToString.h"
#include "DWIMetaDataDictionaryValidator.h"

#include <cmath>

using PixelValueType = short;
/*
 * The Volume4DType is a scalar image with diemensions cols,rows,3Dslices,NumGradients
 */
using Volume4DType = itk::Image< PixelValueType, 4 >;

template < typename TArg >
int
CheckArg( const char * argName, const TArg & argVal, const TArg & emptyVal );

template < typename TImage >
int
WriteVolume( const TImage * img, const std::string & fname );

namespace itk
{
// helper specialization to figure out the component type of a VectorImage.
template < typename TPixel >
struct itk::ImageIOBase::MapPixelType< itk::VariableLengthVector< TPixel > >
{
  static const IOComponentType CType = itk::ImageIOBase::MapPixelType< TPixel >::CType;
};
} // namespace itk
template < typename TImage >
int
ReadScalarVolume( typename TImage::Pointer & img, const std::string & fname, bool allowLossyConversion );

template < typename TImage >
int
ReadVectorVolume( typename TImage::Pointer & img, const std::string & fname, bool allowLossyConversion );

template < typename TImage >
int
RecoverMeasurementFrame( const TImage * img, DWIMetaDataDictionaryValidator::RotationMatrixType & MeasurementFrame );

template < typename TImage >
int
RecoverBVectors( const TImage * img, DWIMetaDataDictionaryValidator::GradientTableType & bVecs );

template < typename TImage >
int
RecoverBValue( const TImage * img, double & val );

template < typename TImage >
int
RecoverBValues( const TImage * inputVol, const DWIMetaDataDictionaryValidator::GradientTableType & bVectors,
                std::vector< double > & bValues );

template < typename TScalar >
inline int
WriteBValues( const std::vector< TScalar > & bValues, const std::string & filename )
{
  // write out in FSL format
  std::ofstream bValFile;

  bValFile.open( filename.c_str(), std::ios::out | std::ios::binary );
  bValFile.precision( 17 );
  if ( !bValFile.is_open() || !bValFile.good() )
  {
    return EXIT_FAILURE;
  }

  for ( unsigned int k = 0; k < bValues.size(); ++k )
  {
    const char * const spacer = ( k == bValues.size() - 1 ) ? "" : " ";
    bValFile << bValues[k] << spacer;
  }
  bValFile << std::endl;
  bValFile.close();

  return EXIT_SUCCESS;
}

extern void
ConvertBvecsToFromFSL( DWIMetaDataDictionaryValidator::GradientTableType & bVecs );


extern void
normalize( const DWIMetaDataDictionaryValidator::GradientDirectionType & vec, double * normedVec );
extern int
WriteBVectors( const DWIMetaDataDictionaryValidator::GradientTableType & bVectors, const std::string & filename );

extern int
ReadBVals( std::vector< double > & bVals, unsigned int & bValCount, const std::string & bValFilename );

extern int
ReadBVecs( DWIMetaDataDictionaryValidator::GradientTableType & bVecs, unsigned int & bVecCount,
           const std::string & bVecFilename );

template < typename TValue >
bool
CloseEnough( const TValue & a, const TValue & b, double magdiv = 100000.0 );

template < typename TVal >
bool
CloseEnough( const std::vector< TVal > & a, const std::vector< TVal > & b, double magdiv = 100000.0 );

extern bool
CloseEnough( const vnl_vector_fixed< double, 3 > & a, const vnl_vector_fixed< double, 3 > & b,
             double magdiv = 100000.0 );

template < typename TVal >
bool
CloseEnough( const itk::VariableLengthVector< TVal > & a, const itk::VariableLengthVector< TVal > & b,
             double magdiv = 100000.0 );

template < typename TVal >
void
PrintVec( const TVal & a );

template < typename TVal >
void
PrintVec( const std::vector< TVal > & vec );

extern void
PrintVec( const vnl_vector_fixed< double, 3 > & vec );
extern void
PrintVec( const DWIMetaDataDictionaryValidator::GradientTableType & vec );

extern int
FSLToNrrd( const std::string & inputVolume, const std::string & outputVolume, const std::string & fslNIFTIFile,
           const std::string & inputBValues, const std::string & inputBVectors, bool horizontalBy3Rows,
           bool allowLossyConversion );

extern int
NrrdToFSL( const std::string & inputVolume, const std::string & outputVolume, const std::string & outputBValues,
           const std::string & outputBVectors, bool allowLossyConversion );


#include "DWIConvertUtils.hxx"

#endif // DWIConvertUtils_h
