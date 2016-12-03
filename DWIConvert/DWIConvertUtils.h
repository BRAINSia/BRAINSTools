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

typedef short                          PixelValueType;
/*
 * The Volume4DType is a scalar image with diemensions cols,rows,3Dslices,NumGradients
 */
typedef itk::Image<PixelValueType, 4>  Volume4DType;



template <typename TArg>
int
CheckArg(const char *argName, const TArg & argVal, const TArg & emptyVal)
{
  if( argVal == emptyVal )
    {
    std::cerr << "Missing argument " << argName << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
WriteVolume( const TImage *img, const std::string & fname )
{
  typename itk::ImageFileWriter<TImage>::Pointer imgWriter =
    itk::ImageFileWriter<TImage>::New();

  imgWriter->SetInput( img );
  imgWriter->SetFileName( fname.c_str() );
  try
    {
    imgWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing "
              << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

namespace itk {
  // helper specialization to figure out the component type of a VectorImage.
  template <typename TPixel> struct itk::ImageIOBase::MapPixelType< itk::VariableLengthVector<TPixel> >
  {
    static const IOComponentType CType = itk::ImageIOBase::MapPixelType<TPixel>::CType;
  };
}

template <typename TImage>
int
ReadVolume( typename TImage::Pointer & img, const std::string & fname, bool allowLossyConversion = false )
{
  typename itk::ImageFileReader<TImage>::Pointer imgReader =
    itk::ImageFileReader<TImage>::New();

  imgReader->SetFileName( fname.c_str() );
  try
    {
    imgReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while reading "
              << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  if (!allowLossyConversion)
    {
    itk::ImageIOBase *imageIO = imgReader->GetImageIO();
    itk::ImageIOBase::IOComponentType ioType =
      itk::ImageIOBase::MapPixelType< typename TImage::PixelType >::CType;

    if (imageIO->GetComponentType() != ioType)
      {
      std::cerr << "Error: ReadVolume: Unsupported source pixel type." << std:: endl
                << "  Input volume:  " << imageIO->GetComponentTypeAsString(imageIO->GetComponentType())
                << std::endl
                << "  Output volume: " << imageIO->GetComponentTypeAsString(ioType)
                << std::endl
                << "The only supported output type is <short>. "
                << "You may consider using allowLossyConversion option."
                << std::endl
                << "However, use this option with caution! "
                << "Conversion from images of a different type may cause data loss due to rounding or truncation."
                << std::endl;
      return EXIT_FAILURE;
      }
    }

  img = imgReader->GetOutput();
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverMeasurementFrame(const TImage *img, DWIMetaDataDictionaryValidator::RotationMatrixType & MeasurementFrame)
{

  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary() );
  MeasurementFrame = myDWIValidator.GetMeasurementFrame();
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBVectors(const TImage *img, DWIMetaDataDictionaryValidator::GradientTableType & bVecs)
{
  bVecs.clear();

  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary() );
  bVecs = myDWIValidator.GetGradientTable();
  if( bVecs.empty() )
  {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBValue(const TImage *img, double & val)
{
  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary() );
  val = myDWIValidator.GetBValue();
  return EXIT_SUCCESS;
}

template <typename TImage>
int RecoverBValues(const TImage *inputVol,
                   const DWIMetaDataDictionaryValidator::GradientTableType & bVectors,
                   std::vector<double> & bValues)
{
  bValues.clear();

  double BValue;
  if( RecoverBValue(inputVol, BValue) != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  for( unsigned i = 0; i < bVectors.size(); ++i )
    {
    double norm = std::sqrt( (bVectors[i][0] * bVectors[i][0])
                            + (bVectors[i][1] * bVectors[i][1])
                            + (bVectors[i][2] * bVectors[i][2]) );
    if( std::abs( 1- norm) < 1e-4 ) // Asssume value very close to 1 are 1
      {
      norm = 1.0;
      }
    // bval_i = (G_norm)^2 * bval_max
    double bval = norm*norm*BValue;
    if( std::abs( bval - itk::Math::Round<double>(bval) ) < 1e-2 )
      {
      bval = itk::Math::Round<double>(bval);
      }
    bValues.push_back(bval);
    }
  return EXIT_SUCCESS;
}

template <typename TScalar>
inline int
WriteBValues(const std::vector<TScalar> & bValues, const std::string & filename)
{
  // write out in FSL format
  std::ofstream bValFile;

  bValFile.open(filename.c_str(), std::ios::out | std::ios::binary);
  bValFile.precision(17);
  if( !bValFile.is_open() || !bValFile.good() )
    {
    return EXIT_FAILURE;
    }

  for( unsigned int k = 0; k < bValues.size(); ++k )
  {
    const char * const spacer = (  k==bValues.size()-1 ) ? "" : " ";
    bValFile << bValues[k] << spacer;
  }
  bValFile << std::endl;
  bValFile.close();

  return EXIT_SUCCESS;
}

extern void normalize(const DWIMetaDataDictionaryValidator::GradientDirectionType &vec,double *normedVec);
extern int WriteBVectors(const DWIMetaDataDictionaryValidator::GradientTableType & bVectors,
              const std::string & filename);

extern int ReadBVals(std::vector<double> & bVals, unsigned int & bValCount,
                     const std::string & bValFilename);

extern int ReadBVecs(DWIMetaDataDictionaryValidator::GradientTableType & bVecs, unsigned int & bVecCount,
                     const std::string & bVecFilename , bool horizontalBy3Rows );



template <typename TValue>
bool
CloseEnough(const TValue & a, const TValue & b, double magdiv = 100000.0)
{
  double averageMag = (std::fabs(static_cast<double>(a) )
                       + std::fabs(static_cast<double>(b) ) ) / 2.0;
  double diff = std::fabs(static_cast<double>(a) - static_cast<double>(b) );

  // case one -- both near zero
  if( averageMag < 0.000001 )
    {
    return true;
    }
  // case 2 -- diff > average / 100000;
  if( diff > (averageMag / magdiv) )
    {
    return false;
    }
  return true;
}

template <typename TVal>
bool
CloseEnough(const std::vector<TVal> & a, const std::vector<TVal> & b, double magdiv = 100000.0)
{
  if( a.size() != b.size() )
    {
    std::cerr << "Vector size mismatch: "
              << a.size() << " "
              << b.size() << std::endl;
    return false;
    }
  for( unsigned i = 0; i < a.size(); ++i )
    {
    if( !CloseEnough(a[i], b[i], magdiv) )
      {
      std::cerr << "Value mismatch" << std::endl;
      return false;
      }
    }
  return true;
}

extern bool CloseEnough(const vnl_vector_fixed<double,3> & a, const vnl_vector_fixed<double,3> & b, double magdiv = 100000.0);

template <typename TVal>
bool
CloseEnough(const itk::VariableLengthVector<TVal> & a,
            const itk::VariableLengthVector<TVal> & b, double magdiv = 100000.0)
{
  if( a.GetSize() != b.GetSize() )
    {
    std::cerr << "Vector size mismatch: "
              << a.GetSize() << " "
              << b.GetSize() << std::endl;
    return false;
    }
  for( unsigned i = 0; i < a.GetSize(); ++i )
    {
    if( !CloseEnough(a[i], b[i], magdiv) )
      {
      std::cerr << "Value mismatch" << std::endl;
      return false;
      }
    }
  return true;
}


template <typename TVal>
void PrintVec(const TVal & a)
{
  std::cerr << a;
}

template <typename TVal>
void PrintVec(const std::vector<TVal> & vec)
{
  std::cerr << "[";

  for( unsigned i = 0; i < vec.size(); ++i )
    {
    PrintVec(vec[i]);
    if( i < vec.size() - 1 )
      {
      std::cerr << " ";
      }
    }
  std::cerr << "]" << std::endl;
}

extern void PrintVec(const vnl_vector_fixed<double,3> & vec);
extern void PrintVec(const DWIMetaDataDictionaryValidator::GradientTableType & vec);

extern int FSLToNrrd(const std::string & inputVolume,
                     const std::string & outputVolume,
                     const std::string & fslNIFTIFile,
                     const std::string & inputBValues,
                     const std::string & inputBVectors,
                     bool horizontalBy3Rows,
                     bool allowLossyConversion);

    extern int NrrdToFSL(const std::string & inputVolume,
                         const std::string & outputVolume,
                         const std::string & outputBValues,
                         const std::string & outputBVectors,
                         bool allowLossyConversion);
#endif // DWIConvertUtils_h
