#ifndef DWIConvertUtils_h
#define DWIConvertUtils_h
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkVariableLengthVector.h"
#include "vcl_cmath.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

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

template <typename TImage>
int
ReadVolume( typename TImage::Pointer & img, const std::string & fname )
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
  img = imgReader->GetOutput();
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBVectors(const TImage *img, std::vector<std::vector<double> > & bVecs)
{
  bVecs.clear();

  const itk::MetaDataDictionary & dict = img->GetMetaDataDictionary();
  for( unsigned curGradientVec = 0; ; ++curGradientVec )
    {
    std::stringstream labelSS;
    labelSS << "DWMRI_gradient_"
            << std::setw(4)
            << std::setfill('0')
            << curGradientVec;
    std::string valString;
    // look for gradients in metadata until none by current name exists
    if( !itk::ExposeMetaData<std::string>(dict, labelSS.str(), valString) )
      {
      break;
      }
    std::stringstream   valSS(valString);
    std::vector<double> vec;
    for( ; ; )
      {
      double curVal;
      valSS >> curVal;
      if( !valSS.fail() )
        {
        vec.push_back(curVal);
        }
      else
        {
        break;
        }
      }
    bVecs.push_back(vec);
    }
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
  std::string valString;

  const itk::MetaDataDictionary & dict = img->GetMetaDataDictionary();

  if( !itk::ExposeMetaData<std::string>(dict, "DWMRI_b-value", valString) )
    {
    return EXIT_FAILURE;
    }

  std::stringstream valSS(valString);

  valSS >> val;

  if( valSS.fail() )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

template <typename TImage>
int RecoverBValues(const TImage *inputVol,
                   const std::vector<std::vector<double> > & bVectors,
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
    double bval = vcl_sqrt( (bVectors[i][0] * bVectors[i][0])
                            + (bVectors[i][1] * bVectors[i][1])
                            + (bVectors[i][2] * bVectors[i][2]) );
    bval *= BValue;
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
    if( !bValFile.good() )
      {
      return EXIT_FAILURE;
      }
    bValFile << bValues[k] << std::endl;
    }

  bValFile.close();

  return EXIT_SUCCESS;
}

template <typename TScalar>
inline int
WriteBVectors(const std::vector<std::vector<TScalar> > & bVectors,
              const std::string & filename)
{
  std::ofstream bVecFile;

  bVecFile.open(filename.c_str(), std::ios::out | std::ios::binary);
  bVecFile.precision(17);
  if( !bVecFile.is_open() || !bVecFile.good() )
    {
    return EXIT_FAILURE;
    }
  for( unsigned int k = 0; k < bVectors.size(); ++k )
    {
    if( !bVecFile.good() )
      {
      return EXIT_FAILURE;
      }
    bVecFile << bVectors[k][0] << " "
             << bVectors[k][1] << " "
             << bVectors[k][2]
             << std::endl;
    }
  bVecFile.close();
  return EXIT_SUCCESS;
}

inline
int
ReadBVals(std::vector<double> & bVals, unsigned int & bValCount, const std::string & bValFilename, double & maxBValue)
{
  std::ifstream bValFile(bValFilename.c_str(), std::ifstream::in);

  if( !bValFile.good() )
    {
    std::cerr << "Failed to open " << bValFilename
              << std::endl;
    return EXIT_FAILURE;
    }
  bVals.clear();
  bValCount = 0;
  while( !bValFile.eof() )
    {
    double x;
    bValFile >> x;
    if( bValFile.fail() )
      {
      break;
      }
    if( x > maxBValue )
      {
      maxBValue = x;
      }
    bValCount++;
    bVals.push_back(x);
    }

  return EXIT_SUCCESS;
}

inline
int
ReadBVecs(std::vector<std::vector<double> > & bVecs, unsigned int & bVecCount, const std::string & bVecFilename)
{
  std::ifstream bVecFile(bVecFilename.c_str(), std::ifstream::in);

  if( !bVecFile.good() )
    {
    std::cerr << "Failed to open " << bVecFilename
              << std::endl;
    return EXIT_FAILURE;
    }
  bVecs.clear();
  bVecCount = 0;
  while( !bVecFile.eof() )
    {
    std::vector<double> x;
    for( unsigned i = 0; i < 3; ++i )
      {
      double val;
      bVecFile >> val;
      if( bVecFile.fail() )
        {
        break;
        }
      x.push_back(val);
      }
    if( bVecFile.fail() )
      {
      break;
      }
    bVecCount++;
    bVecs.push_back(x);
    }

  return EXIT_SUCCESS;
}

template <typename TValue>
bool
CloseEnough(const TValue & a, const TValue & b, double magdiv = 100000.0)
{
  double averageMag = (vcl_fabs(static_cast<double>(a) )
                       + vcl_fabs(static_cast<double>(b) ) ) / 2.0;
  double diff = vcl_fabs(static_cast<double>(a) - static_cast<double>(b) );

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

#endif // DWIConvertUtils_h
