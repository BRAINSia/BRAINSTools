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
/**
 * This file is to make some smaller compilation units to help improve compilation performance.
 */
#ifndef __BRAINSABCUtilities__h__
#define __BRAINSABCUtilities__h__

#include "Log.h"

#include <AtlasDefinition.h>
#include <BRAINSFitHelper.h>

#include <itkImage.h>
#include <itkArray.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>

#include <vector>
#include <map>
#include <csignal>

// The 200805 OpenMPv3.0 specificaiton allows unsigned iterators
#if defined(_OPENMP)
#define LOCAL_USE_OPEN_MP
#endif

#if defined(LOCAL_USE_OPEN_MP) && (_OPENMP < 200805)
typedef int LOOPITERTYPE;
#else
typedef unsigned int LOOPITERTYPE;
#endif

//  vnl_math_isnan(value) || vnl_math_isinf(value) )
#if 1  // Runtime performance penalty that can be used to find faulty code
       // during debugging.
#define CHECK_NAN( XXXTESTXXX, srcfile, srcline, extra_print ) \
    { \
    if( !vnl_math_isfinite( XXXTESTXXX ) ) \
      { \
      std::cout << "Found " << XXXTESTXXX << " at " << srcfile << " " << srcline << extra_print << std::endl; \
      raise(SIGSEGV); \
      } \
    }
#else
#define CHECK_NAN( XXXTESTXXX, srcfile, srcline, extra_print ) \
{ }
#endif

typedef double                          FloatingPrecision;
typedef itk::Image<unsigned char, 3>    ByteImageType;
typedef itk::Image<float, 3>            FloatImageType;
typedef FloatImageType::Pointer         FloatImagePointerType;
typedef itk::Image<signed short int, 3> ShortImageType;

typedef itk::Image<float, 3> CorrectIntensityImageType;

typedef std::map<std::string,std::string> ImageByTypeMap;

/** A utiliy class for holding statistical information
 * for all image channels for a given tissue class type
 */
class RegionStats
{
public:
  typedef vnl_matrix<FloatingPrecision>         MatrixType;
  typedef vnl_matrix_inverse<FloatingPrecision> MatrixInverseType;
  typedef std::map<std::string,double>          MeanMapType;;
  RegionStats() : m_Means(), m_Covariance(), m_Weighting(0.0)
  {
  }

  void resize(const unsigned int numModalities)
  {
    this->m_Covariance = MatrixType(numModalities, numModalities);
    this->m_Means.clear();
  }

  MeanMapType m_Means;             // One measure per image channel type;
  MatrixType  m_Covariance;        // Matrix of covariances of class by image
                                  // channel
  FloatingPrecision m_Weighting;  // The strength of this class.
};

#include "BRAINSABCUtilities.hxx"

// External Templates to improve compilation times.
extern std::vector<CorrectIntensityImageType::Pointer>
CorrectBias(const unsigned int degree,
            const unsigned int CurrentEMIteration,
            const std::vector<ByteImageType::Pointer> & CandidateRegions,
            const std::vector<CorrectIntensityImageType::Pointer> & inputImages,
            const ByteImageType::Pointer currentBrainMask,
            const ByteImageType::Pointer currentForegroundMask,
            const std::vector<FloatImageType::Pointer> & probImages,
            const std::vector<bool> & probUseForBias,
            const FloatingPrecision sampleSpacing, const int DebugLevel,
            const std::string& OutputDebugDir);

extern template std::vector<FloatImagePointerType> DuplicateImageList<FloatImageType>(
  const std::vector<FloatImagePointerType> & );

extern template std::vector<ShortImageType::Pointer> DuplicateImageList<ShortImageType>(
  const std::vector<ShortImageType::Pointer> & );

extern template void ComputeLabels<FloatImageType,
                                   ByteImageType,
                                   double>( std::vector<FloatImageType::Pointer> &, std::vector<bool> &,
                                            const vnl_vector<unsigned int> &, ByteImageType::Pointer &,
                                            ByteImageType::Pointer &,
                                            ByteImageType::Pointer &,
                                            FloatingPrecision );

extern template void NormalizeProbListInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> & );

extern template void ZeroNegativeValuesInPlace<FloatImageType>(  std::vector<FloatImageType::Pointer> & );

template <class TMap>
unsigned int TotalMapSize(const TMap &map)
{
  unsigned int rval = 0;
  for(typename TMap::const_iterator mapIt = map.begin();
      mapIt != map.end(); ++mapIt)
    {
    for(typename TMap::mapped_type::const_iterator listIt =
          mapIt->second.begin(); listIt != mapIt->second.end(); ++listIt)
      {
      ++rval;
      }
    }
  return rval;
}

template <class TMap>
typename TMap::mapped_type::value_type &
GetMapVectorNthElement(TMap &map, int n)
{
  typename TMap::mapped_type::value_type returnElement;
  if( map.size() < n )
  {
  returnElement = NULL;
  }
  else
  {
  typename TMap::iterator it = map.begin();
  for( int i = 0;
       i < n, it != map.end();
       i++, ++ it)
    {
    if( i == n-1)
      {
      returnElement = *(map.begin()->second.begin());
      }
    }
  }
  return returnElement;
}

template <class TMap>
typename TMap::mapped_type::value_type &
GetMapVectorFirstElement(TMap &map)
{
  return *(map.begin()->second.begin());
}

template <class TMap>
const typename TMap::mapped_type::value_type &
GetMapVectorFirstElement(const TMap &map)
{
  return *(map.begin()->second.begin());
}

template <class ImageType>
typename ImageType::Pointer
CopyImage(const typename ImageType::Pointer & input )
{
  typedef itk::ImageDuplicator<ImageType> ImageDupeType;
  typename ImageDupeType::Pointer MyDuplicator = ImageDupeType::New();
  MyDuplicator->SetInputImage(input);
  MyDuplicator->Update();
  return MyDuplicator->GetOutput();
}


#endif // __BRAINSABCUtilities__h__
