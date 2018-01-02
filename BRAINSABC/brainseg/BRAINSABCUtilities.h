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
#include <algorithm>

typedef unsigned int LOOPITERTYPE;

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

/** A utility class for holding ordered maps */

/**
 * Make ordered map by choosing less than operator based on order of inputs
 * This comparitor function will order the maps based on a first in
 * priority. Elements that are not already in the list of keys
 * are added at the end.  Elements that show up in the list first are considered
 * less than elements that show up in the list second.
 **/
class firstInOrderingOfStrings
{
public:
  firstInOrderingOfStrings() {
  //HACK: This is hard coding that should not be here.
  //      A better insertion ordered map is needed.
  this->m_firstInOrdering.push_back("T1");
  this->m_firstInOrdering.push_back("T2");
  this->m_firstInOrdering.push_back("PD");
  this->m_firstInOrdering.push_back("FLAIR");
  this->m_firstInOrdering.push_back("OTHER");
  }
  bool operator() (const std::string& lhs, const std::string& rhs) const
  {
  add_keys_internally(lhs); //For new items, add lhs first
  add_keys_internally(rhs); //Also need to add rhs if it is not in list yet
  for(auto &elem : m_firstInOrdering)
    {
    if(elem == lhs && (lhs != rhs ))
      {
      return true;
      }
    else if(elem == rhs)
      {
      return false;
      }
    }
    throw std::invalid_argument( "Values need to be added before comparisions" );
    return false;
  }
private:
  /*
   *  Need to ensure that items are in list before comparisons.
   */
  void add_keys_internally(const std::string & testValue) const
  {
      const bool found = (
              std::find(m_firstInOrdering.begin(), m_firstInOrdering.end(),
                  testValue) != m_firstInOrdering.end()
              );
      if ( ! found )
      {
          this->m_firstInOrdering.push_back(testValue);
      }
  }
  /* This needs to be mutable. The function signature needs to be
   * const, but this internal list needs to be mutable for the case
   * of new elements being added*/
  mutable std::list<std::string> m_firstInOrdering;
};

template<class Key, class T>
class orderedmap: public std::map<Key, T, firstInOrderingOfStrings>
{
public:
  orderedmap():std::map<Key, T, firstInOrderingOfStrings>()
  {};
  using std::map<Key, T, firstInOrderingOfStrings>::operator[];
};

typedef double                          FloatingPrecision;
typedef itk::Image<unsigned char, 3>    ByteImageType;
typedef itk::Image<float, 3>            FloatImageType;
typedef FloatImageType::Pointer         FloatImagePointerType;
typedef itk::Image<signed short int, 3> ShortImageType;

typedef itk::Image<float, 3> CorrectIntensityImageType;

typedef orderedmap<std::string,std::string> ImageByTypeMap;

typedef std::vector<FloatImagePointerType>        FloatImageVector;
typedef orderedmap<std::string, FloatImageVector> MapOfFloatImageVectors;

typedef itk::Transform<double, 3, 3>               GenericTransformType;
typedef std::vector<GenericTransformType::Pointer> TransformList;
typedef orderedmap<std::string,TransformList>      MapOfTransformLists;

/** A utiliy class for holding statistical information
 * for all image channels for a given tissue class type
 */
class RegionStats
{
public:
  typedef vnl_matrix<FloatingPrecision>         MatrixType;
  typedef vnl_matrix_inverse<FloatingPrecision> MatrixInverseType;
  typedef orderedmap<std::string,double>        MeanMapType;
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

/*
 * This function gets a map of input images, finds its first key image,
 * and resamples all images to the first key image lattice using identity
 * transform and passed interpolation type.
 * Note that it is assumed that all input intensity images are already aligned
 * in physical space.
 */
extern MapOfFloatImageVectors
ResampleImageListToFirstKeyImage(const std::string & resamplerInterpolatorType,
                                 const MapOfFloatImageVectors & inputImageMap);

/*
 * This function, first, transforms all inputImageMap to the space of the first image of the map
 * using rigid transforms (intraSubjectTransforms) and Resampling InPlace interoplation.
 * Then, it resamples all images within one modality to the voxel lattice of the fist image of that modality channel
 * using resamplerInterpolatorType and Identity transform.
 */
extern MapOfFloatImageVectors
ResampleInPlaceImageList(const std::string & resamplerInterpolatorType,
                         const MapOfFloatImageVectors & inputImageMap,
                         MapOfTransformLists & intraSubjectTransforms);

extern template std::vector<FloatImagePointerType> DuplicateImageList<FloatImageType>(
  const std::vector<FloatImagePointerType> & );

extern template std::vector<ShortImageType::Pointer> DuplicateImageList<ShortImageType>(
  const std::vector<ShortImageType::Pointer> & );

extern template void NormalizeProbListInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> & );

extern template void ZeroNegativeValuesInPlace<FloatImageType>(  std::vector<FloatImageType::Pointer> & );

template <class ImageType>
typename ImageType::Pointer
NormalizeInputIntensityImage(const typename ImageType::Pointer inputImage)
{
  muLogMacro(<< "\nNormalize input intensity images..." << std::endl);

  typedef typename itk::Statistics::ImageToHistogramFilter<ImageType>     HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType            InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType                 HistogramSizeType;

  HistogramSizeType histogramSize( 1 );
  histogramSize[0] = 256;

  typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
  autoMinMaxInputObject->Set( true );

  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput( inputImage );
  histogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
  histogramFilter->SetHistogramSize( histogramSize );
  histogramFilter->SetMarginalScale( 10.0 );
  histogramFilter->Update();

  float lowerValue = histogramFilter->GetOutput()->Quantile( 0, 0 );
  float upperValue = histogramFilter->GetOutput()->Quantile( 0, 1 );

  typedef typename itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;
  typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
  windowingFilter->SetInput( inputImage );
  windowingFilter->SetWindowMinimum( lowerValue );
  windowingFilter->SetWindowMaximum( upperValue );
  windowingFilter->SetOutputMinimum( 0 );
  windowingFilter->SetOutputMaximum( 1 );
  windowingFilter->Update();

  typename ImageType::Pointer outputImage = nullptr;
  outputImage = windowingFilter->GetOutput();
  outputImage->Update();
  outputImage->DisconnectPipeline();

  return outputImage;
}

// debug output for map of vector structure
template <class TMap>
void
PrintMapOfImageVectors(const TMap &map)
{
  muLogMacro(<< "Map size: " << map.size() << std::endl);
  for(typename TMap::const_iterator mapIt = map.begin();
      mapIt != map.end(); ++mapIt)
    {
    muLogMacro(<< "  " << mapIt->first << "(" << mapIt->second.size() <<"):" << std::endl);
    for(unsigned i = 0; i < mapIt->second.size(); ++i)
      {
      muLogMacro( << "    " << mapIt->second[i].GetPointer()
                 << mapIt->second[i]->GetLargestPossibleRegion()
                 << " " << mapIt->second[i]->GetBufferedRegion()
                 << std::endl );
      }
    }
}

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
