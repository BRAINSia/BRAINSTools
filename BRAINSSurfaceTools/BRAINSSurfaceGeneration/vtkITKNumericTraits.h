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
///  vtkITKNumericTraits - Extra itk::NumericTraits instantiations for VTK
///
/// vtkITKNumericTraits provides extra instantiations for itk::NumericTraits for VTK scalar types.
#ifndef __vtkITKNumericTraits_h
#define __vtkITKNumericTraits_h

#include "vtkITK.h"
#include "vtkSystemIncludes.h"

#include "itkNumericTraits.h"

namespace itk
{
#if defined(VTK_TYPE_USE___INT64)
template <>
class NumericTraits<__int64> : public vcl_numeric_limits<__int64>
{
public:
  typedef __int64          ValueType;
  typedef __int64          PrintType;
  typedef unsigned __int64 AbsType;
  typedef __int64          AccumulateType;
  typedef double           RealType;
  typedef RealType         ScalarRealType;
  typedef float            FloatType;
  static const __int64 VTK_ITK_EXPORT Zero;
  static const __int64 VTK_ITK_EXPORT One;

  static __int64 min()
  {
    return vcl_numeric_limits<__int64>::min();
  }

  static __int64 max()
  {
    return vcl_numeric_limits<__int64>::max();
  }

  static __int64 min( __int64 )
  {
    return vcl_numeric_limits<__int64>::min();
  }

  static __int64 max( __int64 )
  {
    return vcl_numeric_limits<__int64>::max();
  }

  static __int64 NonpositiveMin()
  {
    return min();
  }

  static bool IsPositive(__int64 val)
  {
    return val > Zero;
  }

  static bool IsNonpositive(__int64 val)
  {
    return val <= Zero;
  }

  static bool IsNegative(__int64 val)
  {
    return val < Zero;
  }

  static bool IsNonnegative(__int64 val)
  {
    return val >= Zero;
  }

  static __int64  ZeroValue()
  {
    return Zero;
  }
};

template <>
class NumericTraits<unsigned __int64> : public vcl_numeric_limits<unsigned __int64>
{
public:
  typedef unsigned __int64 ValueType;
  typedef unsigned __int64 PrintType;
  typedef unsigned __int64 AbsType;
  typedef unsigned __int64 AccumulateType;
  typedef double           RealType;
  typedef RealType         ScalarRealType;
  typedef float            FloatType;
  static const unsigned __int64 VTK_ITK_EXPORT Zero;
  static const unsigned __int64 VTK_ITK_EXPORT One;

  static unsigned __int64 min()
  {
    return vcl_numeric_limits<unsigned __int64>::min();
  }

  static unsigned __int64 max()
  {
    return vcl_numeric_limits<unsigned __int64>::max();
  }

  static unsigned __int64 min( unsigned __int64 )
  {
    return vcl_numeric_limits<unsigned __int64>::min();
  }

  static unsigned __int64 max( unsigned __int64 )
  {
    return vcl_numeric_limits<unsigned __int64>::max();
  }

  static unsigned __int64 NonpositiveMin()
  {
    return min();
  }

  static bool IsPositive(unsigned __int64 val)
  {
    return val != Zero;
  }

  static bool IsNonpositive(unsigned __int64 val)
  {
    return val == Zero;
  }

  static bool IsNegative(unsigned __int64)
  {
    return false;
  }

  static bool IsNonnegative(unsigned __int64)
  {
    return true;
  }

  static unsigned __int64  ZeroValue()
  {
    return Zero;
  }
};
#endif

}

#endif /// namespace
