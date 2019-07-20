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
/*=========================================================================

Program:   Hierarchical Attribute Matching Mechanism for Elastic Registration
Module:    $RCSfile: itkHammerAttributeVectorBase.h,v $
Language:  C++
Date:      $Date: 2009/01/13 20:19:20 $
Version:   $Revision: 1.4 $

Copyright (c)

This program is developed under NIH NCBC collaboration grant
R01 EB006733, "Development and Dissemination of Robust Brain MRI
Measurement Tools".

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHammerAttributeVectorBase_h
#define __itkHammerAttributeVectorBase_h

#include "itkFixedArray.h"
#include "itkNumericTraitsFixedArrayPixel.h"

namespace itk
{
/** \class HammerAttributeVectorBase
 * \brief Abstract base class for attribute vector used by Hammer
 *
 */
template < typename TValueType, unsigned int VLength = 3 >
class HammerAttributeVectorBase : public FixedArray< TValueType, VLength >
{
public:
  using Self = HammerAttributeVectorBase;
  using SuperClass = typename itk::FixedArray< TValueType, VLength >;
  using VectorType = FixedArray< TValueType, VLength >;

  /** Length constant */
  unsigned int
  GetLength() const
  {
    return SuperClass::GetLength();
  }

  /** Dimension constant */
  unsigned int
  GetDimension() const
  {
    return SuperClass::GetDimension();
  }

  /** define interface for computing similarity/difference between two
   * attribute vectors */
  virtual double
  ComputeSimilarity( const VectorType & vec2 ) const
  {
    double d = 1;

    for ( unsigned int k = 0; k < this->Size(); k++ )
    {
      const double c = static_cast< double >( this->GetElement( k ) ) - static_cast< double >( vec2[k] );
      d *= exp( -c * c / 2 );
    }
    return d;
  }

  virtual double
  ComputeDifference( const VectorType & vec2 ) const
  {
    double d = 0;

    for ( unsigned k = 0; k < this->Size(); k++ )
    {
      const double c = static_cast< double >( this->GetElement( k ) ) - static_cast< double >( vec2[k] );
      d += c * c;
    }
    return d;
  }

  virtual bool
  IsQualifiedDrivingVoxel( std::vector< float > & /* qualifier */ )
  {
    return false;
  }

  HammerAttributeVectorBase() {}

  virtual ~HammerAttributeVectorBase() {}

  void
  PrintSelf( std::ostream & /*NOT IMPLEMENTED os*/, Indent /*NOT IMPLEMENTED indent*/ ) const
  {}
};
} // end namespace itk

#endif
