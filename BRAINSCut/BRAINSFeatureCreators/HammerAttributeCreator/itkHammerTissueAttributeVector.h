/*=========================================================================

Program:   Hierarchical Attribute Matching Mechanism for Elastic Registration
Module:    $RCSfile: itkHammerTissueAttributeVector.h,v $
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
#ifndef __itkHammerTissueAttributeVector_h
#define __itkHammerTissueAttributeVector_h

#include "itkHammerAttributeVectorBase.h"

namespace itk
{
/** \class HammerTissueAttributeVector
  * \brief Abstract base class for attribute vector used by Hammer
  *
  */
class ITK_EXPORT HammerTissueAttributeVector :
  public         HammerAttributeVectorBase<unsigned char, 5>
{
public:

  typedef HammerTissueAttributeVector Self;

  /** Length constant */
  unsigned int GetLength() const
  {
    return 5;
  }

  /** Dimension constant */
  unsigned int GetDimension() const
  {
    return 5;
  }

  /** define interface for computing similarity/difference between two
    * attribute vectors */
  virtual double ComputeSimilarity(Self & vec2)
  {
    if( this->GetEdge() != vec2.GetEdge() )
      {
      return 0;
      }
    else
      {
      double sim = 1;
      for( int k = 1; k < 5; k++ )
        {
        sim *= ( 1.0 - fabs( static_cast<double>( this->operator[](k) ) - static_cast<double>( vec2[k] ) ) / 255.0 );
        }
      return sim;
      }
  }

  virtual double ComputeDifference(Self & vec2)
  {
    double diff = ComputeSimilarity(vec2);

    if( diff == 0 )
      {
      return -1;
      }
    else
      {
      return 1 / diff;
      }
  }

  unsigned char GetEdge()
  {
    return this->operator[](0);
  }

  unsigned char GetTissueType()
  {
    return this->operator[](1);
  }

  unsigned char GetGeometry()
  {
    return this->operator[](2);
  }

  unsigned char GetVentricleVolume()
  {
    return this->operator[](3);
  }

  unsigned char GetCSFBackground()
  {
    return this->operator[](4);
  }

  virtual bool IsQualifiedDrivingVoxel(std::vector<float> & qualifier)
  {
    if( this->GetEdge() == 0 )
      {
      return false;
      }

    if( this->GetGeometry() <= qualifier[0]
        && ( this->GetGeometry() >= qualifier[1] || this->GetGeometry() <= qualifier[2] )
        && this->GetVentricleVolume() <= qualifier[3]
        && this->GetCSFBackground() <= qualifier[5] )
      {
      return false;
      }

    if( this->GetVentricleVolume() > qualifier[4] )
      {
      return false;
      }

    return true;
  }

  virtual bool IsQualifiedDrivingVoxel_GR(std::vector<float> & qualifier)
  {
    if( this->GetEdge() == 0 )
      {
      return false;
      }

    if( this->GetGeometry() > qualifier[0]
        || ( this->GetGeometry() < qualifier[1] && this->GetGeometry() > qualifier[2] )
        || this->GetVentricleVolume() > qualifier[3]
        || this->GetCSFBackground() > qualifier[5] )
      {
      return true;
      }
    else
      {
      return false;
      }
  }

  HammerTissueAttributeVector()
  {
  }

  ~HammerTissueAttributeVector()
  {
  }

  void PrintSelf(std::ostream & /* os */, Indent /* indent */) const
  {
  }
};
}   // end namespace itk

#endif
