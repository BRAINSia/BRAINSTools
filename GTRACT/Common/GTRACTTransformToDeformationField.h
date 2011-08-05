/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// ##############HACK:  This sould use the version out of ITK
// Review/itkTransformToDeformationFieldSource.h
//
//

#ifndef TransformToDeformationField_h
#define TransformToDeformationField_h
#include "itkSpatialOrientationAdapter.h"
// #include "itkIO.h"
// #include "itkUtil.h"
// #include "CrossOverAffineSystem.h"

/**
 * Go from any subclass of Transform, to the corresponding deformation field
 */
template <typename DeformationFieldType, typename TransformType>
typename DeformationFieldType::Pointer TransformToDeformationField( typename DeformationFieldType::SizeType size,
                                                                    double *spacing,
                                                                    double *origin,
                                                                    typename DeformationFieldType::DirectionType Dir,
                                                                    typename TransformType::Pointer & xfrm)
{
  typedef typename DeformationFieldType::PixelType DeformationPixelType;

  typename DeformationFieldType::Pointer deformation = DeformationFieldType::New();
  typename DeformationFieldType::RegionType region;

  region.SetSize(size);

  deformation->SetSpacing( spacing );
  deformation->SetRegions( region );
  deformation->SetOrigin(origin);
  deformation->SetDirection( Dir );
  deformation->Allocate();
  for( unsigned int z = 0; z < size[2]; z++ )
    {
    for( unsigned int y = 0; y < size[1]; y++ )
      {
      for( unsigned int x = 0; x < size[0]; x++ )
        {
        itk::Point<double, 3> a, b;
        a[0] = spacing[0] * x + origin[0];
        a[1] = spacing[1] * y + origin[1];
        a[2] = spacing[2] * z + origin[2];

        b = xfrm->TransformPoint(a);
        DeformationPixelType p;
        p[0] = b[0] - a[0];
        p[1] = b[1] - a[1];
        p[2] = b[2] - a[2];
        // std::cout<<"p[0]=" << p[0]<< "p[1]= "<< p[1] << "p[2]=" <<p[2]
        // <<std::endl;

        typename DeformationFieldType::IndexType ind;
        ind[0] = x;
        ind[1] = y;
        ind[2] = z;
        // std::cout<<"i[0]=" << ind[0]<< "i[1]= "<< ind[1] << "i[2]=" <<ind[2]
        // <<std::endl;
        deformation->SetPixel(ind, p);
        }
      }
    }
  return deformation;
}

#endif /* TransformToDeformationField_h*/
