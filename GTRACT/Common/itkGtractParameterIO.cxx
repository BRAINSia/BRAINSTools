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

#ifndef __itkGtractParameterIO_cxx
#define __itkGtractParameterIO_cxx

#include "itkGtractParameterIO.h"

namespace itk
{
GtractParameterIO
::GtractParameterIO() :
  m_Directions(),
  m_Bvalue(0.0),
  m_NumberOfDirections(0),
  m_FileName("")
{
}

void GtractParameterIO::Update()
{
  std::cout << std::endl;
  std::cout << "Loading parameter file...." << std::endl;
  std::cout << m_FileName << std::endl;

  std::ifstream fin( m_FileName.c_str() );

  if( fin )
    {
    // Read in the Nominal b values
    fin >> m_Bvalue;

    // Read in the direction profile
    fin >> m_NumberOfDirections;
      {
      TMatrix temp(m_NumberOfDirections, 6);
      temp.fill(0.0);
      m_Directions = temp;
      }
    for( int i = 0; i < m_NumberOfDirections; i++ )
      {
      fin >> m_Directions[i][0] >> m_Directions[i][1] >> m_Directions[i][2];
      }

    fin.close();

    std::cout << "Done!" << std::endl;
    }
}
} // end namespace itk
#endif
