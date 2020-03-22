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
GtractParameterIO ::GtractParameterIO()
  : m_Directions()
  , m_FileName("")
{}

void
GtractParameterIO::Update()
{
  std::cout << std::endl;
  std::cout << "Loading parameter file...." << std::endl;
  std::cout << m_FileName << std::endl;

  std::ifstream fin(m_FileName.c_str());

  if (fin)
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
    for (int i = 0; i < m_NumberOfDirections; i++)
    {
      fin >> m_Directions[i][0] >> m_Directions[i][1] >> m_Directions[i][2];
    }

    fin.close();

    std::cout << "Done!" << std::endl;
  }
}
} // end namespace itk
#endif
