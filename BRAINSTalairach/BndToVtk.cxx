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
#include "BndToVtk.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "itkMacro.h" //Needed for ITK_NULLPTR

#define PR(x) std::cout << #x " = " << x << "\n"; // print macro

vtkStandardNewMacro(BndToVtk);

void BndToVtk::ProcessBND(std::string bndFile)
{
  /* SLA, IRP, AC and PC points are defined in the .bnd file */

  /* Read in .bnd file and move to relevant section at file's end */
  std::ifstream fin( bndFile.c_str() );
  std::string   line;

  getline(fin, line);
  while( line.compare("IPL_HEADER_END") != 0 )
    {
    getline(fin, line);
    }

  /* Data needed is one line past "IPL_HEADER_END" */
  getline(fin, line);

  /* Requested information is 4 {x,y,z} coordinates given in the
   * form of three whitespace-delimited tokens per line; this
   * vector will store all 12 */
  std::vector<std::string> tokens;
  /* Add values to token list */
  for( int i = 0; i < 4; i++ )
    {
    std::string       buf;
    std::stringstream ss(line);
    while( ss >> buf )
      {
      tokens.push_back(buf);
      }

    getline(fin, line);
    }
  /* fill in AC point values */
  for( int i = 0; i < 3; i++ )
    {
    const std::string fillVals = tokens[i];
    AC[i] = strtod(fillVals.c_str(), ITK_NULLPTR);
    }
  /* fill in PC point values */
  for( int i = 3; i < 6; i++ )
    {
    const std::string fillVals = tokens[i];
    PC[i - 3] = strtod(fillVals.c_str(), ITK_NULLPTR);
    }
  /* fill in SLA point values */
  for( int i = 6; i < 9; i++ )
    {
    const std::string fillVals = tokens[i];
    SLA[i - 6] = strtod(fillVals.c_str(), ITK_NULLPTR);
    }
  /* fill in IRP point values */
  for( int i = 9; i < 12; i++ )
    {
    const std::string fillVals = tokens[i];
    IRP[i - 9] = strtod(fillVals.c_str(), ITK_NULLPTR);
    }

  /* Y and Z are inverted in the BRAINS representation vis a vis VTK */

  inverseAC[0] = AC[0];
  inverseAC[1] = AC[2];
  inverseAC[2] = AC[1];

  inversePC[0] = PC[0];
  inversePC[1] = PC[2];
  inversePC[2] = PC[1];

  inverseIRP[0] = IRP[0];
  inverseIRP[1] = IRP[2];
  inverseIRP[2] = IRP[1];

  inverseSLA[0] = SLA[0];
  inverseSLA[1] = SLA[2];
  inverseSLA[2] = SLA[1];
}

double * BndToVtk::GetAC()
{
  return AC;
}

double * BndToVtk::GetPC()
{
  return PC;
}

double * BndToVtk::GetIRP()
{
  return IRP;
}

double * BndToVtk::GetSLA()
{
  return SLA;
}

double * BndToVtk::GetInverseAC()
{
  return inverseAC;
}

double * BndToVtk::GetInversePC()
{
  return inversePC;
}

double * BndToVtk::GetInverseSLA()
{
  return inverseSLA;
}

double * BndToVtk::GetInverseIRP()
{
  return inverseIRP;
}

void BndToVtk::SetAC(double pnt[3])
{
  AC[0] = pnt[0];
  AC[1] = pnt[1];
  AC[2] = pnt[2];
}

void BndToVtk::SetPC(double pnt[3])
{
  PC[0] = pnt[0];
  PC[1] = pnt[1];
  PC[2] = pnt[2];
}

void BndToVtk::SetIRP(double pnt[3])
{
  IRP[0] = pnt[0];
  IRP[1] = pnt[1];
  IRP[2] = pnt[2];
}

void BndToVtk::SetSLA(double pnt[3])
{
  SLA[0] = pnt[0];
  SLA[1] = pnt[1];
  SLA[2] = pnt[2];
}

void BndToVtk::SetInverseAC(double pnt[3])
{
  inverseAC[0] = pnt[0];
  inverseAC[1] = pnt[1];
  inverseAC[2] = pnt[2];
}

void BndToVtk::SetInversePC(double pnt[3])
{
  inversePC[0] = pnt[0];
  inversePC[1] = pnt[1];
  inversePC[2] = pnt[2];
}

void BndToVtk::SetInverseSLA(double pnt[3])
{
  inverseSLA[0] = pnt[0];
  inverseSLA[1] = pnt[1];
  inverseSLA[2] = pnt[2];
}

void BndToVtk::SetInverseIRP(double pnt[3])
{
  inverseIRP[0] = pnt[0];
  inverseIRP[1] = pnt[1];
  inverseIRP[2] = pnt[2];
}
