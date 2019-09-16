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
#ifndef __BndToVtk_h
#define __BndToVtk_h

#include "vtkDataObjectAlgorithm.h"
#include "vtkObjectFactory.h"
#include <string>

class BndToVtk : public vtkDataObjectAlgorithm
{
public:
  static BndToVtk *
  New();

  /* Proces the BND file to obtain AC, PC, IRP and SLA points */
  void
  ProcessBND(std::string bndFile);

  /* Get the AC, PC, IRP and SLA points */
  double *
  GetAC();

  double *
  GetPC();

  double *
  GetIRP();

  double *
  GetSLA();

  /* Get the inverted AC, PC, IRP and SLA points */
  double *
  GetInverseAC();

  double *
  GetInversePC();

  double *
  GetInverseSLA();

  double *
  GetInverseIRP();

  /* Set the inverted AC, PC, IRP and SLA points */
private:
  /* AC, PC, IRP and SLA points */
  double AC[3];
  double PC[3];
  double IRP[3];
  double SLA[3];

  /* inverted AC, PC, IRP and SLA points */
  double inverseAC[3];
  double inversePC[3];
  double inverseSLA[3];
  double inverseIRP[3];

  /* Set the AC, PC, IRP and SLA points */
  void
  SetAC(double pnt[3]);
  void
  SetPC(double pnt[3]);
  void
  SetIRP(double pnt[3]);
  void
  SetSLA(double pnt[3]);

  /* Set the inverted AC, PC, IRP and SLA points */
  void
  SetInverseAC(double pnt[3]);
  void
  SetInversePC(double pnt[3]);
  void
  SetInverseSLA(double pnt[3]);
  void
  SetInverseIRP(double pnt[3]);
};

#endif
