#ifndef __BndToVtk_h
#define __BndToVtk_h

#include "vtkDataObjectAlgorithm.h"
#include "vtkObjectFactory.h"
#include <string>

class BndToVtk : public vtkDataObjectAlgorithm
{
public:
  static BndToVtk * New();

  /* Proces the BND file to obtain AC, PC, IRP and SLA points */
  void ProcessBND(std::string bndFile);

  /* Get the AC, PC, IRP and SLA points */
  double * GetAC();

  double * GetPC();

  double * GetIRP();

  double * GetSLA();

  /* Get the inverted AC, PC, IRP and SLA points */
  double * GetInverseAC();

  double * GetInversePC();

  double * GetInverseSLA();

  double * GetInverseIRP();

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
  void SetAC(double pnt[3]);
  void SetPC(double pnt[3]);
  void SetIRP(double pnt[3]);
  void SetSLA(double pnt[3]);

  /* Set the inverted AC, PC, IRP and SLA points */
  void SetInverseAC(double pnt[3]);
  void SetInversePC(double pnt[3]);
  void SetInverseSLA(double pnt[3]);
  void SetInverseIRP(double pnt[3]);
};

#endif
