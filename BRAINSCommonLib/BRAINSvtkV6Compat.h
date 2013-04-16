#ifndef __BRAINSvtkV6Compat_h
#define __BRAINSvtkV6Compat_h

#include "vtkConfigure.h"
#include "vtkSmartPointer.h"
#if (VTK_MAJOR_VERSION < 6)
//
// SetInput() replaced with SetInputData or SetInputConnection.
// see:
// http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
template <class TObject, class TParm>
void
BRAINSvtkV6_SetInputData(TObject *obj, TParm *arg)
{
  obj->SetInput(arg);
}

template <class TObject, class TParm>
void
BRAINSvtkV6_AddInputData(TObject *obj, TParm *arg)
{
  obj->AddInput(arg);
}

//
// SetWholeExtent replaced by SetExtent
// see:
// http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_SetWholeExtent
template <class TObject, class TParm>
void
BRAINSvtkV6_SetExtent(TObject *obj, TParm *arg)
{
  obj->SetWholeExtent(arg);
}

#else

/* #define BRAINSvtkV6_SetInputData(obj,arg) (obj)->SetInputData((arg)) */
template <class TObject, class TParm>
void
BRAINSvtkV6_SetInputData(TObject *obj, TParm *arg)
{
  obj->SetInputData(arg);
}

template <class TObject, class TParm>
void
BRAINSvtkV6_AddInputData(TObject *obj, TParm *arg)
{
  obj->AddInputData(arg);
}

template <class TObject, class TParm>
void
BRAINSvtkV6_SetExtent(TObject *obj, TParm *arg)
{
  obj->SetExtent(arg);
}

#endif

// Additional template methods to
// help compiler with different combinations of smart pointers
// and regular pointers
#define BRAINSvtkV6_TemplateHelpers(FuncName) \
template <class TObject, class TParm> \
void \
BRAINSvtkV6_##FuncName(vtkSmartPointer<TObject> &obj, vtkSmartPointer<TParm> &arg) \
{  BRAINSvtkV6_##FuncName(obj.GetPointer(),arg.GetPointer()); } \
template <class TObject, class TParm> \
void \
BRAINSvtkV6_##FuncName(vtkSmartPointer<TObject> &obj, TParm *arg) \
{ BRAINSvtkV6_##FuncName(obj.GetPointer(),arg); } \
template <class TObject, class TParm> \
void \
BRAINSvtkV6_##FuncName(TObject *obj, vtkSmartPointer<TParm> &arg) \
{  BRAINSvtkV6_##FuncName(obj,arg.GetPointer()); }

BRAINSvtkV6_TemplateHelpers(SetInputData)
BRAINSvtkV6_TemplateHelpers(AddInputData)
BRAINSvtkV6_TemplateHelpers(SetExtent)


#endif // __BRAINSvtkV6Compat_h
