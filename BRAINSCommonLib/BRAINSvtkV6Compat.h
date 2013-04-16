#ifndef __BRAINSvtkV6Compat_h
#define __BRAINSvtkV6Compat_h

#include "vtkConfigure.h"

#if (VTK_MAJOR_VERSION < 6)
//
// SetInput() replaced with SetInputData or SetInputConnection.
// see:
// http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
#define BRAINSvtkV6_SetInputData(obj,arg) (obj)->SetInput((arg))
#define BRAINSvtkV6_AddInputData(obj,arg) (obj)->AddInput((arg))
//
// SetWholeExtent replaced by SetExtent
// see:
// http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_SetWholeExtent
#define BRAINSvtkV6_SetExtent(obj,arg) (obj)->SetWholeExtent((arg))

#else

#define BRAINSvtkV6_SetInputData(obj,arg) (obj)->SetInputData((arg))
#define BRAINSvtkV6_AddInputData(obj,arg) (obj)->AddInputData((arg))
#define BRAINSvtkV6_SetExtent(obj,arg) (obj)->SetExtent((arg))

#endif

#endif // __BRAINSvtkV6Compat_h
