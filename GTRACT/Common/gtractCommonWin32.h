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

#ifndef __gtractCommonWin32_h
#define __gtractCommonWin32_h

#include "gtractConfigure.h"

#if defined( WIN32 ) && !defined( GTRACTSTATIC )
#if defined( GTRACTCommon_EXPORTS )
#define GTRACT_COMMON_EXPORT __declspec( dllexport )
#else
#define GTRACT_COMMON_EXPORT __declspec( dllimport )
#endif
#else
#define GTRACT_COMMON_EXPORT
#endif

#endif   // __gtractCommonWin32_h
