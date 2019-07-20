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
// /  BRAINSResampleCommonLibWin32Header - manage Windows system differences
// /
// / The BRAINSResampleCommonLibWin32Header captures some system differences
// between Unix
// / and Windows operating systems.

#ifndef __BRAINSResampleCommonLibWin32Header_h
#define __BRAINSResampleCommonLibWin32Header_h

#include <BRAINSResampleCommonLib.h>

#if defined( WIN32 ) && !defined( BRAINSResampleCommonLib_STATIC )
#  if defined( BRAINSResampleCommonLib_EXPORTS )
#    define BRAINSResampleCommonLib_EXPORT __declspec( dllexport )
#  else
#    define BRAINSResampleCommonLib_EXPORT __declspec( dllimport )
#  endif
#else
#  define BRAINSResampleCommonLib_EXPORT
#endif

#endif
