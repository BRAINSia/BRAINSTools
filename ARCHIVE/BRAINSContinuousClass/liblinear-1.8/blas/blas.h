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
/* blas.h  --  C header file for BLAS                         Ver 1.0 */
/* Jesse Bennett                                       March 23, 2000 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

  - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef BLAS_INCLUDE
#define BLAS_INCLUDE

/* Data types specific to BLAS implementation */
typedef struct { float r, i; }  fcomplex;
typedef struct { double r, i; } dcomplex;
typedef int                     blasbool;

#include "blasp.h"    /* Prototypes for all BLAS functions */

#define FALSE 0
#define TRUE  1

/* Macro functions */
#define MIN(a, b) ( (a) <= (b) ? (a) : (b) )
#define MAX(a, b) ( (a) >= (b) ? (a) : (b) )

#endif
