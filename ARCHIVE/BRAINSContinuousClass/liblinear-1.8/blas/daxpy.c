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
#include "blas.h"

int
daxpy_(int * n, double * sa, double * sx, int * incx, double * sy, int * incy)
{
  long int        nn, iincx, iincy;
  register double ssa;

  /* constant times a vector plus a vector.
     uses unrolled loop for increments equal to one.
     jack dongarra, linpack, 3/11/78.
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  ssa = *sa;
  iincx = *incx;
  iincy = *incy;

  if (nn > 0 && ssa != 0.0)
  {
    long int i;
    if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
    {
      long int m = nn - 3;
      for (i = 0; i < m; i += 4)
      {
        sy[i] += ssa * sx[i];
        sy[i + 1] += ssa * sx[i + 1];
        sy[i + 2] += ssa * sx[i + 2];
        sy[i + 3] += ssa * sx[i + 3];
      }
      for (; i < nn; ++i) /* clean-up loop */
        sy[i] += ssa * sx[i];
    }
    else /* code for unequal increments or equal increments not equal to 1 */
    {
      long int ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
      long int iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
      for (i = 0; i < nn; i++)
      {
        sy[iy] += ssa * sx[ix];
        ix += iincx;
        iy += iincy;
      }
    }
  }

  return 0;
} /* daxpy_ */
