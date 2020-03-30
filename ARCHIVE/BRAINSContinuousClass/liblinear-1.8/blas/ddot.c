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

double
ddot_(int * n, double * sx, int * incx, double * sy, int * incy)
{
  long int nn, iincx, iincy;
  double   stemp;
  long int ix, iy;

  /* forms the dot product of two vectors.
     uses unrolled loops for increments equal to one.
     jack dongarra, linpack, 3/11/78.
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;
  iincy = *incy;

  stemp = 0.0;
  if (nn > 0)
  {
    long int i;
    if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
    {
      long int m = nn - 4;
      for (i = 0; i < m; i += 5)
        stemp +=
          sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * sy[i + 2] + sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];

      for (; i < nn; i++) /* clean-up loop */
        stemp += sx[i] * sy[i];
    }
    else /* code for unequal increments or equal increments not equal to 1 */
    {
      ix = 0;
      iy = 0;
      if (iincx < 0)
        ix = (1 - nn) * iincx;
      if (iincy < 0)
        iy = (1 - nn) * iincy;
      for (i = 0; i < nn; i++)
      {
        stemp += sx[ix] * sy[iy];
        ix += iincx;
        iy += iincy;
      }
    }
  }

  return stemp;
} /* ddot_ */
