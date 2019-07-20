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
dscal_( int * n, double * sa, double * sx, int * incx )
{
  long int i, nn, iincx;
  double   ssa;

  /* scales a vector by a constant.
     uses unrolled loops for increment equal to 1.
     jack dongarra, linpack, 3/11/78.
     modified 3/93 to return if incx .le. 0.
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;
  ssa = *sa;

  if ( nn > 0 && iincx > 0 )
  {
    if ( iincx == 1 ) /* code for increment equal to 1 */
    {
      long int m = nn - 4;
      for ( i = 0; i < m; i += 5 )
      {
        sx[i] = ssa * sx[i];
        sx[i + 1] = ssa * sx[i + 1];
        sx[i + 2] = ssa * sx[i + 2];
        sx[i + 3] = ssa * sx[i + 3];
        sx[i + 4] = ssa * sx[i + 4];
      }
      for ( ; i < nn; ++i ) /* clean-up loop */
        sx[i] = ssa * sx[i];
    }
    else /* code for increment not equal to 1 */
    {
      long int nincx = nn * iincx;
      for ( i = 0; i < nincx; i += iincx )
        sx[i] = ssa * sx[i];
    }
  }

  return 0;
} /* dscal_ */
