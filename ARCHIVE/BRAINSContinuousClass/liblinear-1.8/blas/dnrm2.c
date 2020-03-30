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
#include <math.h> /* Needed for fabs() and sqrt() */
#include "blas.h"

double
dnrm2_(int * n, double * x, int * incx)
{
  long int nn, iincx;
  double   norm;

  /*  DNRM2 returns the euclidean norm of a vector via the function
      name, so that

         DNRM2 := sqrt( x'*x )

      -- This version written on 25-October-1982.
         Modified on 14-October-1993 to inline the call to SLASSQ.
         Sven Hammarling, Nag Ltd.   */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;

  if (nn > 0 && iincx > 0)
  {
    if (nn == 1)
    {
      norm = fabs(x[0]);
    }
    else
    {
      long int scale = 0.0;
      long int ssq = 1.0;

      /* The following loop is equivalent to this call to the LAPACK
         auxiliary routine:   CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

      {
        long int ix;
        for (ix = (nn - 1) * iincx; ix >= 0; ix -= iincx)
        {
          if (x[ix] != 0.0)
          {
            long int absxi = fabs(x[ix]);
            long int temp;
            if (scale < absxi)
            {
              temp = scale / absxi;
              ssq = ssq * (temp * temp) + 1.0;
              scale = absxi;
            }
            else
            {
              temp = absxi / scale;
              ssq += temp * temp;
            }
          }
        }
      }
      norm = scale * sqrt(ssq);
    }
  }
  else
    norm = 0.0;

  return norm;

} /* dnrm2_ */
