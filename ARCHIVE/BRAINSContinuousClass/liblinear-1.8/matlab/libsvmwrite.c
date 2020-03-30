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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"

#if MX_API_VER < 0x07030000
using mwIndex = int;
#endif

void
exit_with_help()
{
  mexPrintf("Usage: libsvmwrite('filename', label_vector, instance_matrix);\n");
}

void
libsvmwrite(const char * filename, const mxArray * label_vec, const mxArray * instance_mat)
{
  FILE *    fp = fopen(filename, "w");
  int       i, k, low, high, l;
  mwIndex * ir, *jc;
  int       label_vector_row_num;
  double *  samples, *labels;
  mxArray * instance_mat_col; // instance sparse matrix in column format

  if (fp == NULL)
  {
    mexPrintf("can't open output file %s\n", filename);
    return;
  }

  // transpose instance matrix
  {
    mxArray *prhs[1], *plhs[1];
    prhs[0] = mxDuplicateArray(instance_mat);
    if (mexCallMATLAB(1, plhs, 1, prhs, "transpose"))
    {
      mexPrintf("Error: cannot transpose instance matrix\n");
      fclose(fp);
      return;
    }
    instance_mat_col = plhs[0];
    mxDestroyArray(prhs[0]);
  }

  // the number of instance
  l = (int)mxGetN(instance_mat_col);
  label_vector_row_num = (int)mxGetM(label_vec);

  if (label_vector_row_num != l)
  {
    mexPrintf("Length of label vector does not match # of instances.\n");
    return;
  }

  // each column is one instance
  labels = mxGetPr(label_vec);
  samples = mxGetPr(instance_mat_col);
  ir = mxGetIr(instance_mat_col);
  jc = mxGetJc(instance_mat_col);

  for (i = 0; i < l; i++)
  {
    fprintf(fp, "%g", labels[i]);

    low = (int)jc[i], high = (int)jc[i + 1];
    for (k = low; k < high; k++)
      fprintf(fp, " %ld:%g", ir[k] + 1, samples[k]);

    fprintf(fp, "\n");
  }

  fclose(fp);
  return;
}
