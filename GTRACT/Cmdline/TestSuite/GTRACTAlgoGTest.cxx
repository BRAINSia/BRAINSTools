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
/// \file GTRACTAlgoGTest.cxx
/// \brief My_lsf rejects mismatched input sizes with a catchable exception.

#include <gtest/gtest.h>

#include "itkMacro.h"
#include "algo.h"

TEST(GTRACTAlgo, LeastSquaresFitSizeMismatchThrows)
{
  TVector x(3, 1.0f);
  TVector y(4, 1.0f);
  EXPECT_THROW(My_lsf(x, y), itk::ExceptionObject);
}
