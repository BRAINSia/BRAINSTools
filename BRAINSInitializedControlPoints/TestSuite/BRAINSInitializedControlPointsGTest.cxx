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
/// \file BRAINSInitializedControlPointsGTest.cxx
/// \brief Invalid command-line arguments produce EXIT_FAILURE rather than a crash.

#include <gtest/gtest.h>
#include <cstdlib>
#include <string>
#include <vector>

extern "C" int
ModuleEntryPoint(int, char *[]);

TEST(BRAINSInitializedControlPoints, WrongPermuteOrderSizeReturnsFailure)
{
  std::vector<std::string> args{ "BRAINSInitializedControlPoints",
                                 "--inputVolume",
                                 "unused.nii.gz",
                                 "--splineGridSize",
                                 "3,3,3",
                                 "--outputLandmarksFile",
                                 "unused.fcsv",
                                 "--permuteOrder",
                                 "0,1" };
  std::vector<char *>      argv;
  for (auto & a : args)
  {
    argv.push_back(a.data());
  }
  argv.push_back(nullptr);
  EXPECT_EQ(ModuleEntryPoint(static_cast<int>(args.size()), argv.data()), EXIT_FAILURE);
}
