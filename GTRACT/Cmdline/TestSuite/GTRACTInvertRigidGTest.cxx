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
/// \file GTRACTInvertRigidGTest.cxx
/// \brief gtractInvertRigidTransform writes the true inverse of a rigid transform.

#include <gtest/gtest.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "itkMacro.h"
#include "itkVersorRigid3DTransform.h"
#include "GenericTransformImage.h"

extern "C" int
ModuleEntryPoint(int, char *[]);

namespace
{
int
RunInvert(const std::string & input, const std::string & output)
{
  std::vector<std::string> args{ "gtractInvertRigidTransform", "--inputTransform", input, "--outputTransform", output };
  std::vector<char *>      argv;
  for (auto & a : args)
  {
    argv.push_back(a.data());
  }
  argv.push_back(nullptr);
  return ModuleEntryPoint(static_cast<int>(args.size()), argv.data());
}
} // namespace

TEST(GTRACTInvertRigidTransform, InverseComposesToIdentity)
{
  using RigidTransformType = itk::VersorRigid3DTransform<double>;
  auto                               forward = RigidTransformType::New();
  RigidTransformType::ParametersType params(6);
  params[0] = 0.1;
  params[1] = 0.2;
  params[2] = 0.05;
  params[3] = 10.0;
  params[4] = -5.0;
  params[5] = 3.0;
  forward->SetParameters(params);
  RigidTransformType::FixedParametersType center(3);
  center[0] = 1.0;
  center[1] = 2.0;
  center[2] = 3.0;
  forward->SetFixedParameters(center);

  const std::string tmp = ::testing::TempDir();
  const std::string forwardFile = tmp + "/gtractInvert_forward.h5";
  const std::string inverseFile = tmp + "/gtractInvert_inverse.h5";
  itk::WriteTransformToDisk<double>(forward.GetPointer(), forwardFile);

  ASSERT_EQ(RunInvert(forwardFile, inverseFile), EXIT_SUCCESS);

  const auto         inverseGeneric = itk::ReadTransformFromDisk(inverseFile);
  const auto * const inverse = dynamic_cast<const RigidTransformType *>(inverseGeneric.GetPointer());
  ASSERT_NE(inverse, nullptr);

  RigidTransformType::InputPointType point;
  point[0] = 4.0;
  point[1] = -7.0;
  point[2] = 12.5;
  const auto roundTrip = inverse->TransformPoint(forward->TransformPoint(point));
  for (unsigned int i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(roundTrip[i], point[i], 1e-9) << "component " << i;
  }

  std::remove(forwardFile.c_str());
  std::remove(inverseFile.c_str());
}
