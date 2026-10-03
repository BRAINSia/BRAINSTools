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
/// \file LandmarkLookupGTest.cxx
/// \brief Missing landmarks and malformed model files raise catchable itk::ExceptionObject.

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "itkMacro.h"
#include "itk_hdf5.h"
#include "itk_H5Cpp.h"

#include "LLSModel.h"
#include "landmarksDataSet.h"
#include "landmarksConstellationModelBase.h"

namespace
{
class TestModel : public landmarksConstellationModelBase
{};
} // namespace

TEST(LandmarkLookup, MissingRadiusThrows)
{
  const TestModel model;
  EXPECT_THROW(model.GetRadius("NoSuchLandmark"), itk::ExceptionObject);
}

TEST(LandmarkLookup, MissingHeightThrows)
{
  const TestModel model;
  EXPECT_THROW(model.GetHeight("NoSuchLandmark"), itk::ExceptionObject);
}

TEST(LandmarkLookup, MissingNamedPointThrows)
{
  const landmarksDataSet dataSet;
  EXPECT_THROW(dataSet.GetNamedPoint("NoSuchLandmark"), itk::ExceptionObject);
}

namespace
{
std::string
WriteModelWithCorruptedDataSet(const std::string & fileName, const char * groupName, bool makeVector)
{
  LLSModel model;
  model.SetFileName(fileName);
  LLSModel::LLSMeansType means;
  means["AC"] = std::vector<double>{ 1.0, 2.0, 3.0 };
  model.SetLLSMeans(means);
  LLSModel::LLSMatricesType matrices;
  matrices["AC"] = LLSModel::MatrixType(3, 3, 0.0);
  model.SetLLSMatrices(matrices);
  LLSModel::LLSSearchRadiiType radii;
  radii["AC"] = 5.0;
  model.SetSearchRadii(radii);
  EXPECT_EQ(model.Write(), 0);

  const std::string dataSetName = std::string(groupName) + "/AC";
  H5::H5File        file(fileName, H5F_ACC_RDWR);
  file.unlink(dataSetName);
  if (makeVector)
  {
    const hsize_t dim = 3;
    file.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, &dim));
  }
  else
  {
    const hsize_t dims[2] = { 3, 3 };
    file.createDataSet(dataSetName, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims));
  }
  file.close();
  return fileName;
}
} // namespace

TEST(LLSModelRead, MatrixWhereVectorExpectedThrows)
{
  const std::string fileName = std::string(::testing::TempDir()) + "/LLSModel_badVector.h5";
  WriteModelWithCorruptedDataSet(fileName, "LLSMeans", false);
  LLSModel model;
  model.SetFileName(fileName);
  EXPECT_THROW(model.Read(), itk::ExceptionObject);
  std::remove(fileName.c_str());
}

TEST(LLSModelRead, VectorWhereMatrixExpectedThrows)
{
  const std::string fileName = std::string(::testing::TempDir()) + "/LLSModel_badMatrix.h5";
  WriteModelWithCorruptedDataSet(fileName, "LLSMatrices", true);
  LLSModel model;
  model.SetFileName(fileName);
  EXPECT_THROW(model.Read(), itk::ExceptionObject);
  std::remove(fileName.c_str());
}
