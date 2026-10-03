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
/// \file SiemensBMatrixGTest.cxx
/// \brief ExtractBMatrix reads the six B-matrix values from private tag 0019|1027.
///
/// Built with hardened standard-library checks so an out-of-bounds vector access
/// in SiemensDWIConverter aborts instead of silently reading reserved memory.

#include <gtest/gtest.h>
#include <cstdio>
#include <string>

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcuid.h"
#include "dcmtk/dcmdata/dcvrfd.h"

#include "brainsDCMTKFileReader.h"
#include "SiemensDWIConverter.h"

TEST(SiemensDWIConverter, ExtractBMatrixFromPrivateTag)
{
  const std::string path = std::string(::testing::TempDir()) + "/siemens_bmatrix.dcm";
  const double      expected[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  {
    DcmFileFormat fileFormat;
    DcmDataset *  dataset = fileFormat.getDataset();
    dataset->putAndInsertString(DCM_SOPClassUID, UID_MRImageStorage);
    dataset->putAndInsertString(DCM_SOPInstanceUID, "1.2.826.0.1.3680043.8.498.1");
    dataset->putAndInsertString(DCM_Manufacturer, "SIEMENS");
    dataset->putAndInsertString(DcmTag(0x0019, 0x0010, EVR_LO), "SIEMENS MR HEADER");
    auto * bMatrixElement = new DcmFloatingPointDouble(DcmTag(0x0019, 0x1027, EVR_FD));
    ASSERT_TRUE(bMatrixElement->putFloat64Array(expected, 6).good());
    ASSERT_TRUE(dataset->insert(bMatrixElement).good());
    ASSERT_TRUE(fileFormat.saveFile(path.c_str(), EXS_LittleEndianExplicit).good());
  }

  auto * header = new brains::DCMTKFileReader;
  header->SetFileName(path);
  header->LoadFile();
  SiemensDWIConverter::DCMTKFileVector headers{ header };
  DWIConverter::FileNamesContainer     names{ path };
  SiemensDWIConverter                  converter(headers, names, true, 0.2);

  vnl_matrix_fixed<double, 3, 3> bMatrix;
  bMatrix.fill(0.0);
  ASSERT_TRUE(converter.ExtractBMatrix(nullptr, 0, bMatrix));
  EXPECT_DOUBLE_EQ(bMatrix[0][0], 1.0);
  EXPECT_DOUBLE_EQ(bMatrix[0][1], 2.0);
  EXPECT_DOUBLE_EQ(bMatrix[0][2], 3.0);
  EXPECT_DOUBLE_EQ(bMatrix[1][1], 4.0);
  EXPECT_DOUBLE_EQ(bMatrix[1][2], 5.0);
  EXPECT_DOUBLE_EQ(bMatrix[2][2], 6.0);
  EXPECT_DOUBLE_EQ(bMatrix[1][0], bMatrix[0][1]);
  EXPECT_DOUBLE_EQ(bMatrix[2][0], bMatrix[0][2]);
  EXPECT_DOUBLE_EQ(bMatrix[2][1], bMatrix[1][2]);

  delete header;
  std::remove(path.c_str());
}
