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
/// \file OutputVolumeTemplateGTest.cxx
/// \brief --outputVolumes is a filename pattern, never a printf format string.

#include <gtest/gtest.h>
#include <string>

#include "itkMacro.h"
#include "OutputVolumeTemplate.h"

TEST(OutputVolumeTemplate, DefaultPatternSubstitutesTypeAndIndex)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("%s_corrected_%d.nii.gz", "T1", 0), "T1_corrected_0.nii.gz");
  EXPECT_EQ(ExpandOutputVolumeTemplate("%s_corrected_%d.nii.gz", "T2", 12), "T2_corrected_12.nii.gz");
}

TEST(OutputVolumeTemplate, ZeroPaddedIndexIsHonored)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("%s_%03d.nrrd", "PD", 7), "PD_007.nrrd");
}

TEST(OutputVolumeTemplate, WideZeroPaddingIsNotTruncated)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("%s_%012d.nii", "T1", 7), "T1_000000000007.nii");
  EXPECT_EQ(ExpandOutputVolumeTemplate("%s_%099d", "T1", 5).size(), std::string("T1_").size() + 99);
}

TEST(OutputVolumeTemplate, DoubledPercentIsLiteral)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("100%%_%s_%d.nii.gz", "FL", 3), "100%_FL_3.nii.gz");
}

TEST(OutputVolumeTemplate, PatternWithoutConversionsIsUnchanged)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("fixed_name.nii.gz", "T1", 4), "fixed_name.nii.gz");
}

TEST(OutputVolumeTemplate, IndexMayPrecedeTypeInThePattern)
{
  EXPECT_EQ(ExpandOutputVolumeTemplate("%d_%s.nii.gz", "T1", 2), "2_T1.nii.gz");
}

TEST(OutputVolumeTemplate, LongTypeNameIsNotTruncatedOrOverrun)
{
  const std::string type(10000, 'x');
  const std::string result = ExpandOutputVolumeTemplate("%s_%d.nii.gz", type, 1);
  EXPECT_EQ(result.size(), type.size() + std::string("_1.nii.gz").size());
}

TEST(OutputVolumeTemplate, SecondStringConversionThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%s_%d.nii.gz", "T1", 0), itk::ExceptionObject);
}

TEST(OutputVolumeTemplate, SecondIndexConversionThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%d_%d.nii.gz", "T1", 0), itk::ExceptionObject);
}

TEST(OutputVolumeTemplate, WriteBackConversionThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%n_%d.nii.gz", "T1", 0), itk::ExceptionObject);
}

TEST(OutputVolumeTemplate, UnsupportedConversionThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%x.nii.gz", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%*d.nii.gz", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%10s_%d.nii.gz", "T1", 0), itk::ExceptionObject);
}

TEST(OutputVolumeTemplate, TrailingPercentThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%d%", "T1", 0), itk::ExceptionObject);
}

TEST(OutputVolumeTemplate, MalformedWidthThrows)
{
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%0d.nii", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%03x.nii", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%03", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%003d.nii", "T1", 0), itk::ExceptionObject);
  EXPECT_THROW(ExpandOutputVolumeTemplate("%s_%3d.nii", "T1", 0), itk::ExceptionObject);
}
