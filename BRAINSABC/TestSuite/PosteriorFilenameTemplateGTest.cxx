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
/// \file PosteriorFilenameTemplateGTest.cxx
/// \brief --posteriorTemplate is a filename pattern, never a printf format string.

#include <gtest/gtest.h>
#include <string>

#include "itkMacro.h"
#include "PosteriorFilenameTemplate.h"

TEST(PosteriorFilenameTemplate, SubstitutesPriorName)
{
  EXPECT_EQ(ExpandPosteriorTemplate("POST_%s.nii.gz", "GM"), "POST_GM.nii.gz");
}

TEST(PosteriorFilenameTemplate, DoubledPercentIsLiteral)
{
  EXPECT_EQ(ExpandPosteriorTemplate("100%%_%s.nrrd", "WM"), "100%_WM.nrrd");
}

TEST(PosteriorFilenameTemplate, LongNameIsNotTruncatedOrOverrun)
{
  const std::string name(10000, 'x');
  const std::string result = ExpandPosteriorTemplate("P_%s.nii.gz", name);
  EXPECT_EQ(result.size(), name.size() + std::string("P_.nii.gz").size());
}

TEST(PosteriorFilenameTemplate, SecondStringConversionThrows)
{
  EXPECT_THROW(ExpandPosteriorTemplate("%s_%s.nii.gz", "GM"), itk::ExceptionObject);
}

TEST(PosteriorFilenameTemplate, WriteBackConversionThrows)
{
  EXPECT_THROW(ExpandPosteriorTemplate("POST_%n_%s.nii.gz", "GM"), itk::ExceptionObject);
}

TEST(PosteriorFilenameTemplate, NumericConversionThrows)
{
  EXPECT_THROW(ExpandPosteriorTemplate("POST_%x_%s.nii.gz", "GM"), itk::ExceptionObject);
}

TEST(PosteriorFilenameTemplate, TrailingPercentThrows)
{
  EXPECT_THROW(ExpandPosteriorTemplate("POST_%s%", "GM"), itk::ExceptionObject);
}
