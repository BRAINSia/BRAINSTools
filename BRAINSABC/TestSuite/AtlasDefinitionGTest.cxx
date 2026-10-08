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
/// \file AtlasDefinitionGTest.cxx
/// \brief AtlasDefinition reports bad atlases and missing entries as itk::ExceptionObject.

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>

#include "itkMacro.h"
#include "AtlasDefinition.h"

namespace
{
const char * const kValidAtlas = R"(<Atlas>
 <AtlasImage>
  <type>T1</type>
  <filename>t1.nii.gz</filename>
 </AtlasImage>
 <BrainMask>
  <filename>mask.nii.gz</filename>
 </BrainMask>
 <HeadRegion>
  <filename>head.nii.gz</filename>
 </HeadRegion>
 <Prior>
  <type>WM</type>
  <filename>wm.nii.gz</filename>
  <Weight>1.5</Weight>
  <GaussianClusterCount>2</GaussianClusterCount>
  <UseForBias>1</UseForBias>
  <IsForegroundPrior>1</IsForegroundPrior>
  <LabelCode>7</LabelCode>
  <bounds>
   <type>T1</type>
   <lower>0.25</lower>
   <upper>0.75</upper>
  </bounds>
 </Prior>
</Atlas>
)";

std::string
WriteAtlasFile(const std::string & name, const std::string & contents)
{
  const std::string path = std::string(::testing::TempDir()) + "/" + name;
  std::ofstream     out(path);
  out << contents;
  return path;
}

std::string
Replace(std::string text, const std::string & from, const std::string & to)
{
  text.replace(text.find(from), from.size(), to);
  return text;
}
} // namespace

TEST(AtlasDefinition, ValidAtlasAnswersEveryAccessor)
{
  const std::string path = WriteAtlasFile("atlas_valid.xml", kValidAtlas);
  AtlasDefinition   atlas;
  atlas.InitFromXML(path);
  EXPECT_EQ(atlas.GetPriorFilename("WM"), "wm.nii.gz");
  EXPECT_DOUBLE_EQ(atlas.GetWeight("WM"), 1.5);
  EXPECT_EQ(atlas.GetGaussianClusterCount("WM"), 2);
  EXPECT_EQ(atlas.GetLabelCode("WM"), 7);
  EXPECT_TRUE(atlas.GetUseForBias("WM"));
  EXPECT_TRUE(atlas.GetIsForegroundPrior("WM"));
  EXPECT_DOUBLE_EQ(atlas.GetLow("WM", "T1"), 0.25);
  EXPECT_DOUBLE_EQ(atlas.GetHigh("WM", "T1"), 0.75);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, MissingTissueTypeThrowsFromEveryAccessor)
{
  const std::string path = WriteAtlasFile("atlas_missing_tissue.xml", kValidAtlas);
  AtlasDefinition   atlas;
  atlas.InitFromXML(path);
  EXPECT_THROW(atlas.GetPriorFilename("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetWeight("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetGaussianClusterCount("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetLabelCode("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetUseForBias("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetIsForegroundPrior("CSF"), itk::ExceptionObject);
  EXPECT_THROW(atlas.GetBounds("CSF", "T1"), itk::ExceptionObject);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, MissingModalityBoundsThrows)
{
  const std::string path = WriteAtlasFile("atlas_missing_modality.xml", kValidAtlas);
  AtlasDefinition   atlas;
  atlas.InitFromXML(path);
  EXPECT_THROW(atlas.GetLow("WM", "T2"), itk::ExceptionObject);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, UnreadableFileThrows)
{
  AtlasDefinition atlas;
  EXPECT_THROW(atlas.InitFromXML("/nonexistent/path/atlas.xml"), itk::ExceptionObject);
}

TEST(AtlasDefinition, MalformedXMLThrows)
{
  const std::string path = WriteAtlasFile("atlas_malformed.xml", "<Atlas><Prior></Atlas>");
  AtlasDefinition   atlas;
  EXPECT_THROW(atlas.InitFromXML(path), itk::ExceptionObject);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, NonNumericWeightThrows)
{
  const std::string path = WriteAtlasFile("atlas_bad_weight.xml", Replace(kValidAtlas, "1.5", "heavy"));
  AtlasDefinition   atlas;
  EXPECT_THROW(atlas.InitFromXML(path), itk::ExceptionObject);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, NonNumericLabelCodeThrows)
{
  const std::string path = WriteAtlasFile("atlas_bad_label.xml", Replace(kValidAtlas, ">7<", ">seven<"));
  AtlasDefinition   atlas;
  EXPECT_THROW(atlas.InitFromXML(path), itk::ExceptionObject);
  std::remove(path.c_str());
}

TEST(AtlasDefinition, UnknownElementThrows)
{
  const std::string path =
    WriteAtlasFile("atlas_unknown_element.xml", Replace(kValidAtlas, "</Atlas>", "<Bogus></Bogus></Atlas>"));
  AtlasDefinition atlas;
  EXPECT_THROW(atlas.InitFromXML(path), itk::ExceptionObject);
  std::remove(path.c_str());
}

namespace
{
std::string
InitFailureMessage(const std::string & path)
{
  AtlasDefinition atlas;
  try
  {
    atlas.InitFromXML(path);
  }
  catch (const itk::ExceptionObject & error)
  {
    return error.GetDescription();
  }
  return {};
}
} // namespace

TEST(AtlasDefinition, MalformedXMLReportsParserErrorAndLine)
{
  const std::string path = WriteAtlasFile("atlas_malformed_detail.xml", "<Atlas>\n <Prior>\n</Atlas>\n");
  const std::string message = InitFailureMessage(path);
  EXPECT_NE(message.find("mismatched tag"), std::string::npos) << message;
  EXPECT_NE(message.find("line 3"), std::string::npos) << message;
  std::remove(path.c_str());
}

TEST(AtlasDefinition, HandlerExceptionMessageIsPreserved)
{
  const std::string path = WriteAtlasFile("atlas_bad_weight_detail.xml", Replace(kValidAtlas, "1.5", "heavy"));
  const std::string message = InitFailureMessage(path);
  EXPECT_NE(message.find("Bad Weight given"), std::string::npos) << message;
  EXPECT_NE(message.find("heavy"), std::string::npos) << message;
  EXPECT_NE(message.find(path), std::string::npos) << message;
  std::remove(path.c_str());
}

TEST(AtlasDefinition, UnknownElementMessageNamesTheElement)
{
  const std::string path =
    WriteAtlasFile("atlas_unknown_element_detail.xml", Replace(kValidAtlas, "</Atlas>", "<Bogus></Bogus></Atlas>"));
  const std::string message = InitFailureMessage(path);
  EXPECT_NE(message.find("Bogus"), std::string::npos) << message;
  std::remove(path.c_str());
}

TEST(AtlasDefinition, EmptyFileThrows)
{
  const std::string path = WriteAtlasFile("atlas_empty.xml", "");
  AtlasDefinition   atlas;
  EXPECT_THROW(atlas.InitFromXML(path), itk::ExceptionObject);
  std::remove(path.c_str());
}
