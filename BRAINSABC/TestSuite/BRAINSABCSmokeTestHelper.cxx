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
/// \file BRAINSABCSmokeTestHelper.cxx
/// \brief Prepares inputs for, and checks the output of, the BRAINSABC smoke test.
///
///   BRAINSABCSmokeTestHelper shrink <factor> <input> <output>
///   BRAINSABCSmokeTestHelper check <baselineLabels> <testLabels> <factor> <tolerance>

#include "itkBinShrinkImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

namespace
{
using FloatImageType = itk::Image<float, 3>;
using LabelImageType = itk::Image<unsigned char, 3>;
using LabelCounts = std::map<unsigned char, itk::SizeValueType>;

constexpr itk::SizeValueType MinimumExpectedCount = 100;

int
Shrink(const unsigned int factor, const std::string & inputFileName, const std::string & outputFileName)
{
  auto reader = itk::ImageFileReader<FloatImageType>::New();
  reader->SetFileName(inputFileName);

  auto shrinker = itk::BinShrinkImageFilter<FloatImageType, FloatImageType>::New();
  shrinker->SetInput(reader->GetOutput());
  shrinker->SetShrinkFactors(factor);

  auto writer = itk::ImageFileWriter<FloatImageType>::New();
  writer->SetInput(shrinker->GetOutput());
  writer->SetFileName(outputFileName);
  writer->Update();
  return EXIT_SUCCESS;
}

LabelImageType::Pointer
ReadLabels(const std::string & fileName)
{
  auto reader = itk::ImageFileReader<LabelImageType>::New();
  reader->SetFileName(fileName);
  reader->Update();
  return reader->GetOutput();
}

LabelCounts
CountLabels(const LabelImageType * image)
{
  LabelCounts counts;
  for (itk::ImageRegionConstIterator<LabelImageType> it(image, image->GetBufferedRegion()); !it.IsAtEnd(); ++it)
  {
    ++counts[it.Get()];
  }
  return counts;
}

// The shrunken test image has about 1/factor^3 of the baseline's voxels, so each label's count should
// stay near the baseline count divided by factor^3. The tolerance is a ratio, so 2 accepts 0.5x to 2x.
int
Check(const std::string & baselineFileName,
      const std::string & testFileName,
      const unsigned int  factor,
      const double        tolerance)
{
  const auto baseline = ReadLabels(baselineFileName);
  const auto test = ReadLabels(testFileName);

  int status = EXIT_SUCCESS;
  for (unsigned int d = 0; d < 3; ++d)
  {
    const auto expectedSize = baseline->GetLargestPossibleRegion().GetSize(d) / factor;
    if (test->GetLargestPossibleRegion().GetSize(d) != expectedSize)
    {
      std::cerr << "FAIL: dimension " << d << " is " << test->GetLargestPossibleRegion().GetSize(d) << ", expected "
                << expectedSize << '\n';
      status = EXIT_FAILURE;
    }
  }

  const LabelCounts baselineCounts = CountLabels(baseline);
  const LabelCounts testCounts = CountLabels(test);
  const double      scale = static_cast<double>(factor) * factor * factor;

  for (const auto & entry : testCounts)
  {
    if (baselineCounts.count(entry.first) == 0)
    {
      std::cerr << "FAIL: label " << static_cast<int>(entry.first) << " does not occur in the baseline\n";
      status = EXIT_FAILURE;
    }
  }

  for (const auto & entry : baselineCounts)
  {
    const double expected = static_cast<double>(entry.second) / scale;
    if (entry.first == 0 || expected < static_cast<double>(MinimumExpectedCount))
    {
      continue;
    }
    const auto   found = testCounts.find(entry.first);
    const double actual = found == testCounts.end() ? 0.0 : static_cast<double>(found->second);
    const double ratio = actual / expected;
    const bool   ok = ratio >= 1.0 / tolerance && ratio <= tolerance;
    std::cout << "label " << static_cast<int>(entry.first) << ": expected about " << expected << ", found " << actual
              << ", ratio " << ratio << (ok ? "" : "  FAIL") << '\n';
    if (!ok)
    {
      status = EXIT_FAILURE;
    }
  }
  return status;
}
} // namespace

int
main(int argc, char * argv[])
{
  try
  {
    const std::string command = argc > 1 ? argv[1] : "";
    if (command == "shrink" && argc == 5)
    {
      return Shrink(static_cast<unsigned int>(std::stoul(argv[2])), argv[3], argv[4]);
    }
    if (command == "check" && argc == 6)
    {
      return Check(argv[2], argv[3], static_cast<unsigned int>(std::stoul(argv[4])), std::stod(argv[5]));
    }
  }
  catch (const itk::ExceptionObject & error)
  {
    std::cerr << error << '\n';
    return EXIT_FAILURE;
  }
  catch (const std::exception & error)
  {
    std::cerr << "Invalid argument: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " shrink <factor> <input> <output>\n"
            << "  " << argv[0] << " check <baselineLabels> <testLabels> <factor> <tolerance>\n";
  return EXIT_FAILURE;
}
