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
/// \file Slicer3LandmarkIOGTest.cxx
/// \brief GTest port of Slicer3LandmarkIOTest.cxx (issue #403 — GTest migration).
///
/// Verifies that ReadSlicer3toITKLmk() and ReadLandmarkWeights() throw
/// itk::ExceptionObject when given non-existent file paths.  The original
/// test used manual try/catch; GTest's EXPECT_THROW makes the intent explicit
/// and produces named CTest entries.

#include <gtest/gtest.h>

#include "Slicer3LandmarkIO.h"

TEST(Slicer3LandmarkIO, ThrowsOnNonExistentLandmarkFile)
{
  // Expect an itk::ExceptionObject (or any derived class) when the file
  // path does not exist.  This mirrors the original catch block that
  // treated the exception as proof of correct error handling.
  EXPECT_THROW(ReadSlicer3toITKLmk("/this/file/does/not/exist.fcsv"), itk::ExceptionObject);
}

TEST(Slicer3LandmarkIO, ThrowsOnNonExistentLandmarkWeightFile)
{
  EXPECT_THROW(ReadLandmarkWeights("/this/file/does/not/exist/either.fcsv"), itk::ExceptionObject);
}
