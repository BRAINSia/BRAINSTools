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
#ifndef BRAINSThreadControl_h
#define BRAINSThreadControl_h

#include <itksys/SystemTools.hxx>
#include <sstream>
#include "itkMultiThreader.h"

namespace BRAINSUtils
{
/**
 * This class is designed so that
 * the ITK number of threads can be
 * adjusted to a different number of threads
 * during the running of this one call.
 * The desired number of threads
 * must be selected as part of the
 * construction process, and at
 * destruction, the original value is
 * restored.
 *
 * This type of functionality is needed
 * so that the shared libary version of
 * BRAINSFit does not globally change
 * the behavior of all other programs
 * in slicer.
 */
class StackPushITKDefaultNumberOfThreads
{
public:
  explicit StackPushITKDefaultNumberOfThreads(const int desiredCount);
  ~StackPushITKDefaultNumberOfThreads();
protected:
  StackPushITKDefaultNumberOfThreads();                                                 // Purposefully not implemented
  StackPushITKDefaultNumberOfThreads & operator=(StackPushITKDefaultNumberOfThreads &); // Purposefully not implemented

private:
  int m_originalThreadValue;
};
}

#endif // BRAINSThreadControl_h
