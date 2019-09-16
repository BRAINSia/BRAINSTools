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
//
// //////////////////////////////////////////////////////////////////////////////
//
// Handles log messages to terminal and disk, follows singleton pattern
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 2/2004

#ifndef __Log_h
#define __Log_h

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

namespace mu
{
/** \class Log
 */
class Log
{
public:
  static Log *
  GetInstance();

  void
  CloseFile();

  // Enable / disable writes to terminal
  inline void
  EchoOn()
  {
    m_EchoFlag = true;
  }

  inline void
  EchoOff()
  {
    m_EchoFlag = false;
  }

  void
  SetOutputFileName(const char * s);

  void
  SetOutputFileName(const std::string & s);

  void
  WriteString(const char * s);

  void
  WriteString(const std::string & s);

  std::ofstream &
  GetFileObject()
  {
    return m_Output;
  }

private:
  // Restrict access to constructors
  Log();
  ~Log();
  Log(const Log & l);

  bool m_EchoFlag;

  std::ofstream m_Output;

  std::string m_OutputFileName;
};
} // namespace mu

// Allows declarations such as: muLogMacro(<< "Message: " << 1.1234);
#define muLogMacro(x)                                                                                                  \
  {                                                                                                                    \
    std::ostringstream outss;                                                                                          \
    outss << "" x << std::ends;                                                                                        \
    (mu::Log::GetInstance())->WriteString(outss.str().c_str());                                                        \
  }

#endif
