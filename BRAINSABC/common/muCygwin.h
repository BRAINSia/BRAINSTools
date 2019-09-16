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
// Include this file before including stuff from mu

#ifndef __muCygwin_h
#define __muCygwin_h

#include "itkMacro.h"
#include <exception>
#include <iostream>

// Cygwin exception handling work-around
#undef itkExceptionMacro
#define itkExceptionMacro(x)                                                                                           \
  {                                                                                                                    \
    std::cerr << "Exception: " x << std::endl;                                                                         \
    std::cerr << "Possibly crashing about now..." << std::endl;                                                        \
    throw "exc";                                                                                                       \
  }

// TODO wrap main so that uncaught exception does not crash program
#define MU_DEFINE_MAIN                                                                                                 \
  int main(int argc, char ** argv)                                                                                     \
  {                                                                                                                    \
    int r = 0;                                                                                                         \
    try                                                                                                                \
    {                                                                                                                  \
      r = _mu_main(argc, argv);                                                                                        \
    }                                                                                                                  \
    catch (itk::ExceptionObject & e)                                                                                   \
    {                                                                                                                  \
      std::cerr << e << std::endl;                                                                                     \
      return -1;                                                                                                       \
    }                                                                                                                  \
    catch (std::exception & e)                                                                                         \
    {                                                                                                                  \
      std::cerr << "Exception: " << e.what() << std::endl;                                                             \
      return -1;                                                                                                       \
    }                                                                                                                  \
    catch (char * s)                                                                                                   \
    {                                                                                                                  \
      std::cerr << "Exception: " << s << std::endl;                                                                    \
      return -1;                                                                                                       \
    }                                                                                                                  \
    return r;                                                                                                          \
  }

#endif
