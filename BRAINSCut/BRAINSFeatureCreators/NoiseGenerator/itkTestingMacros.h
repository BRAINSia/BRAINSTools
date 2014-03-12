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
#ifndef __itkTestingMacros_h
#define __itkTestingMacros_h

#define TRY_EXPECT_EXCEPTION( command ) \
  try \
    {  \
    std::cout << "Trying " << #command << std::endl; \
    command;  \
    std::cerr << "Failed to catch expected exception" << std::endl;  \
    return EXIT_FAILURE;  \
    }  \
  catch( itk::ExceptionObject & excp )  \
    {  \
    std::cout << "Caught expected exception" << std::endl;  \
    std::cout << excp << std::endl; \
    }

#define TRY_EXPECT_NO_EXCEPTION( command ) \
  try \
    {  \
    std::cout << "Trying " << #command << std::endl; \
    command;  \
    }  \
  catch( itk::ExceptionObject & excp )  \
    {  \
    std::cerr << excp << std::endl; \
    return EXIT_FAILURE;  \
    }

#define TEST_SET_GET( variable, command ) \
  if( variable.GetPointer() != command )   \
    {   \
    std::cerr << "Error in " << #command << std::endl; \
    std::cerr << "Expected " << variable.GetPointer() << std::endl; \
    std::cerr << "but got  " << command << std::endl; \
    return EXIT_FAILURE; \
    }

#define TEST_SET_GET_VALUE( variable, command ) \
  if( variable != command )   \
    {   \
    std::cerr << "Error in " << #command << std::endl; \
    std::cerr << "Expected " << variable << std::endl; \
    std::cerr << "but got  " << command << std::endl; \
    return EXIT_FAILURE; \
    }

#endif
