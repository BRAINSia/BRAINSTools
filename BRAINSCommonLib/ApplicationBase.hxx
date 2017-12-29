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
/*=========================================================================
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile$
 *  Language:  C++
 *
 *  Copyright (c) 2002 Insight Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/
#ifndef __ApplicationBase_hxx
#define __ApplicationBase_hxx

#include "ApplicationBase.h"

namespace itk
{
template <typename TParser, typename TPreprocessor, typename TRegistrator>
ApplicationBase<TParser, TPreprocessor, TRegistrator>
::ApplicationBase()
{
  m_Parser       = ParserType::New();
  m_Preprocessor = PreprocessorType::New();
  m_Registrator  = RegistratorType::New();
  m_OutDebug  = false;
}

template <typename TParser, typename TPreprocessor, typename TRegistrator>
void
ApplicationBase<TParser, TPreprocessor, TRegistrator>
::Execute()
{
  /**************************
    * Parse input
    */
  if( this->GetOutDebug() )
    {
    std::cout << "Parsing input ... " << std::endl;
    }

  try
    {
    this->InitializeParser();
    m_Parser->Execute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
    }
  catch( ... )
    {
    std::cout << "Error occurred during input parsing." << std::endl;
    throw;
    }

  /**************************
    * Preprocess the images before registration
    */

  if( this->GetOutDebug() )
    {
    std::cout << "Preprocess the images ... " << std::endl;
    }

  try
    {
    this->InitializePreprocessor();
    m_Preprocessor->Execute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
    }
  catch( ... )
    {
    std::cout << "Error occured during preprocessing." << std::endl;
    throw;
    }

  /**************************
    * Registered the processed images
    */
  if( this->GetOutDebug() )
    {
    std::cout << "Register the images ... " << std::endl;
    }

  try
    {
    this->InitializeRegistrator();
    m_Preprocessor = nullptr;
    m_Parser = nullptr;
    m_Registrator->Execute();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cout << "Caught an ITK exception: " << std::endl;
    std::cout << err << " " << __FILE__ << " " << __LINE__ << std::endl;
    throw;
    }
  catch( ... )
    {
    std::cout << "Error occured during registration" << std::endl;
    throw;
    }
}
}   // namespace itk

#endif
