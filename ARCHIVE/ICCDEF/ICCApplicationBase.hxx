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
#ifndef _ICCApplicationBase_txx
#define _ICCApplicationBase_txx

#include "ICCApplicationBase.h"

namespace itk
{
template <typename TPreprocessor, typename TRegistrator>
ICCApplicationBase<TPreprocessor, TRegistrator>
::ICCApplicationBase()
{
  m_Preprocessor = PreprocessorType::New();
  m_Registrator  = RegistratorType::New();
  m_OutDebug  = false;
}

template <typename TPreprocessor, typename TRegistrator>
void
ICCApplicationBase<TPreprocessor, TRegistrator>
::Execute()
{
  /**************************
    * Preprocess the images before registration
    **************************/

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
   **************************/
  if( this->GetOutDebug() )
    {
    std::cout << "Register the images ... " << std::endl;
    }

  try
    {
    this->InitializeRegistrator();
    m_Preprocessor = nullptr;
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
