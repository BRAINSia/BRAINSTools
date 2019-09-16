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
#ifndef _ICCApplicationBase_h
#define _ICCApplicationBase_h

#include "itkObjectFactory.h"
#include "itkObject.h"

namespace itk
{
/** \class ICCApplicationBase
 *
 * This class ties together an input parser, a preprocessor,
 * a registrator components to
 * form a deformable registration/atlas segmentation application.
 *
 */
template <typename TPreprocessor, typename TRegistrator>
class ICCApplicationBase : public Object
{
public:
  /** Standard class type alias. */
  using Self = ICCApplicationBase;
  using Superclass = Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MIMApplication, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Preprocessor type. */
  using PreprocessorType = TPreprocessor;
  using PreprocessorPointer = typename PreprocessorType::Pointer;

  /** Registrator type. */
  using RegistratorType = TRegistrator;
  using RegistratorPointer = typename RegistratorType::Pointer;

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  RegistratorType *
  GetRegistratorType(void)
  {
    return m_Registrator;
  }

  /** Execute the application. */
  virtual void
  Execute();

protected:
  ICCApplicationBase();
  ~ICCApplicationBase() override {}

  /*** Initialize the preprocessor */
  virtual void
  InitializePreprocessor()
  {}

  /*** Initialize the registrator  */
  virtual void
  InitializeRegistrator()
  {}

  PreprocessorPointer m_Preprocessor;
  RegistratorPointer  m_Registrator;
  bool                m_OutDebug;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "ICCApplicationBase.hxx"
#endif

#endif
