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
#ifndef __itkMeshToMeshRegistrationMethod_hxx
#define __itkMeshToMeshRegistrationMethod_hxx

#include "itkMeshToMeshRegistrationMethod.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedMesh, typename TMovingMesh>
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::MeshToMeshRegistrationMethod()
{
  this->SetNumberOfRequiredOutputs(1); // for the Transform

  m_FixedMesh = nullptr;  // has to be provided by the user.
  m_MovingMesh = nullptr; // has to be provided by the user.
  m_Transform = nullptr;  // has to be provided by the user.
  m_Metric = nullptr;     // has to be provided by the user.
  m_Optimizer = nullptr;  // has to be provided by the user.

  m_InitialTransformParameters = ParametersType(1);
  m_LastTransformParameters = ParametersType(1);

  m_InitialTransformParameters.Fill(0.0f);
  m_LastTransformParameters.Fill(0.0f);

  TransformOutputPointer transformDecorator = static_cast<TransformOutputType *>(this->MakeOutput(0).GetPointer());

  this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

/**
 * Set the initial transform parameters
 */
template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::SetInitialTransformParameters(const ParametersType & param)
{
  m_InitialTransformParameters = param;
  this->Modified();
}

/**
 * Initialize by setting the interconnects between components.
 */
template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::Initialize() throw(ExceptionObject)
{
  if (!m_FixedMesh)
  {
    itkExceptionMacro(<< "FixedMesh is not present");
  }

  if (!m_MovingMesh)
  {
    itkExceptionMacro(<< "MovingMesh is not present");
  }

  if (!m_Metric)
  {
    itkExceptionMacro(<< "Metric is not present");
  }

  if (!m_Optimizer)
  {
    itkExceptionMacro(<< "Optimizer is not present");
  }

  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform is not present");
  }

  if (!m_Interpolator)
  {
    itkExceptionMacro(<< "Interpolator is not present");
  }

  // Setup the metric
  m_Metric->SetMovingMesh(m_MovingMesh);
  m_Metric->SetFixedMesh(m_FixedMesh);
  m_Metric->SetTransform(m_Transform);
  m_Metric->SetInterpolator(m_Interpolator);

  m_Metric->Initialize();

  // Setup the optimizer
  m_Optimizer->SetCostFunction(m_Metric);

  // Validate initial transform parameters
  if (m_InitialTransformParameters.Size() != m_Transform->GetNumberOfParameters())
  {
    itkExceptionMacro(<< "Size mismatch between initial parameter and transform");
  }

  m_Optimizer->SetInitialPosition(m_InitialTransformParameters);

  //
  // Connect the transform to the Decorator.
  //
  TransformOutputType * transformOutput = static_cast<TransformOutputType *>(this->ProcessObject::GetOutput(0));

  transformOutput->Set(m_Transform.GetPointer());
}

/**
 * PrintSelf
 */
template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Metric: " << m_Metric.GetPointer() << std::endl;
  os << indent << "Optimizer: " << m_Optimizer.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Fixed Mesh: " << m_FixedMesh.GetPointer() << std::endl;
  os << indent << "Moving Mesh: " << m_MovingMesh.GetPointer() << std::endl;
  os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
  os << indent << "Last    Transform Parameters: " << m_LastTransformParameters << std::endl;
}

template <typename TFixedMesh, typename TMovingMesh>
void
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::GenerateData()
{
  try
  {
    // initialize the interconnects between components
    this->Initialize();
  }
  catch (ExceptionObject & err)
  {
    m_LastTransformParameters = ParametersType(1);
    m_LastTransformParameters.Fill(0.0f);

    // pass exception to caller
    throw err;
  }

  try
  {
    // do the optimization
    m_Optimizer->StartOptimization();
  }
  catch (ExceptionObject & err)
  {
    // An error has occurred in the optimization.
    // Update the parameters
    m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

    // Pass exception to caller
    throw err;
  }

  // get the results
  m_LastTransformParameters = m_Optimizer->GetCurrentPosition();

  m_Transform->SetParameters(m_LastTransformParameters);
}

/**
 *  Get Output
 */
template <typename TFixedMesh, typename TMovingMesh>
const typename MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::TransformOutputType *
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::GetOutput() const
{
  return static_cast<const TransformOutputType *>(this->ProcessObject::GetOutput(0));
}

template <typename TFixedMesh, typename TMovingMesh>
DataObject::Pointer
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::MakeOutput(size_t output)
{
  switch (output)
  {
    case 0:
    {
      return static_cast<DataObject *>(TransformOutputType::New().GetPointer());
    }
    break;
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
      return nullptr;
  }
}

/**
 *
 */
template <typename TFixedMesh, typename TMovingMesh>
unsigned long
MeshToMeshRegistrationMethod<TFixedMesh, TMovingMesh>::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  // Some of the following should be removed once ivars are put in the
  // input and output lists

  if (m_Transform)
  {
    m = m_Transform->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_Metric)
  {
    m = m_Metric->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_Optimizer)
  {
    m = m_Optimizer->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_FixedMesh)
  {
    m = m_FixedMesh->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_MovingMesh)
  {
    m = m_MovingMesh->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  return mtime;
}
} // end namespace itk

#endif
