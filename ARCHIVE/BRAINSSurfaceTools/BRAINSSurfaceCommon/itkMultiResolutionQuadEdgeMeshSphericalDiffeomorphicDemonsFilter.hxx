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
#ifndef __itkMultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_hxx
#define __itkMultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter_hxx

#include "itkMultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkProgressReporter.h"

#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkLinearInterpolateMeshFunction.h"

#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkDeformationFieldFromTransformMeshFilter.h"
#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkIdentityTransform.h"
#include "itkReplaceDestinationPointsQuadEdgeMeshFilter.h"
#include "itkDeformQuadEdgeMeshFilter.h"

namespace itk
{
template <typename TMesh>
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
  TMesh>::MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter()
{
  this->SetNumberOfIndexedInputs(8); // four resolution levels, two meshes on each
  this->SetNumberOfIndexedOutputs(1);

  this->SetNumberOfRequiredInputs(8); // four resolution levels, two meshes on each
  this->SetNumberOfRequiredOutputs(1);

  this->SetNthOutput(0, TMesh::New());

  this->m_SphereCenter.Fill(0.0);
  this->m_SphereRadius = 1.0;

  this->m_RigidTransform = TransformType::New();

  this->m_RigidOptimizer = RigidOptimizerType::New();

  this->m_CurrentLevelFixedMesh = MeshType::New();
  this->m_CurrentLevelMovingMesh = MeshType::New();
  this->m_CurrentLevelInitialFixedMesh = MeshType::New();

  this->m_CurrentLevelRigidlyMappedFixedMesh = MeshType::New();

  this->m_FinalDestinationPoints = DestinationPointSetType::New();

  this->m_DemonsRegistrationFilter = DemonsRegistrationFilterType::New();

  this->m_CurrentResolutionLevel = 0;
  this->m_NumberOfResolutionLevels = 4;

  this->m_SelfRegulatedMode = false;
  this->m_SelfStopMode = false;
}

template <typename TMesh>
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
  TMesh>::~MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter()
{}

template <typename TMesh>
DataObject::Pointer
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::MakeOutput(size_t idx)
{
  DataObject::Pointer output;

  switch (idx)
  {
    case 0:
    {
      output = (TMesh::New()).GetPointer();
    }
    break;
  }
  return output.GetPointer();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::SetFixedMesh(const MeshType * fixedMesh)
{
  itkDebugMacro("setting Fixed Mesh to " << fixedMesh);

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, const_cast<MeshType *>(fixedMesh));

  this->Modified();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::SetMovingMesh(const MeshType * movingMesh)
{
  itkDebugMacro("setting Moving Mesh to " << movingMesh);

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(1, const_cast<MeshType *>(movingMesh));

  this->Modified();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::SetFixedMesh(unsigned int     level,
                                                                                   const MeshType * fixedMesh)
{
  itkDebugMacro("setting Fixed Mesh to " << fixedMesh);

  const unsigned int inputNumber = 2 * level;

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(inputNumber, const_cast<MeshType *>(fixedMesh));

  this->Modified();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::SetMovingMesh(unsigned int     level,
                                                                                    const MeshType * movingMesh)
{
  itkDebugMacro("setting Moving Mesh to " << movingMesh);

  const unsigned int inputNumber = 2 * level + 1;

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(inputNumber, const_cast<MeshType *>(movingMesh));

  this->Modified();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::GenerateData()
{
  this->InitializeRigidRegistrationParameters();
  this->InitializeDemonsRegistrationParameters();
  this->PrepareCoarsestResolutionMeshes();

  while (this->m_CurrentResolutionLevel < this->m_NumberOfResolutionLevels)
  {
    std::cout << "RESOLUTION LEVEL = " << this->m_CurrentResolutionLevel << std::endl;
    this->ComputeRigidRegistration();
    this->RigidlyTransformPointsOfFixedMesh();
    this->ComputeDemonsRegistration();
    this->PrepareNextResolutionLevelMeshes();
    this->DeformNextResolutionLevelFixedMesh();
    this->m_CurrentResolutionLevel++;
  }
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::InitializeRigidRegistrationParameters()
{
  using ScalesType = RigidOptimizerType::ScalesType;

  const unsigned int numberOfTransformParameters = this->m_RigidTransform->GetNumberOfParameters();

  ScalesType parametersScale(numberOfTransformParameters);
  parametersScale[0] = 1.0;
  parametersScale[1] = 1.0;
  parametersScale[2] = 1.0;

  this->m_RigidOptimizer->SetScales(parametersScale);

  this->m_RigidOptimizer->MinimizeOn();
  this->m_RigidOptimizer->SetGradientMagnitudeTolerance(1e-6);
  this->m_RigidOptimizer->SetMinimumStepLength(1e-9);
  this->m_RigidOptimizer->SetRelaxationFactor(0.9);

  if (this->m_RigidRegistrationIterations.size() < this->m_NumberOfResolutionLevels)
  {
    itkExceptionMacro("Rigid registration iterations array size is smaller than number of iteration levels");
  }

  if (this->m_RigidRegistrationStepLength.size() < this->m_NumberOfResolutionLevels)
  {
    itkExceptionMacro("Rigid registration step length array size is smaller than number of iteration levels");
  }
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::InitializeDemonsRegistrationParameters()
{
  this->m_DemonsRegistrationFilter->SetSphereCenter(this->m_SphereCenter);
  this->m_DemonsRegistrationFilter->SetSphereRadius(this->m_SphereRadius);

  this->m_DemonsRegistrationFilter->SetEpsilon(0.016);
  this->m_DemonsRegistrationFilter->SetSigmaX(8.0);

  this->m_DemonsRegistrationFilter->SetLambda(1.0);
  this->m_DemonsRegistrationFilter->SetMetricSignificance(1.0);

  if (this->m_SmoothingIterations.size() < this->m_NumberOfResolutionLevels)
  {
    itkExceptionMacro("Smoothing iterations array size is smaller than number of iteration levels");
  }

  if (this->m_DemonsIterations.size() < this->m_NumberOfResolutionLevels)
  {
    itkExceptionMacro("Demons iterations array size is smaller than number of iteration levels");
  }

  if (this->m_SelfRegulatedMode)
  {
    this->m_DemonsRegistrationFilter->SelfRegulatedModeOn();
  }
  else
  {
    this->m_DemonsRegistrationFilter->SelfRegulatedModeOff();

    if (this->m_EpsilonValues.size() < this->m_NumberOfResolutionLevels)
    {
      itkExceptionMacro("Demons Epsilon array size is smaller than number of iteration levels");
    }
    if (this->m_SigmaXValues.size() < this->m_NumberOfResolutionLevels)
    {
      itkExceptionMacro("Demons SigmaX array size is smaller than number of iteration levels");
    }
  }

  if (this->m_SelfStopMode)
  {
    this->m_DemonsRegistrationFilter->SelfStopModeOn();
  }
  else
  {
    this->m_DemonsRegistrationFilter->SelfStopModeOff();
  }

  this->m_FinalDestinationPoints = DestinationPointSetType::New();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::PrepareCoarsestResolutionMeshes()
{
  this->m_CurrentResolutionLevel = 0;

  const TMesh * fixedMesh = dynamic_cast<MeshType *>(this->ProcessObject::GetInput(0));
  const TMesh * movingMesh = dynamic_cast<MeshType *>(this->ProcessObject::GetInput(1));

  this->m_CurrentLevelFixedMesh = fixedMesh;
  this->m_CurrentLevelMovingMesh = movingMesh;
  this->m_CurrentLevelInitialFixedMesh = fixedMesh;
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::PrepareNextResolutionLevelMeshes()
{
  if (this->m_CurrentResolutionLevel + 1 == this->m_NumberOfResolutionLevels)
  {
    return;
  }

  const unsigned int fixedInput = (this->m_CurrentResolutionLevel + 1) * 2;
  const unsigned int movingInput = (this->m_CurrentResolutionLevel + 1) * 2 + 1;

  //
  //   Prepare meshes for next resolution level
  //
  const TMesh * fixedMesh = dynamic_cast<MeshType *>(this->ProcessObject::GetInput(fixedInput));
  const TMesh * movingMesh = dynamic_cast<MeshType *>(this->ProcessObject::GetInput(movingInput));

  this->m_NextLevelFixedMesh = fixedMesh;
  this->m_NextLevelMovingMesh = movingMesh;
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::ComputeRigidRegistration()
{
  this->SetRigidTransformToIdentity();

  this->m_RegistrationMode = RIGID;
  std::cout << "RIGID" << std::endl;

  this->m_RigidOptimizer->SetNumberOfIterations(this->m_RigidRegistrationIterations[this->m_CurrentResolutionLevel]);

  this->m_RigidOptimizer->SetMaximumStepLength(this->m_RigidRegistrationStepLength[this->m_CurrentResolutionLevel]);

#ifdef USE_VTK
  this->m_RegistrationMonitor->SetResolutionLevel(this->m_CurrentResolutionLevel);
  this->m_RegistrationMonitor->Observe(this->GetRigidOptimizer());
  this->m_RegistrationMonitor->ObserveData(this->GetRigidTransform(), this->m_CurrentLevelFixedMesh);

#endif

  using RegistrationType = MeshToMeshRegistrationMethod<MeshType, MeshType>;

  using MetricType = MeanSquaresMeshToMeshMetric<MeshType, MeshType>;

  typename RegistrationType::Pointer registration = RegistrationType::New();
  // registration->InPlaceOn();

  typename MetricType::Pointer metric = MetricType::New();

  registration->SetMetric(metric);

  registration->SetFixedMesh(this->m_CurrentLevelFixedMesh);
  registration->SetMovingMesh(this->m_CurrentLevelMovingMesh);

  registration->SetTransform(this->m_RigidTransform);

  using InterpolatorType = LinearInterpolateMeshFunction<MeshType>;

  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  registration->SetInterpolator(interpolator);

  const unsigned int numberOfTransformParameters = this->m_RigidTransform->GetNumberOfParameters();

  using ParametersType = typename TransformType::ParametersType;
  ParametersType parameters(numberOfTransformParameters);

  parameters = this->m_RigidTransform->GetParameters();

  registration->SetInitialTransformParameters(parameters);

  registration->SetOptimizer(this->m_RigidOptimizer);

  try
  {
    registration->Update();
  }
  catch (ExceptionObject & e)
  {
    std::cerr << "Registration failed" << std::endl;
    std::cerr << "Reason " << e << std::endl;
    throw e;
  }

  ParametersType finalParameters = registration->GetLastTransformParameters();

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << m_RigidOptimizer->GetValue() << std::endl;

  this->m_RigidTransform->SetParameters(finalParameters);
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::RigidlyTransformPointsOfFixedMesh()
{
  this->m_CurrentLevelRigidlyMappedFixedMesh = MeshType::New();

  CopyMeshToMesh<MeshType, MeshType>(this->m_CurrentLevelFixedMesh, this->m_CurrentLevelRigidlyMappedFixedMesh);

  PointsContainer * fixedPoints = this->m_CurrentLevelRigidlyMappedFixedMesh->GetPoints();

  PointsContainerIterator fixedPointItr = fixedPoints->Begin();
  PointsContainerIterator fixedPointEnd = fixedPoints->End();

  while (fixedPointItr != fixedPointEnd)
  {
    fixedPointItr.Value().SetPoint(this->m_RigidTransform->TransformPoint(fixedPointItr.Value()));
    ++fixedPointItr;
  }

  this->m_CurrentLevelRigidlyMappedFixedMesh->Modified();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::SetRigidTransformToIdentity()
{
  this->m_RigidTransform->SetIdentity();
}

template <typename TMesh>
typename MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::RegistrationModeType
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::GetRegistrationMode() const
{
  return this->m_RegistrationMode;
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::ComputeDemonsRegistration()
{
  // This is needed for the proper visual monitoring of the Demons registration
  this->SetRigidTransformToIdentity();

  this->m_RegistrationMode = DEFORMABLE;
  std::cout << "DEMONS" << std::endl;

  this->m_DemonsRegistrationFilter->SetMaximumNumberOfSmoothingIterations(
    this->m_SmoothingIterations[this->m_CurrentResolutionLevel]);

  this->m_DemonsRegistrationFilter->SetMaximumNumberOfIterations(
    this->m_DemonsIterations[this->m_CurrentResolutionLevel]);

  this->m_DemonsRegistrationFilter->SetEpsilon(this->m_EpsilonValues[this->m_CurrentResolutionLevel]);
  std::cout << "Epsilon: " << this->m_EpsilonValues[this->m_CurrentResolutionLevel] << std::endl;

  this->m_DemonsRegistrationFilter->SetSigmaX(this->m_SigmaXValues[this->m_CurrentResolutionLevel]);
  std::cout << "SigmaX: " << this->m_SigmaXValues[this->m_CurrentResolutionLevel] << std::endl;

  this->m_DemonsRegistrationFilter->SetMetricSignificance(this->m_MetricSignificances[this->m_CurrentResolutionLevel]);
  std::cout << "Metric Significance: " << this->m_MetricSignificances[this->m_CurrentResolutionLevel] << std::endl;

#ifdef USE_VTK
  this->m_RegistrationMonitor->Observe(this->GetDemonsRegistrationFilter());
  this->m_RegistrationMonitor->ObserveData(this->GetRigidTransform(), this->GetCurrentDestinationPoints());
#endif

  this->m_DemonsRegistrationFilter->SetFixedMesh(this->m_CurrentLevelRigidlyMappedFixedMesh);
  this->m_DemonsRegistrationFilter->SetMovingMesh(this->m_CurrentLevelMovingMesh);

  try
  {
    this->m_DemonsRegistrationFilter->Update();
  }
  catch (ExceptionObject & exp)
  {
    std::cerr << exp << std::endl;
    throw exp;
  }
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::DeformNextResolutionLevelFixedMesh()
{
  typename DestinationPointSetType::Pointer currentDestinationPoints =
    this->m_DemonsRegistrationFilter->GetFinalDestinationPoints();

  currentDestinationPoints->DisconnectPipeline();

  // Create the new destination points for the following iteration
  this->m_DemonsRegistrationFilter->MakeOutput(2);

  this->m_FinalDestinationPoints = currentDestinationPoints;

  if (this->m_CurrentResolutionLevel + 1 == this->m_NumberOfResolutionLevels)
  {
    return;
  }

  //
  //   Deform fixed new mesh using current fixed mesh and destination points
  //
  using DeformFilterType = DeformQuadEdgeMeshFilter<MeshType, MeshType, DestinationPointSetType>;

  typename DeformFilterType::Pointer deformFilter = DeformFilterType::New();

  deformFilter->SetInput(this->m_NextLevelFixedMesh);
  deformFilter->SetReferenceMesh(this->m_CurrentLevelInitialFixedMesh);

  deformFilter->SetDestinationPoints(currentDestinationPoints);

  deformFilter->SetSphereRadius(this->m_SphereRadius);
  deformFilter->SetSphereCenter(this->m_SphereCenter);

  try
  {
    deformFilter->Update();
  }
  catch (ExceptionObject & excp)
  {
    std::cerr << excp << std::endl;
    throw excp;
  }

  this->m_CurrentLevelDemonsMappedFixedMesh = deformFilter->GetOutput();

  this->m_CurrentLevelDemonsMappedFixedMesh->DisconnectPipeline();

  this->m_CurrentLevelInitialFixedMesh = this->m_NextLevelFixedMesh;

  this->m_CurrentLevelFixedMesh = this->m_CurrentLevelDemonsMappedFixedMesh;

  this->m_CurrentLevelMovingMesh = this->m_NextLevelMovingMesh;
}

template <typename TMesh>
const typename MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::DestinationPointSetType *
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::GetCurrentDestinationPoints() const
{
  return this->m_DemonsRegistrationFilter->GetFinalDestinationPoints();
}

template <typename TMesh>
void
MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<TMesh>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Sphere center: " << this->m_SphereCenter << std::endl;
  os << "Sphere radius: " << this->m_SphereRadius << std::endl;
}
} // namespace itk

#endif
