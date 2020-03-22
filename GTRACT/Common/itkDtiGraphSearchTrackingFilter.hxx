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

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDtiGraphSearchTrackingFilter_hxx
#define __itkDtiGraphSearchTrackingFilter_hxx

#include "vtkFloatArray.h"
#include "vtkPolyLine.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include <itkIOCommon.h>
#include "itkMetaDataObject.h"
#include "itkProgressAccumulator.h"
#include "itkImageMomentsCalculator.h"
#include "itkDtiGraphSearchTrackingFilter.h"
#include "algo.h"

#include <iostream>

namespace itk
{
template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
DtiGraphSearchTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>::DtiGraphSearchTrackingFilter()
  : DtiTrackingFilterBase<TTensorImageType, TAnisotropyImageType, TMaskImageType>::DtiTrackingFilterBase()
{
  this->m_UseRandomWalk = false;
  this->m_MaximumBranches = 5;
  this->m_CurvatureBranchAngle = 60.0;
  this->m_AnisotropyBranchingValue = 0.0;
  this->m_RandomWalkAngle = 45.0;
  this->m_RandomSeed = -1;
  this->m_RandomGenerator = RandomGeneratorType::New();
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
typename itk::Point<double, 3>
DtiGraphSearchTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>::InitializeCenterOfMask()
{
  // ////////////////////////////////////////////////////////////////////////
  // Get the Center point for the Ending Region
  // ///////////////////////////////////////////////////////////////////////

  using MomentsCalculatorType = itk::ImageMomentsCalculator<typename Self::MaskImageType>;
  typename MomentsCalculatorType::Pointer centerOfMassFilter = MomentsCalculatorType::New();
  centerOfMassFilter->SetImage(this->m_EndingRegion);
  centerOfMassFilter->Compute();

  typename MomentsCalculatorType::VectorType com = centerOfMassFilter->GetCenterOfGravity();

  typename itk::Point<double, 3> sum;
  sum[0] = com[0];
  sum[1] = com[1];
  sum[2] = com[2];
  return sum;
}

template <typename TTensorImageType, typename TAnisotropyImageType, typename TMaskImageType>
void
DtiGraphSearchTrackingFilter<TTensorImageType, TAnisotropyImageType, TMaskImageType>::Update()
{
  using EigenValuesArrayType = typename Self::TensorImageType::PixelType::EigenValuesArrayType;
  using EigenVectorsMatrixType = typename Self::TensorImageType::PixelType::EigenVectorsMatrixType;

  float   anisotropy, anisotropySum;
  TVector vin(3), vout(3);

  typename Self::ContinuousIndexType index, tmpIndex;
  bool                               stop;
  typename Self::BranchListType      branchList;
  typename Self::SeedListType        newSeeds;
  typename Self::DirectionListType   newDirections;

  const double inRadians = this->pi / 180.0;
  double       curvatureBranchAngle = std::cos(this->m_CurvatureBranchAngle * inRadians);
  double       randomWalkAngle = std::cos(this->m_RandomWalkAngle / 2.0 * inRadians);

  // newSeeds.clear();
  // newDirections.clear();
  this->m_Output = vtkPolyData::New();
  this->m_TrackingDirections.clear();
  this->m_Seeds.clear();

  /* Initialize the random number generator */

  if (this->m_RandomSeed == -1)
  {
    this->m_RandomGenerator->Initialize();
  }
  else
  {
    this->m_RandomGenerator->Initialize(this->m_RandomSeed);
  }

  // std::cout << "AnisotropyImageRegion: " <<
  // this->m_AnisotropyImage->GetLargestPossibleRegion() << std::endl;
  // std::cout << "TensorImageRegion: " <<
  // this->m_TensorImage->GetLargestPossibleRegion() << std::endl;
  // std::cout << "EndingMaskRegion: " <<
  // this->m_EndingRegion->GetLargestPossibleRegion() << std::endl;;
  // std::cout << "StartingMaskRegion: " <<
  // this->m_StartingRegion->GetLargestPossibleRegion() << std::endl;;

  this->m_ScalarIP->SetInputImage(this->m_AnisotropyImage);
  this->m_VectorIP->SetInputImage(this->m_TensorImage);
  this->m_EndIP->SetInputImage(this->m_EndingRegion);
  typename Self::AnisotropyImageRegionType ImageRegion = this->m_AnisotropyImage->GetLargestPossibleRegion();

  this->m_StartIP->SetInputImage(this->m_StartingRegion);
  Self::InitializeSeeds();

  constexpr int maxFibers = 25;
  int           sizeOfOutput = 0;
  int           oldsizeOfOutput = 0;
  int           originalSeedCount = 0;

  // ////////////////////////////////////////////////////////////////////////
  // Get the Center Of Mass for the Ending Region
  // ///////////////////////////////////////////////////////////////////////
  typename itk::Point<double, 3> midPoint = this->InitializeCenterOfMask();
  double                         tmpPoint[3];
  tmpPoint[0] = midPoint[0];
  tmpPoint[1] = midPoint[1];
  tmpPoint[2] = midPoint[2];
  typename Self::ContinuousIndexType endP;

  this->MMToContinuousIndex(tmpPoint, endP);
  // std::cout << "Ending Bound Box: " << bb << std::endl;
  // std::cout << "Ending Center Point: " << endP << std::endl;
  while (!(this->m_Seeds.empty()))
  {
    originalSeedCount++;
    if (oldsizeOfOutput != sizeOfOutput)
    {
      oldsizeOfOutput = sizeOfOutput;
      std::cout << std::endl;
    }
    if (sizeOfOutput > maxFibers) // useful for some cases of debugging:
    {
      // break;
    }

    anisotropySum = 0.0;

    index = this->m_Seeds.back();
    this->m_Seeds.pop_back();
    vin = vout = this->m_TrackingDirections.back();
    this->m_TrackingDirections.pop_back();

    // std::cout << "Seed Index: " << index << " Seed Size: " <<
    // this->m_Seeds.size() << std::endl;
    // std::cout << "Direction: " << vin << std::endl;
    // std::cout << "Image Region: " << ImageRegion << std::endl;

    stop = false;
    // int order = 1;

    vtkPoints *     fiber = vtkPoints::New();
    vtkFloatArray * fiberTensors = vtkFloatArray::New();
    fiberTensors->SetName("Tensors");
    fiberTensors->SetNumberOfComponents(9);
    vtkFloatArray * fiberAnisotropy = vtkFloatArray::New();
    fiberAnisotropy->SetName("Anisotropy");
    vtkFloatArray * fiberAnisotropySum = vtkFloatArray::New();
    fiberAnisotropySum->SetName("Anisotropy-Sum");
    int currentPointId = 0;

    branchList.clear();

    // ////////////////////////////////////////////////////////////////////////
    // Tracking start from given 'index' and 'vout'
    while (!stop)
    {
      if (ImageRegion.IsInside(index))
      {
        anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex(index);
      }
      else
      {
        anisotropy = -1;
      }

      //
      // ////////////////////////////////////////////////////////////////////////
      // Evaluate the stopping criteria
      //
      // ////////////////////////////////////////////////////////////////////////
      bool isLoop = false;
      if (this->m_UseLoopDetection)
      {
        isLoop = Self::IsLoop(fiber);
      }
      bool outOfBounds = false;

      if ((currentPointId > (this->m_MaximumLength / this->m_StepSize)) || (anisotropy < this->m_AnisotropyThreshold) ||
          (isLoop) || (outOfBounds))
      // || ( branchList.size() > this->m_MaximumBranches) ) - Removed as a
      // stopping criteria
      {
        // std::cout << "Stop Track" << std::endl;
        //
        // ////////////////////////////////////////////////////////////////////////
        // some other conditions: (avrAI<this->m_MeanAI)
        //
        // ////////////////////////////////////////////////////////////////////////
        // Backup to the previous branch restart tracking
        if (!branchList.empty())
        {
          BranchPointType bp = branchList.back();
          branchList.pop_back();
          // fiber.resize(bp.m_DivergePoint);
          while (currentPointId > bp.m_DivergePoint)
          {
            fiberTensors->RemoveLastTuple();
            fiberAnisotropy->RemoveLastTuple();
            fiberAnisotropySum->RemoveLastTuple();
            currentPointId--;
          }

          // fiber->SetNumberOfPoints( currentPointId );

          vtkPoints * newfiber = vtkPoints::New();
          newfiber->SetNumberOfPoints(currentPointId);
          for (int i = 0; i < currentPointId; i++)
          {
            newfiber->SetPoint(i, fiber->GetPoint(i));
          }
          fiber->Delete();
          fiber = newfiber;

          vout = bp.m_Direction;
          double * p = fiber->GetPoint(currentPointId - 1);
          this->MMToContinuousIndex(p, index);
        }
        else
        {
          stop = true;
        }
      }
      else
      {
        //
        // ////////////////////////////////////////////////////////////////////////
        // forward propagating
        //
        // ////////////////////////////////////////////////////////////////////////
        anisotropy = this->m_ScalarIP->EvaluateAtContinuousIndex(index);
        if (currentPointId)
        {
          anisotropySum += anisotropy;
        }
        else
        {
          anisotropySum = anisotropy;
        }
        fiberAnisotropy->InsertNextValue(anisotropy);
        fiberAnisotropySum->InsertNextValue(anisotropySum);

        typename Self::PointType t;
        this->ContinuousIndexToMM(index, t);
        fiber->InsertNextPoint(t.GetDataPointer());
        currentPointId++;

        EigenValuesArrayType                eigenValues;
        EigenVectorsMatrixType              eigenVectors;
        typename Self::TensorImagePixelType tensorPixel = this->m_VectorIP->EvaluateAtContinuousIndex(index);

        TMatrix fullTensorPixel(3, 3);

        fullTensorPixel = Tensor2Matrix(tensorPixel);
        fiberTensors->InsertNextTypedTuple(fullTensorPixel.data_block());

        tensorPixel.ComputeEigenAnalysis(eigenValues, eigenVectors);

        //
        // ////////////////////////////////////////////////////////////////////////
        // Get two tracking vectors - Primary and Secondary Eigen Value
        //
        // ////////////////////////////////////////////////////////////////////////
        TVector e2(3);

        e2[0] = eigenVectors[2][0];
        e2[1] = eigenVectors[2][1];
        e2[2] = eigenVectors[2][2];
        if (dot_product(vin, e2) < 0)
        {
          e2 *= -1;
        }
        TVector e1(3);

        e1[0] = eigenVectors[1][0];
        e1[1] = eigenVectors[1][1];
        e1[2] = eigenVectors[1][2];
        if (dot_product(vin, e1) < 0)
        {
          e1 *= -1;
        }

        // std::cout << "Forward Track Directions: (" << e2 << ") and (" << e1
        // << ")" << std::endl;

        //
        // ////////////////////////////////////////////////////////////////////////
        // Add a branch points - Check Criteria for Branching
        //
        // ////////////////////////////////////////////////////////////////////////

        if (((anisotropy < this->m_AnisotropyBranchingValue) || (dot_product(e2, vin) < curvatureBranchAngle)) &&
            (branchList.size() <= this->m_MaximumBranches))
        {
          // std::cout << "Branch Point" << std::endl;

          BranchPointType bp;
          // bp.m_AI = 0;        bp.m_Length = 0;
          bp.m_DivergePoint = currentPointId;

          if (this->m_UseRandomWalk)
          {
            TVector v(3);

            v[0] = endP[0] - index[0];
            v[1] = endP[1] - index[1];
            v[2] = endP[2] - index[2];
            v.normalize();

            // std::cout << "Random Walk Direction: " << v << std::endl;

            double x, y, z;
            x = (0.5 - this->m_RandomGenerator->GetVariateWithOpenRange()) * 2.0;
            y = (0.5 - this->m_RandomGenerator->GetVariateWithOpenRange()) * 2.0;
            z = (0.5 - this->m_RandomGenerator->GetVariateWithOpenRange()) * 2.0;
            double m = std::sqrt((x * x) + (y * y) + (z * z));
            x /= m;
            y /= m;
            z /= m;

            // Scale the angle in radians 0...pi/2 to the range 0...1
            // for scaling of the random direction
            x *= (randomWalkAngle / (this->pi / 2.0));
            y *= (randomWalkAngle / (this->pi / 2.0));
            z *= (randomWalkAngle / (this->pi / 2.0));

            // gsl_ran_dir_3d(r,&x,&y,&z);
            TVector randDir(3);

            randDir[0] = x;
            randDir[1] = y;
            randDir[2] = z;
            if (dot_product(v, randDir) < 0)
            {
              randDir *= -1;
            }
            v += randDir;
            v.normalize();
            vout = v;
            bp.m_Direction = e2;
            branchList.push_back(bp);
            bp.m_Direction = e1;
            branchList.push_back(bp);
          }
          else
          {
            bp.m_Direction = e1;
            branchList.push_back(bp);
            vout = e2;
          }
        }
        else
        {
          // Using TEND????
          if (this->m_UseTend)
          {
            this->ApplyTensorDeflection(vin, fullTensorPixel, e2, vout);
          }
          else
          {
            vout = e2;
          }
        }
      }

      //
      // ////////////////////////////////////////////////////////////////////////
      // Calculate the new index
      this->StepIndex(tmpIndex, index, vout);
      bool backTrack = false;
      if (ImageRegion.IsInside(tmpIndex))
      {
        //
        // ////////////////////////////////////////////////////////////////////////
        // Check if we are in the ending region?
        //
        // ////////////////////////////////////////////////////////////////////////
        if ((this->m_EndIP->EvaluateAtContinuousIndex(tmpIndex) >= 0.5) &&
            (fiber->GetNumberOfPoints() / this->m_StepSize >= this->m_MinimumLength))
        {
          // Add Fiber to the Current Fiber Track //
          sizeOfOutput++;
          std::cerr << "Fiber (" << originalSeedCount << "  " << sizeOfOutput << "  " << currentPointId << "); ";
          this->AddFiberToOutput(fiber, fiberTensors);
          // // newSeeds.push_back(index);
          // // newDirections.push_back(vout);
          // std::cout << "Add Track" << std::endl;

          backTrack = true;
        }
      }
      else
      {
        backTrack = true;
        outOfBounds = true; // back up to a previous branch point, if any.
      }

      if (backTrack)
      {
        //
        // ////////////////////////////////////////////////////////////////////////
        // back tracking
        if (!branchList.empty())
        {
          BranchPointType bp = branchList.back();
          branchList.pop_back();
          // fiber.resize(bp.m_DivergePoint);
          while (currentPointId > bp.m_DivergePoint)
          {
            fiberTensors->RemoveLastTuple();
            fiberAnisotropy->RemoveLastTuple();
            fiberAnisotropySum->RemoveLastTuple();
            currentPointId--;
          }

          // fiber->SetNumberOfPoints( currentPointId );

          vtkPoints * newfiber = vtkPoints::New();
          newfiber->SetNumberOfPoints(currentPointId);
          for (int i = 0; i < currentPointId; i++)
          {
            newfiber->SetPoint(i, fiber->GetPoint(i));
          }
          fiber->Delete();
          fiber = newfiber;

          vout = bp.m_Direction;
          double *                           p = fiber->GetPoint(currentPointId - 1);
          typename Self::ContinuousIndexType prevIndex;
          this->MMToContinuousIndex(p, prevIndex);
          this->StepIndex(tmpIndex, prevIndex, vout);
        }
        else
        {
          stop = true;
        }
      }

      //
      // ////////////////////////////////////////////////////////////////////////
      // Reset the current index
      index = tmpIndex;
      vin = vout;

      // std::cout << "New Index: " << index << std::endl;
      // std::cout << "Vin: " << vin << std::endl;
    } // End Stop
  }   // End - Seeds not empty

  //  fiber.clear();
  //  branchList.clear();
  //  seeds               = newSeeds;
  //  directions  = newDirections;
}
} // end namespace itk
#endif
