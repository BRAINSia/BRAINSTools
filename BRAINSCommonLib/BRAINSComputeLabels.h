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
#ifndef BRAINSComputeLabels_h
#define BRAINSComputeLabels_h

#include <iostream>
#include <vector>
#include <itkImage.h>
#include <vnl/vnl_vector.h>
#include "ExtractSingleLargestRegion.h"
#include "itkMultiplyImageFilter.h"

#include "itkLabelStatisticsImageFilter.h"

using LabelCountMapType = std::map<size_t, size_t>;
using ByteImageType = itk::Image<unsigned char, 3>;
extern LabelCountMapType
GetMinLabelCount(ByteImageType::Pointer & labelsImage, const vnl_vector<unsigned int> & PriorLabelCodeVector);
// Labeling using maximum a posteriori, also do brain stripping using
// mathematical morphology and connected component
template <typename TProbabilityImage, typename TByteImage, typename TFloatingPrecision>
void
ComputeLabels(std::vector<typename TProbabilityImage::Pointer> & Posteriors,
              std::vector<bool> &                                PriorIsForegroundPriorVector,
              const vnl_vector<unsigned int> &                   PriorLabelCodeVector,
              typename TByteImage::Pointer &                     NonAirRegion,
              typename TByteImage::Pointer &                     DirtyLabels,
              typename TByteImage::Pointer &                     CleanedLabels,
              TFloatingPrecision                                 InclusionThreshold, // No thresholding = 0.0F
              const size_t                                       minLabelSizeAllowed)                                      // Allow zero sized labels = 0
{
  std::cout << "\nComputing labels..." << std::endl;

  std::map<size_t, size_t> reverseLabelMap;
  for (size_t i = 0; i < PriorLabelCodeVector.size(); ++i)
  {
    reverseLabelMap[PriorLabelCodeVector[i]] = i;
  }

  const unsigned int                           numClasses = Posteriors.size();
  const typename TProbabilityImage::RegionType region = Posteriors[0]->GetLargestPossibleRegion();
  DirtyLabels = TByteImage::New();
  DirtyLabels->CopyInformation(Posteriors[0]);
  DirtyLabels->SetRegions(region);
  DirtyLabels->Allocate();
  typename TByteImage::Pointer foregroundMask = TByteImage::New();
  foregroundMask->CopyInformation(Posteriors[0]);
  foregroundMask->SetRegions(region);
  foregroundMask->Allocate();

  size_t                   currentMinLabelSize = 0;
  constexpr unsigned short max_iterations = 10; // Prevent infinite looping, just fail
  unsigned short           current_iteration = 0;
  do
  {
    current_iteration++;
    if (current_iteration > max_iterations)
    {
      std::cout << "ERROR:  Infinite loop detected for auto-correction" << std::endl;
      std::cout << "        Check input images to ensure proper intializaiton was completed." << std::endl;
      exit(-1);
    }
    DirtyLabels->FillBuffer(0);
    foregroundMask->FillBuffer(0);
    using LocalLOOPITERTYPE = unsigned int;
    const typename TByteImage::SizeType size = DirtyLabels->GetLargestPossibleRegion().GetSize();
    {
      for (LocalLOOPITERTYPE kk = 0; kk < (LocalLOOPITERTYPE)size[2]; kk++)
      {
        for (LocalLOOPITERTYPE jj = 0; jj < (LocalLOOPITERTYPE)size[1]; jj++)
        {
          for (LocalLOOPITERTYPE ii = 0; ii < (LocalLOOPITERTYPE)size[0]; ii++)
          {
            const typename TProbabilityImage::IndexType currIndex = { { ii, jj, kk } };
            if (NonAirRegion->GetPixel(currIndex) == 0) // If outside the tissue
            // region, then set to
            // zero vIndex!
            {
              // INFO:  May want to specify this explicitly in the XML file for
              // the proper background value
              DirtyLabels->SetPixel(currIndex, 0); // This is implied by the
              // FillBuffer(0) above;
              continue;
            }

            TFloatingPrecision maxPosteriorClassValue = Posteriors[0]->GetPixel(currIndex);
            unsigned int       indexMaxPosteriorClassValue = 0;
            for (unsigned int iclass = 1; iclass < numClasses; iclass++)
            {
              const TFloatingPrecision currentPosteriorClassValue = Posteriors[iclass]->GetPixel(currIndex);
              if (currentPosteriorClassValue > maxPosteriorClassValue)
              {
                maxPosteriorClassValue = currentPosteriorClassValue;
                indexMaxPosteriorClassValue = iclass;
              }
            }

            {
              bool         fgflag = PriorIsForegroundPriorVector[indexMaxPosteriorClassValue];
              unsigned int label = 99;
              if (maxPosteriorClassValue > InclusionThreshold)
              {
                label = PriorLabelCodeVector[indexMaxPosteriorClassValue];
              }

              // Only use non-zero probabilities and foreground classes
              if (!fgflag || (maxPosteriorClassValue < 0.001))
              {
                fgflag = false; // If priors are zero or negative, then set the
                // fgflag back to false
              }
              DirtyLabels->SetPixel(currIndex, label);
              foregroundMask->SetPixel(currIndex, fgflag);
            }
          }
        }
      }
    }
    //
    LabelCountMapType currentLabelsMapCounts = GetMinLabelCount(DirtyLabels, PriorLabelCodeVector);
    currentMinLabelSize = currentLabelsMapCounts.begin()->second;
    for (typename LabelCountMapType::const_iterator it = currentLabelsMapCounts.begin();
         it != currentLabelsMapCounts.end();
         ++it)
    {
      const size_t currentLabelCount = it->second;
      currentMinLabelSize = std::min<size_t>(currentMinLabelSize, currentLabelCount);
      if (currentLabelCount < minLabelSizeAllowed)
      {
        std::cout << "\n\nWARNING:  Increasing importance of label " << it->first
                  << "\n            because too few samples found."
                  << "\n           " << currentLabelCount << " < " << minLabelSizeAllowed << std::endl;
        // Multiply this prior by 1.1 to increase it's importance.
        using MultiplyFilterType = itk::MultiplyImageFilter<TProbabilityImage, TProbabilityImage>;
        typename MultiplyFilterType::Pointer filter = MultiplyFilterType::New();
        filter->SetInput1(Posteriors[reverseLabelMap[it->first]]);
        filter->SetInput2(1.1); // Multiply by 1.1 to increase it's importance
        filter->Update();
        Posteriors[reverseLabelMap[it->first]] = filter->GetOutput();
      }
    }
  } while (currentMinLabelSize < minLabelSizeAllowed);
  // CleanedLabels=ExtractSingleLargestRegionFromMask(foregroundMask,2,2,1,DirtyLabels);
  CleanedLabels = ExtractSingleLargestRegionFromMask(foregroundMask, 0, 0, 0, DirtyLabels);
}

#endif // BRAINSComputeLabels_h
