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
 *  Program:   BRAINS (Brain Research: Analysis of Images, Networks, and
 * Systems)
 *  Module:    $RCSfile: $
 *  Language:  TCL
 *  Date:      $Date: 2006/03/29 14:53:40 $
 *  Version:   $Revision: 1.9 $
 *
 *  Copyright (c) Iowa Mental Health Clinical Research Center. All rights
 * reserved.
 *  See BRAINSCopyright.txt or
 * http://www.psychiatry.uiowa.edu/HTML/Copyright.html
 *  for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/

#ifndef __BrainMUSH_h__
#define __BrainMUSH_h__

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "itkImage.h"

#include "itkImageFileReader.h"

#include <itkLabelStatisticsImageFilter.h>
#include "itkMatrix.h"
#include "itkImageFileWriter.h"

#include "itkImageIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkLevenbergMarquardtOptimizer.h"
#include "itkMixtureStatisticCostFunction.h"

namespace BRAINSMush
{
constexpr int Dimension = 3;
}

namespace
{
using InputPixelType = unsigned char;
using PixelType = float;
using ImageType = itk::Image<PixelType, 3>;
using MaskPixelType = signed short;
using MaskImageType = itk::Image<MaskPixelType, 3>;
using MaskIndexType = MaskImageType::IndexType;

using MaskImageWriterType = itk::ImageFileWriter<MaskImageType>;

using ReaderType = itk::ImageFileReader<ImageType>;
using MaskReaderType = itk::ImageFileReader<MaskImageType>;

using ConstIteratorType = itk::ImageRegionConstIterator<ImageType>;
using MaskIteratorType = itk::ImageRegionIterator<MaskImageType>;
using ConstMaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;

using StructuringElementType = itk::BinaryBallStructuringElement<InputPixelType, BRAINSMush::Dimension>;
} // namespace

MaskImageType::Pointer
GenerateInitializerRegion(ImageType::Pointer & referenceImage,
                          std::vector<int>     boundingBoxSize,
                          std::vector<int>     boundingBoxStart);

ImageType::Pointer
MixtureOptimizer(ImageType::Pointer &     firstImage,
                 ImageType::Pointer &     secondImage,
                 MaskImageType::Pointer & maskImage,
                 double                   desiredMean,
                 double                   desiredVariance,
                 std::string              outputWeightsFile);

void
GenerateBrainVolume(ImageType::Pointer &     firstImage,
                    ImageType::Pointer &     secondImage,
                    MaskImageType::Pointer & maskImage,
                    std::string              inputMaskVolume,
                    double                   desiredMean,
                    double                   desiredVariance,
                    double                   lowerThresholdFactor,
                    double                   upperThresholdFactor,
                    std::vector<int>         boundingBoxSize,
                    std::vector<int>         boundingBoxStart,
                    //  std::vector<int> seed,
                    std::string outputVolume,
                    //  std::string outputMask,
                    std::string              outputWeightsFile,
                    MaskImageType::Pointer & resultImage);

#endif /* __BrainMUSH_h__ */
