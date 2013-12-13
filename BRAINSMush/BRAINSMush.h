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
const int Dimension = 3;
}

namespace
{
typedef unsigned char                InputPixelType;
typedef float                        PixelType;
typedef itk::Image<PixelType, 3>     ImageType;
typedef signed short                 MaskPixelType;
typedef itk::Image<MaskPixelType, 3> MaskImageType;
typedef MaskImageType::IndexType     MaskIndexType;

typedef itk::ImageFileWriter<MaskImageType> MaskImageWriterType;

typedef itk::ImageFileReader<ImageType>     ReaderType;
typedef itk::ImageFileReader<MaskImageType> MaskReaderType;

typedef itk::ImageRegionConstIterator<ImageType>     ConstIteratorType;
typedef itk::ImageRegionIterator<MaskImageType>      MaskIteratorType;
typedef itk::ImageRegionConstIterator<MaskImageType> ConstMaskIteratorType;

typedef itk::BinaryBallStructuringElement<InputPixelType,
                                          BRAINSMush::Dimension> StructuringElementType;
}

ImageType::Pointer LoadImage(std::string);

MaskImageType::Pointer LoadMaskImage(std::string);

MaskImageType::Pointer GenerateInitializerRegion(ImageType::Pointer & referenceImage, std::vector<int> boundingBoxSize,
                                                 std::vector<int> boundingBoxStart);

ImageType::Pointer MixtureOptimizer(ImageType::Pointer & firstImage, ImageType::Pointer & secondImage,
                                    MaskImageType::Pointer & maskImage, double desiredMean, double desiredVariance,
                                    std::string outputWeightsFile);

void GenerateBrainVolume(ImageType::Pointer & firstImage, ImageType::Pointer & secondImage,
                         MaskImageType::Pointer & maskImage, std::string inputMaskVolume, double desiredMean,
                         double desiredVariance, double lowerThresholdFactor, double upperThresholdFactor,
                         std::vector<int> boundingBoxSize, std::vector<int> boundingBoxStart,
                         //  std::vector<int> seed,
                         std::string outputVolume,
                         //  std::string outputMask,
                         std::string outputWeightsFile, MaskImageType::Pointer & resultImage);

#endif /* __BrainMUSH_h__ */
