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
Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)

Copyright (c) University of Iowa Department of Radiology. All rights reserved.
See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#include <math.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include "itkImage.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImageFileReader.h"
#include "itkDOMNodeXMLReader.h"
#include "itkDOMNode.h"

#include "BRAINSLabelStatsCLP.h"

std::string GetXmlLabelName( std::string fileName, int label )
{
  // Set default return value
  std::string labelName = "Error";

  // Load the XML File
  itk::DOMNode::Pointer          labelXmlInfo;
  itk::DOMNodeXMLReader::Pointer xmlReader = itk::DOMNodeXMLReader::New();

  xmlReader->SetFileName( fileName );
  xmlReader->Update();

  // Find the labels and extract name of specified label
  labelXmlInfo = xmlReader->GetOutput();
  itk::DOMNode::ChildrenListType elementList;
  labelXmlInfo->GetAllChildren(elementList);
  for( std::vector<itk::DOMNode *>::iterator it = elementList.begin(); it != elementList.end(); ++it )
    {
    // std::cout << *it << std::endl;
    // std::cout << (*it)->GetName() << std::endl;
    // std::cout << (*it)->GetPath() << std::endl;
    if( (*it)->GetName() == "data" )
      {
      itk::DOMNode::ChildrenListType labelList;
      (*it)->GetAllChildren(labelList);
      for( std::vector<itk::DOMNode *>::iterator lt = labelList.begin(); lt != labelList.end(); ++lt )
        {
        // std::cout << *lt << std::endl;
        // std::cout << (*lt)->GetName() << std::endl;
        // std::cout << (*lt)->GetPath() << std::endl;
        std::string attributeValue = (*lt)->GetAttribute("index");
        int         currentLabel = atoi( attributeValue.c_str() );
        if( currentLabel == label )
          {
          itk::DOMTextNode::Pointer textNode = (*lt)->GetTextChild(0);
          labelName = textNode->GetText();
          }
        }
      }
    }

  return labelName;
}

std::string GetAntsLabelName( std::string fileName, int label )
{
  // Set default return value
  std::string labelName = "Error";

  std::string value, txtLabel;
  std::string x1, x2, y1, y2, z1, z2;

  std::ifstream labelFile( fileName.c_str() );

  if( labelFile.is_open() )
    {
    while( !labelFile.eof() )
      {
      labelFile >> value >> x1 >> y1 >> z1;
      labelFile >> x2 >> y2 >> z2;
      std::getline(labelFile, txtLabel);

      int currentLabel = atoi( value.c_str() );
      if( currentLabel == label )
        {
        unsigned first = txtLabel.find('"');
        unsigned last = txtLabel.rfind('"');
        labelName = txtLabel.substr(first + 1, last - first - 1);
        labelFile.close();
        return labelName;
        }
      }

    labelFile.close();
    }

  return labelName;
}

std::string GetLabelName( int mode, std::string fileName, int label )
{
  switch( mode )
    {
    case 1:
      {
      return GetXmlLabelName(fileName, label);
      }
      break;
    case 3:
      {
      return GetAntsLabelName(fileName, label);
      }
      break;
    default:
      return "UNKNOWN";
    }

  return "UNKNOWN";
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  if( ( imageVolume.length() == 0 ) &&
      ( labelVolume.length() == 0 ) )
    {
    std::cout << "Error: Both the image and label must be specified" << std::endl;
    return EXIT_FAILURE;
    }

  if( outputPrefixColumnNames.size() != outputPrefixColumnValues.size() )
    {
    std::cout << "Error: Number of column names must match number of values" << std::endl;
    return EXIT_FAILURE;
    }

  int mode = 0;
  if( labelFileType == "fslxml" )
    {
    mode = 1;
    }
  else if( labelFileType == "csv" )
    {
    mode = 2;
    }
  else if( labelFileType == "ants" )
    {
    mode = 3;
    }

  if( echoSwitch )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Image: " <<   imageVolume << std::endl;
    std::cout << "Label Map: " <<   labelVolume << std::endl;
    std::cout << "Label Name File: " <<   labelNameFile << std::endl;
    std::cout << "Column Prefix Names: ";
    for( size_t i = 0; i < outputPrefixColumnNames.size(); ++i )
      {
      std::cout << outputPrefixColumnNames[i] << ", ";
      }
    std::cout << std::endl;
    std::cout << "Column Prefix Values: ";
    for( size_t i = 0; i < outputPrefixColumnValues.size(); ++i )
      {
      std::cout << outputPrefixColumnValues[i] << ", ";
      }
    std::cout << std::endl;
    std::cout << "Label File Type: " <<   labelFileType << std::endl;
    std::cout << "Number of Histogram Bins: " <<   numberOfHistogramBins << std::endl;
    std::cout << "Define Min/Max Method: " <<   minMaxType << std::endl;
    std::cout << "User define min: " <<   userDefineMinimum << std::endl;
    std::cout << "User defined max: " <<   userDefineMaximum << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  typedef itk::Image<float, 3>            ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( imageVolume );
  imageReader->Update();

  ImageType::RegionType  imageRegion;
  ImageType::SizeType    imageSize;
  ImageType::SpacingType imageSpacing;
  ImageType::PointType   imageOrigin;
  imageRegion = imageReader->GetOutput()->GetLargestPossibleRegion();
  imageSize = imageRegion.GetSize();
  imageSpacing = imageReader->GetOutput()->GetSpacing();
  imageOrigin = imageReader->GetOutput()->GetOrigin();

  typedef itk::Image<short, 3>            LabelType;
  typedef itk::ImageFileReader<LabelType> LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( labelVolume );
  labelReader->Update();

  LabelType::RegionType  labelRegion;
  LabelType::SizeType    labelSize;
  LabelType::SpacingType labelSpacing;
  LabelType::PointType   labelOrigin;
  labelRegion = labelReader->GetOutput()->GetLargestPossibleRegion();
  labelSize = labelRegion.GetSize();
  labelSpacing = labelReader->GetOutput()->GetSpacing();
  labelOrigin = labelReader->GetOutput()->GetOrigin();
  // Check the Image and Label Map to Make sure they define the same space
  for( size_t i = 0; i < 3; ++i )
    {
    if( imageSize[i] != labelSize[i] )
      {
      std::cout << "Error: Image and label size do not match" << std::endl;
      std::cout << "Image: " << imageSize << std::endl;
      std::cout << "Label: " << labelSize << std::endl;
      return EXIT_FAILURE;
      }
    if( fabs(labelSpacing[i] - imageSpacing[i]) > 0.01 )
      {
      std::cout << "Error: Image and label spacing do not match" << std::endl;
      std::cout << "Image: " << imageSpacing << std::endl;
      std::cout << "Label: " << labelSpacing << std::endl;
      return EXIT_FAILURE;
      }
    if( fabs(labelOrigin[i] - imageOrigin[i]) > 0.01 )
      {
      std::cout << "Error: Image and label origin do not match" << std::endl;
      std::cout << "Image: " << imageOrigin << std::endl;
      std::cout << "Label: " << labelOrigin << std::endl;
      return EXIT_FAILURE;
      }
    }

  float minValue;
  float maxValue;
  bool  computeGlobalHistogram = false;

  if( minMaxType == "manual" )
    {
    minValue = userDefineMinimum;
    maxValue = userDefineMaximum;
    computeGlobalHistogram = true;
    }
  else if( minMaxType == "image" )
    {
    typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput( imageReader->GetOutput() );
    minMaxFilter->Update();
    minValue = minMaxFilter->GetMinimum();
    maxValue = minMaxFilter->GetMaximum();
    computeGlobalHistogram = true;
    }

  typedef itk::LabelStatisticsImageFilter<ImageType, LabelType> StatsFilterType;
  StatsFilterType::Pointer statsFilter = StatsFilterType::New();
  statsFilter->SetInput( imageReader->GetOutput() );
  statsFilter->SetLabelInput( labelReader->GetOutput() );
  if( computeGlobalHistogram )
    {
    statsFilter->UseHistogramsOn();
    statsFilter->SetHistogramParameters(numberOfHistogramBins, minValue, maxValue);
    }
  statsFilter->Update();

  typedef StatsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  typedef StatsFilterType::LabelPixelType                LabelPixelType;
  for( size_t i = 0; i < outputPrefixColumnNames.size(); ++i )
    {
    std::cout << outputPrefixColumnNames[i] << ", ";
    }
  std::cout << "Name, label, min, max, median, mean, stddev, var, sum, count" << std::endl;
  for( ValidLabelValuesType::const_iterator vIt = statsFilter->GetValidLabelValues().begin();
       vIt != statsFilter->GetValidLabelValues().end();
       ++vIt )
    {
    if( statsFilter->HasLabel(*vIt) )
      {
      LabelPixelType labelValue = *vIt;

      std::string labelName = GetLabelName(mode, labelNameFile, labelValue);
      for( size_t i = 0; i < outputPrefixColumnValues.size(); ++i )
        {
        std::cout << outputPrefixColumnValues[i] << ", ";
        }
      std::cout << labelName << ", ";
      std::cout << labelValue << ", ";
      std::cout << statsFilter->GetMinimum( labelValue ) << ", ";
      std::cout << statsFilter->GetMaximum( labelValue ) << ", ";
      float medianValue;

      if( !computeGlobalHistogram )
        {
        minValue = statsFilter->GetMinimum( labelValue );
        maxValue = statsFilter->GetMaximum( labelValue );
        StatsFilterType::Pointer labelStatsFilter = StatsFilterType::New();
        labelStatsFilter->SetInput( imageReader->GetOutput() );
        labelStatsFilter->SetLabelInput( labelReader->GetOutput() );
        labelStatsFilter->UseHistogramsOn();
        labelStatsFilter->SetHistogramParameters(numberOfHistogramBins, minValue, maxValue);
        labelStatsFilter->Update();
        medianValue = labelStatsFilter->GetMedian( labelValue );
        }
      else
        {
        medianValue = statsFilter->GetMedian( labelValue );
        }
      std::cout << medianValue << ", ";
      std::cout << statsFilter->GetMean( labelValue ) << ", ";
      std::cout << statsFilter->GetSigma( labelValue ) << ", ";
      std::cout << statsFilter->GetVariance( labelValue ) << ", ";
      std::cout << statsFilter->GetSum( labelValue ) << ", ";
      std::cout << statsFilter->GetCount( labelValue ) << std::endl;
      }
    }
  return EXIT_SUCCESS;
}
