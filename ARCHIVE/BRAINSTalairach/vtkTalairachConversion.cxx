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
#include "vtkTalairachConversion.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "itkPoint.h"
#include "vtkStructuredGrid.h"

#include "vtkPoints.h"

#define PR(x) \
  std::cout << #x " = " << x << "\n"; // a simple print macro for
                                      // use when debugging

#define TALAIRACH_X_POINTS 9
#define TALAIRACH_Y_POINTS 12
#define TALAIRACH_Z_POINTS 15

vtkStandardNewMacro(vtkTalairachConversion);

vtkTalairachConversion::vtkTalairachConversion()
{
  this->MaskImage = ImageType::New();
  this->TalairachGrid = nullptr;
  this->SegmentationMode = false;
  this->HemisphereMode = both;
}

vtkTalairachConversion::~vtkTalairachConversion()
{
  if (this->MaskImage)
  {
    this->MaskImage->Delete();
    this->MaskImage = nullptr;
  }
  if (this->TalairachGrid)
  {
    this->TalairachGrid->Delete();
    this->TalairachGrid = nullptr;
  }
}

void
vtkTalairachConversion::Initialize()
{
  this->TalairachBoxList.clear();
  this->MaskImage->Initialize();
  if (this->TalairachGrid)
  {
    this->TalairachGrid->Delete();
    this->TalairachGrid = nullptr;
  }
  this->SegmentationMode = false;
  this->HemisphereMode = both;
}

int
vtkTalairachConversion::AddTalairachBox(std::string talairachBox)
{
  int index = TalairachBoxList.size();

  TalairachBoxList.push_back(talairachBox);
  return index;
}

void
vtkTalairachConversion::RemoveTalairachBox(std::string talairachBox)
{
  TalairachBoxList.remove(talairachBox);
}

void
vtkTalairachConversion::RemoveTalairachBox(int index)
{
  std::list<std::string>::iterator it;
  int                              currentIndex = 0;

  for (it = TalairachBoxList.begin(); it != TalairachBoxList.end(); ++it)
  {
    if (currentIndex == index)
    {
      TalairachBoxList.remove(*it);
      break;
    }
    ++currentIndex;
  }
}

void
vtkTalairachConversion::EraseTalairachBoxList()
{
  this->TalairachBoxList.clear();
}

void
vtkTalairachConversion::SetHemisphereMode(int mode)
{
  this->HemisphereMode = mode;
}

void
vtkTalairachConversion::SetHemisphereModeBoth()
{
  this->HemisphereMode = both;
}

void
vtkTalairachConversion::SetHemisphereModeLeft()
{
  this->HemisphereMode = left;
}

void
vtkTalairachConversion::SetHemisphereModeRight()
{
  this->HemisphereMode = right;
}

int
vtkTalairachConversion::GetHemisphereMode() const
{
  return this->HemisphereMode;
}

void
vtkTalairachConversion::SetSegmentationMode(bool mode)
{
  this->SegmentationMode = mode;
}

void
vtkTalairachConversion::SetSegmentationModeOn()
{
  this->SegmentationMode = true;
}

void
vtkTalairachConversion::SetSegmentationModeOff()
{
  this->SegmentationMode = false;
}

bool
vtkTalairachConversion::GetSegmentationMode() const
{
  return this->SegmentationMode;
}

void
vtkTalairachConversion::PrintSelf(ostream & os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Talairach Grid: " << std::endl;
  // this->TalairachGrid->PrintSelf(os, indent);
  os << indent << "Mask Image: " << std::endl;
  // this->MaskImage->PrintSelf(os, indent);
  os << indent << "Segmentation Mode: " << this->SegmentationMode << std::endl;
  os << indent << "Hemisphere Mode: ";

  switch (this->HemisphereMode)
  {
    case right:
      os << "Right" << std::endl;
      break;
    case left:
      os << "Left" << std::endl;
      break;
    case both:
      os << "Both" << std::endl;
      break;
  }
}

void
vtkTalairachConversion::SetTalairachGrid(vtkStructuredGrid * grid)
{
  this->TalairachGrid = grid;
}

vtkStructuredGrid *
vtkTalairachConversion::GetTalairachGrid()
{
  return this->TalairachGrid;
}

void
vtkTalairachConversion::ProcessBOX(bool _left)
{
  std::list<std::string>::iterator it;

  for (it = TalairachBoxList.begin(); it != TalairachBoxList.end(); ++it)
  {
    /* Requested information is 3 alphanumeric coordinate pairs
     * given in the form of six whitespace-delimited tokens per
     * line; this vector will store one line's worth at a time */
    std::vector<std::string> tokens;

    tokens.clear();
    for (int i = 0; i < 6; i++)
    {
      std::string       buf;
      std::stringstream ss(*it);
      while (ss >> buf)
      {
        tokens.push_back(buf);
      }
    }

    // std::cout << "============================================" << std::endl;
    // std::cout << "Box Line: " << *it << std::endl;
    // std::cout << "Hemisphere: " << _left << std::endl;

    /***************************************************/
    /*          START PROCESSING Y DIRECTION           */
    /***************************************************/
    double yStart1;

    double yStart2;
    double yEnd1;

    double yEnd2;
    int    yGridStartIndex(-1);

    int yGridEndIndex(-1);

    /* Determine y direction grid indices */
    std::string yStart = tokens[0].substr(0, 1);

    // std::cout << "YStart: " << yStart << std::endl;

    if (yStart.compare("E") == 0)
    {
      yStart = tokens[0].substr(0, 2);
    }

    if (yStart.compare("A") == 0)
    {
      yGridStartIndex = 0;
    }
    else if (yStart.compare("B") == 0)
    {
      yGridStartIndex = 1;
    }
    else if (yStart.compare("C") == 0)
    {
      yGridStartIndex = 2;
    }
    else if (yStart.compare("D") == 0)
    {
      yGridStartIndex = 3;
    }
    else if (yStart.compare("E1") == 0)
    {
      yGridStartIndex = 4;
    }
    else if (yStart.compare("E2") == 0)
    {
      yGridStartIndex = 5;
    }
    else if (yStart.compare("E3") == 0)
    {
      yGridStartIndex = 6;
    }
    else if (yStart.compare("F") == 0)
    {
      yGridStartIndex = 7;
    }
    else if (yStart.compare("G") == 0)
    {
      yGridStartIndex = 8;
    }
    else if (yStart.compare("H") == 0)
    {
      yGridStartIndex = 9;
    }
    else if (yStart.compare("I") == 0)
    {
      yGridStartIndex = 10;
    }

    /* Compute distance into desired ybox */
    std::string distance = tokens[0].substr(yStart.size());

    double boxDistancePercentage = 0.0;
    if (distance.empty())
    {
      boxDistancePercentage = 0.0;
    }
    else
    {
      boxDistancePercentage = std::stod(distance.c_str());
    }

    /* Base the Distance on the High Resolution Grid */
    yStart1 = this->TalairachGrid->GetPoint(yGridStartIndex * TALAIRACH_X_POINTS)[1];
    yStart2 = this->TalairachGrid->GetPoint((yGridStartIndex + 1) * TALAIRACH_X_POINTS)[1];
    vtkTalairachConversion::ImageType::PointType regionStart;
    regionStart.Fill(0.0);
    regionStart[1] = yStart1 + boxDistancePercentage * (yStart2 - yStart1);

    if ((this->SegmentationMode) && (yGridStartIndex == 0) && (boxDistancePercentage == 0.0))
    {
      yStart1 = this->TalairachGrid->GetPoint(yGridStartIndex * TALAIRACH_X_POINTS)[1];
      yStart2 = this->TalairachGrid->GetPoint((yGridStartIndex + 1) * TALAIRACH_X_POINTS)[1];
      regionStart[1] -= (yStart2 - yStart1);
    }

    /* Determine y direction grid indices */
    std::string yEnd = tokens[1].substr(0, 1);
    if (yEnd.compare("E") == 0)
    {
      yEnd = tokens[1].substr(0, 2);
    }

    // std::cout << "YEnd: " << yEnd << std::endl;

    if (yEnd.compare("A") == 0)
    {
      yGridEndIndex = 0;
    }
    else if (yEnd.compare("B") == 0)
    {
      yGridEndIndex = 1;
    }
    else if (yEnd.compare("C") == 0)
    {
      yGridEndIndex = 2;
    }
    else if (yEnd.compare("D") == 0)
    {
      yGridEndIndex = 3;
    }
    else if (yEnd.compare("E1") == 0)
    {
      yGridEndIndex = 4;
    }
    else if (yEnd.compare("E2") == 0)
    {
      yGridEndIndex = 5;
    }
    else if (yEnd.compare("E3") == 0)
    {
      yGridEndIndex = 6;
    }
    else if (yEnd.compare("F") == 0)
    {
      yGridEndIndex = 7;
    }
    else if (yEnd.compare("G") == 0)
    {
      yGridEndIndex = 8;
    }
    else if (yEnd.compare("H") == 0)
    {
      yGridEndIndex = 9;
    }
    else if (yEnd.compare("I") == 0)
    {
      yGridEndIndex = 10;
    }

    /* Compute distance into desired ybox */
    distance = tokens[1].substr(yEnd.size());

    if (distance.empty())
    {
      boxDistancePercentage = 1.0;
    }
    else
    {
      boxDistancePercentage = std::stod(distance.c_str());
    }

    /* Base the Distance on the High Resolution Grid */
    yEnd1 = this->TalairachGrid->GetPoint(yGridEndIndex * TALAIRACH_X_POINTS)[1];
    yEnd2 = this->TalairachGrid->GetPoint((yGridEndIndex + 1) * TALAIRACH_X_POINTS)[1];

    ImageType::PointType regionEnd;
    regionEnd.Fill(0.0);
    regionEnd[1] = yEnd1 + boxDistancePercentage * (yEnd2 - yEnd1);

    if ((this->SegmentationMode) && (yGridEndIndex == 10) && (boxDistancePercentage == 1.0))
    {
      yEnd1 = this->TalairachGrid->GetPoint(yGridEndIndex * TALAIRACH_X_POINTS)[1];
      yEnd2 = this->TalairachGrid->GetPoint((yGridEndIndex + 1) * TALAIRACH_X_POINTS)[1];
      regionEnd[1] += (yEnd2 - yEnd1);
    }

    /***************************************************/
    /*          START PROCESSING X DIRECTION           */
    /***************************************************/
    double xStart1;

    double xStart2;
    double xEnd1;

    double xEnd2;
    int    xGridStartIndex(-1);

    int xGridEndIndex(-1);

    std::string xStart = tokens[2].substr(0, 1);
    // std::cout << "XStart: " << xStart << std::endl;

    if (_left)
    {
      if (xStart.compare("a") == 0)
      {
        xGridStartIndex = 4;
      }
      else if (xStart.compare("b") == 0)
      {
        xGridStartIndex = 5;
      }
      else if (xStart.compare("c") == 0)
      {
        xGridStartIndex = 6;
      }
      else if (xStart.compare("d") == 0)
      {
        xGridStartIndex = 7;
      }

      distance = tokens[2].substr(xStart.size());
      if (distance.empty())
      {
        boxDistancePercentage = 0.0;
      }
      else
      {
        boxDistancePercentage = std::stod(distance.c_str());
      }

      /* Base the Distance on the High Resolution Grid */
      xStart1 = this->TalairachGrid->GetPoint(xGridStartIndex)[0];
      xStart2 = this->TalairachGrid->GetPoint(xGridStartIndex + 1)[0];
      regionStart[0] = xStart1 + boxDistancePercentage * (xStart2 - xStart1);
    }
    else
    {
      if (xStart.compare("a") == 0)
      {
        xGridEndIndex = 4;
      }
      else if (xStart.compare("b") == 0)
      {
        xGridEndIndex = 3;
      }
      else if (xStart.compare("c") == 0)
      {
        xGridEndIndex = 2;
      }
      else if (xStart.compare("d") == 0)
      {
        xGridEndIndex = 1;
      }

      distance = tokens[2].substr(xStart.size());
      if (distance.empty())
      {
        boxDistancePercentage = 0.0;
      }
      else
      {
        boxDistancePercentage = std::stod(distance.c_str());
      }

      /* Base the Distance on the High Resolution Grid */
      xEnd1 = this->TalairachGrid->GetPoint(xGridEndIndex)[0];
      xEnd2 = this->TalairachGrid->GetPoint(xGridEndIndex - 1)[0];
      regionEnd[0] = xEnd1 + boxDistancePercentage * (xEnd2 - xEnd1);
    }

    std::string xEnd = tokens[3].substr(0, 1);
    // std::cout << "xEnd: " << xEnd << std::endl;

    if (_left)
    {
      if (xEnd.compare("a") == 0)
      {
        xGridEndIndex = 4;
      }
      else if (xEnd.compare("b") == 0)
      {
        xGridEndIndex = 5;
      }
      else if (xEnd.compare("c") == 0)
      {
        xGridEndIndex = 6;
      }
      else if (xEnd.compare("d") == 0)
      {
        xGridEndIndex = 7;
      }

      distance = tokens[3].substr(xEnd.size());
      if (distance.empty())
      {
        boxDistancePercentage = 1.0;
      }
      else
      {
        boxDistancePercentage = std::stod(distance.c_str());
      }

      /* Base the Distance on the High Resolution Grid */
      xEnd1 = this->TalairachGrid->GetPoint(xGridEndIndex)[0];
      xEnd2 = this->TalairachGrid->GetPoint(xGridEndIndex + 1)[0];
      regionEnd[0] = xEnd1 + boxDistancePercentage * (xEnd2 - xEnd1);
      if ((this->SegmentationMode) && (xGridEndIndex == 7) && (boxDistancePercentage == 1.0))
      {
        // std::cout << "Expand X axis - END" << std::endl;
        regionEnd[0] += (xEnd2 - xEnd1);
      }
    }
    else
    {
      if (xEnd.compare("a") == 0)
      {
        xGridStartIndex = 4;
      }
      else if (xEnd.compare("b") == 0)
      {
        xGridStartIndex = 3;
      }
      else if (xEnd.compare("c") == 0)
      {
        xGridStartIndex = 2;
      }
      else if (xEnd.compare("d") == 0)
      {
        xGridStartIndex = 1;
      }

      distance = tokens[3].substr(xEnd.size());

      if (distance.empty())
      {
        boxDistancePercentage = 1.0;
      }
      else
      {
        boxDistancePercentage = std::stod(distance.c_str());
      }

      /* Base the Distance on the High Resolution Grid */
      xStart1 = this->TalairachGrid->GetPoint(xGridStartIndex)[0];
      xStart2 = this->TalairachGrid->GetPoint(xGridStartIndex - 1)[0];
      regionStart[0] = xStart1 + boxDistancePercentage * (xStart2 - xStart1);
      if ((this->SegmentationMode) && (xGridStartIndex == 1) && (boxDistancePercentage == 1.0))
      {
        // std::cout << "Expand X axis - Start" << std::endl;
        regionStart[0] += (xStart2 - xStart1);
      }
    }

    /***************************************************/
    /*          START PROCESSING Z DIRECTION           */
    /***************************************************/
    double zStart1;

    double zStart2;
    double zEnd1;

    double zEnd2;
    int    zGridStartIndex;

    int zGridEndIndex;

    size_t      pos = tokens[4].find(".");
    std::string zEnd = tokens[4].substr(0, pos);
    // std::cout << "zEnd: " << zEnd << std::endl;

    zGridEndIndex = TALAIRACH_Z_POINTS - static_cast<int>(floor(std::stod(zEnd.c_str())));

    if (pos > 0 && pos < 3)
    {
      distance = tokens[4].substr(pos);
    }
    else
    {
      distance = "";
    }

    if (distance.empty())
    {
      boxDistancePercentage = 0.0;
    }
    else
    {
      boxDistancePercentage = std::stod(distance.c_str());
    }

    // std::cout << "zEndIndex: " << zGridEndIndex << std::endl;

    /* Base the Distance on the High Resolution Grid */
    zEnd1 = this->TalairachGrid->GetPoint(zGridEndIndex * TALAIRACH_X_POINTS * TALAIRACH_Y_POINTS)[2];
    zEnd2 = this->TalairachGrid->GetPoint((zGridEndIndex - 1) * TALAIRACH_X_POINTS * TALAIRACH_Y_POINTS)[2];
    // std::cout << "zEnd1: " << zEnd1 << " " << zEnd2 << std::endl;
    regionEnd[2] = zEnd1 + boxDistancePercentage * (zEnd2 - zEnd1);
    if ((this->SegmentationMode) && (zGridEndIndex == 14) && (boxDistancePercentage == 0.0))
    {
      // std::cout << "Expand z axis - END" << std::endl;
      regionEnd[2] += (zEnd1 - zEnd2);
    }

    pos = tokens[5].find(".");
    std::string zStart = tokens[5].substr(0, pos);
    // std::cout << "zStart: " << zStart << std::endl;
    zGridStartIndex = TALAIRACH_Z_POINTS - static_cast<int>(floor(std::stod(zStart.c_str())));

    // std::cout << "zGridStartIndex: " << zGridStartIndex << std::endl;

    if (pos > 0 && pos < 3)
    {
      distance = tokens[5].substr(pos);
    }
    else
    {
      distance = "";
    }

    if (distance.empty())
    {
      boxDistancePercentage = 1.0;
    }
    else
    {
      boxDistancePercentage = std::stod(distance.c_str());
    }

    // std::cout << "boxDistancePercentage1: " << boxDistancePercentage <<
    // std::endl;
    // std::cout << "distance1: " << distance << std::endl;

    /* Base the Distance on the High Resolution Grid */
    zStart1 = this->TalairachGrid->GetPoint(zGridStartIndex * TALAIRACH_X_POINTS * TALAIRACH_Y_POINTS)[2];
    zStart2 = this->TalairachGrid->GetPoint((zGridStartIndex - 1) * TALAIRACH_X_POINTS * TALAIRACH_Y_POINTS)[2];
    // std::cout << "zStart1: " << zStart1 << " " << zStart2 << std::endl;
    regionStart[2] = zStart1 + boxDistancePercentage * (zStart2 - zStart1);
    if ((this->SegmentationMode) && (zGridStartIndex == 1) && (boxDistancePercentage == 1.0))
    {
      // std::cout << "Expand z axis - START" << std::endl;
      regionStart[2] += (zStart2 - zStart1);
    }

    // std::cout << "Start Location : " << regionStart << std::endl;
    // std::cout << "End Location : " << regionEnd << std::endl;

    ImageType::IndexType gridStart;
    ImageType::IndexType gridEnd;

    this->MaskImage->TransformPhysicalPointToIndex(regionStart, gridStart);
    this->MaskImage->TransformPhysicalPointToIndex(regionEnd, gridEnd);

    /* Make sure that the boxes are within the image space */
    ImageType::RegionType::SizeType imageSize = this->MaskImage->GetLargestPossibleRegion().GetSize();
    for (int i = 0; i < 3; i++)
    {
      if (gridStart[i] < 0)
      {
        gridStart[i] = 0;
        std::cout << "WARNING:" << std::endl;
        std::cout << "WARNING: Adjusting box bounds to fit inside of image. Please check your Talairach parameters."
                  << std::endl;
        std::cout << "WARNING:" << std::endl;
      }
      if (gridStart[i] >= static_cast<ImageType::IndexValueType>(imageSize[i]))
      {
        gridStart[i] = imageSize[i] - 1;
        std::cout << "WARNING:" << std::endl;
        std::cout << "WARNING: Adjusting box bounds to fit inside of image. Please check your Talairach parameters."
                  << std::endl;
        std::cout << "WARNING:" << std::endl;
      }
      if (gridEnd[i] < 0)
      {
        gridEnd[i] = 0;
        std::cout << "WARNING:" << std::endl;
        std::cout << "WARNING: Adjusting box bounds to fit inside of image. Please check your Talairach parameters."
                  << std::endl;
        std::cout << "WARNING:" << std::endl;
      }
      if (gridEnd[i] >= static_cast<ImageType::IndexValueType>(imageSize[i]))
      {
        gridEnd[i] = imageSize[i] - 1;
        std::cout << "WARNING:" << std::endl;
        std::cout << "WARNING: Adjusting box bounds to fit inside of image. Please check your Talairach parameters."
                  << std::endl;
        std::cout << "WARNING:" << std::endl;
      }
    }
    /* SOME POST PROCESSING TO ENSURE THAT THE BOXES ARE NOT EMPTY */
    for (int i = 0; i < 3; i++)
    {
      if (gridStart[i] > gridEnd[i])
      {
        long int tmp = gridEnd[i];
        gridEnd[i] = gridStart[i];
        gridStart[i] = tmp;
      }
    }

    // std::cout << "Start Index : " << gridStart << std::endl;
    // std::cout << "End Index : " << gridEnd << std::endl;

    ImageType::RegionType::SizeType gridSize;

    gridSize[0] = gridEnd[0] - gridStart[0];
    gridSize[1] = gridEnd[1] - gridStart[1];
    gridSize[2] = gridEnd[2] - gridStart[2];

    ImageType::RegionType testRegion;
    testRegion.SetSize(gridSize);
    testRegion.SetIndex(gridStart);

    // std::cout << "Region Size: " << gridSize << std::endl;
    // std::cout << "============================================" << std::endl;
    using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
    IteratorType itr(this->MaskImage, testRegion);
    for (itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
    {
      itr.Set(1);
    }
  }
}

void
vtkTalairachConversion::SetImageInformation(ImageType::Pointer exampleImage)
{
  this->MaskImage->SetOrigin(exampleImage->GetOrigin());
  this->MaskImage->SetSpacing(exampleImage->GetSpacing());
  this->MaskImage->SetDirection(exampleImage->GetDirection());
  this->MaskImage->SetRegions(exampleImage->GetLargestPossibleRegion());
}

vtkTalairachConversion::ImageType::Pointer
vtkTalairachConversion::GetImage()
{
  return this->MaskImage;
}

void
vtkTalairachConversion::Update()
{
  this->MaskImage->Allocate();
  this->MaskImage->FillBuffer(0);

  if (this->HemisphereMode == right || this->HemisphereMode == both)
  {
    ProcessBOX(false);
  }

  if (this->HemisphereMode == left || this->HemisphereMode == both)
  {
    ProcessBOX(true);
  }
}
