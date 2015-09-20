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
/*
 * Author: Hans J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "landmarkIO.h"
//#include "itk_hdf5.h"
//#include "itk_H5Cpp.h"
#include "itkNumberToString.h"

RGBImageType::Pointer ReturnOrientedRGBImage(SImageType::Pointer inputImage)
{
  RGBImageType::Pointer orientedImage;

  // std::cout << "inputImage information:\n" << inputImage << std::endl;

  // StatisticsImageFilter clears the buffer of the input image, so we have to pass a copy
  // of the input image through that
  SImageType::Pointer inputStatsImage;
    {
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(inputImage);
    duplicator->Update();
    inputStatsImage = duplicator->GetModifiableOutput();
    }

  itk::StatisticsImageFilter<SImageType>::Pointer stats = itk::StatisticsImageFilter<SImageType>::New();

  stats->SetInput(inputStatsImage);
  stats->Update();

  SImageType::PixelType minPixel( stats->GetMinimum() );
  SImageType::PixelType maxPixel( stats->GetMaximum() );

  // std::cout << "size of inputImage: " << inputImage->GetLargestPossibleRegion().GetSize()[0] << ","
  //    << inputImage->GetLargestPossibleRegion().GetSize()[1] << "," <<
  // inputImage->GetLargestPossibleRegion().GetSize()[2] << std::endl;

  //  itkUtil::WriteImage<SImageType>(inputImage, "inputImage.nii.gz");

  RGBImageType::Pointer rgbImage = RGBImageType::New();
  rgbImage->CopyInformation(inputImage);
  rgbImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  rgbImage->Allocate();

  // First just make RGB Image with greyscale values.
  itk::ImageRegionIterator<RGBImageType> rgbIt( rgbImage, rgbImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<SImageType>   sIt( inputImage, inputImage->GetLargestPossibleRegion() );
  for( ; !sIt.IsAtEnd(); ++rgbIt, ++sIt )
    {
    unsigned char charVal( ShortToUChar(sIt.Value(), minPixel, maxPixel) );
    RGBPixelType  pixel;
    pixel.SetRed(charVal);
    pixel.SetGreen(charVal);
    pixel.SetBlue(charVal);
    rgbIt.Set(pixel);
    }
  return orientedImage = rgbImage;
}

RGB2DImageType::Pointer GenerateRGB2DImage( RGBImageType::Pointer orientedImage)
{
  // Alocate 2DImage
  RGB2DImageType::Pointer   TwoDImage = RGB2DImageType::New();
  RGB2DImageType::IndexType TwoDIndex;

  TwoDIndex[0] = 0;
  TwoDIndex[1] = 0;
  RGB2DImageType::SizeType TwoDSize;
  TwoDSize[0] = orientedImage->GetLargestPossibleRegion().GetSize()[1];
  TwoDSize[1] = orientedImage->GetLargestPossibleRegion().GetSize()[2];
  RGB2DImageType::RegionType TwoDImageRegion;
  TwoDImageRegion.SetIndex(TwoDIndex);
  TwoDImageRegion.SetSize(TwoDSize);
  TwoDImage->SetRegions(TwoDImageRegion);
  TwoDImage->Allocate();

  // Fill 2DImage
  RGBImageType::IndexType ThreeDIndex;
  ThreeDIndex[0] = ( orientedImage->GetLargestPossibleRegion().GetSize()[0] ) / 2;
  for( TwoDIndex[1] = 0; TwoDIndex[1] < static_cast<signed int>( TwoDSize[1] ); ( TwoDIndex[1] )++ )
    {
    ThreeDIndex[2] = TwoDSize[1] - 1 - TwoDIndex[1];
    for( TwoDIndex[0] = 0; TwoDIndex[0] < static_cast<signed int>( TwoDSize[0] ); ( TwoDIndex[0] )++ )
      {
      ThreeDIndex[1] = TwoDIndex[0];
      TwoDImage->SetPixel( TwoDIndex, orientedImage->GetPixel(ThreeDIndex) );
      }
    }

  return TwoDImage;
}

static bool IsOnCylinder(const SImageType::PointType & curr_point,
                         const SImageType::PointType & center_point,
                         const SImageType::PointType & center_point2,
                         const double radius,
                         const double thickness)
{
  // const double cylinder_end=vcl_abs(height - vcl_abs(curr_point[0]-center_point[0]));
  const double APdist = curr_point[1] - center_point[1];
  const double ISdist = curr_point[2] - center_point[2];
  const double cylinder_side_squared =
    vcl_abs( radius * radius - ( APdist * APdist + ISdist * ISdist ) );
  const SImageType::PointType::VectorType PointDist2 =
    curr_point.GetVectorFromOrigin() - center_point2.GetVectorFromOrigin();

  return cylinder_side_squared < thickness * thickness
         || PointDist2.GetNorm() < 2;
}

static bool IsOnSphere(const SImageType::PointType & curr_point,
                       const SImageType::PointType & center_point,
                       const double radius)
{
  const SImageType::PointType::VectorType PointDist =
    curr_point.GetVectorFromOrigin() - center_point.GetVectorFromOrigin();

  return PointDist.GetNorm() < radius;
}

// TODO:  BrandedImages should be from the modelFile instead of the mDef.
void
MakeBrandeddebugImage(SImageType::ConstPointer in,
                      const landmarksConstellationModelIO & mDef,
                      const SImageType::PointType & RP,
                      const SImageType::PointType & AC,
                      const SImageType::PointType & PC,
                      const SImageType::PointType & VN4,
                      const std::string & fname,
                      const SImageType::PointType & RP2,
                      const SImageType::PointType & AC2,
                      const SImageType::PointType & PC2,
                      const SImageType::PointType & VN42)
{
  SImageType::Pointer inputImage = itkUtil::OrientImage<SImageType>(in,
                                                                    itk::SpatialOrientation::
                                                                    ITK_COORDINATE_ORIENTATION_RAI);

  RGBImageType::Pointer orientedImage = ReturnOrientedRGBImage( inputImage );

  for( unsigned int which = 0; which < 4; which++ )
    {
    SImageType::PointType pt = RP;
    SImageType::PointType pt2 = RP2;
    double                radius = 0.0;
    double                height = 0.0;

    switch( which )
      {
      case 0:
        {
        height = mDef.GetHeight("RP");
        radius = 4 * mDef.GetRadius("RP");
        pt = RP;
        pt2 = RP2;
        }
        break;
      case 1:
        {
        height = mDef.GetHeight("VN4");
        radius = 1.6 * mDef.GetRadius("VN4");
        pt = VN4;
        pt2 = VN42;
        }
        break;
      case 2:
        {
        height = mDef.GetHeight("AC");
        radius = 1.6 * mDef.GetRadius("AC");
        pt = AC;
        pt2 = AC2;
        }
        break;
      case 3:
        {
        height = mDef.GetHeight("PC");
        radius = 4 * mDef.GetRadius("PC");
        pt = PC;
        pt2 = PC2;
        }
        break;
      }

    itk::ImageRegionIterator<RGBImageType> rgbIt( orientedImage, orientedImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<SImageType>   sIt( inputImage, inputImage->GetLargestPossibleRegion() );
    for( ; !rgbIt.IsAtEnd() && !sIt.IsAtEnd(); ++sIt, ++rgbIt )
      {
      SImageType::IndexType index = rgbIt.GetIndex();
      SImageType::PointType p;
      orientedImage->TransformIndexToPhysicalPoint(index, p);

      if( IsOnCylinder(p, pt, pt2, radius, height) )
        {
        RGBPixelType pixel = rgbIt.Value();

        switch( which )
          {
          case 0:
            {
            pixel.SetRed(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 1:
            {
            pixel.SetGreen(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 2:
            {
            pixel.SetBlue(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 3:
            {
            pixel.SetRed(255);
            pixel.SetBlue(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          }
        }
      }
    }

  RGB2DImageType::Pointer TwoDImage = GenerateRGB2DImage( orientedImage );

  itkUtil::WriteImage<RGB2DImageType>(TwoDImage, fname);
  itkUtil::WriteImage<RGBImageType>(orientedImage, fname + "_ThreeDRGB.nii.gz");
  itkUtil::WriteImage<SImageType>(inputImage, fname + "_ThreeD.nii.gz");
}

void
MakePointBranded3DImage(SImageType::ConstPointer in,
                        const SImageType::PointType & CenterPoint,
                        const std::string & fname)
{
  SImageType::Pointer inputImage = itkUtil::OrientImage<SImageType>(in,
                                                                    itk::SpatialOrientation::
                                                                    ITK_COORDINATE_ORIENTATION_RAI);
  SImageType::Pointer inputStatsImage;
    {
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(inputImage);
    duplicator->Update();
    inputStatsImage = duplicator->GetModifiableOutput();
    }

  itk::StatisticsImageFilter<SImageType>::Pointer stats = itk::StatisticsImageFilter<SImageType>::New();

  stats->SetInput(inputStatsImage);
  stats->Update();
  SImageType::PixelType maxPixel( stats->GetMaximum() );

  SImageType::PointType pt = CenterPoint;
  double                radius = 3.0;

  itk::ImageRegionIterator<SImageType> sIt( inputImage, inputImage->GetLargestPossibleRegion() );
  for( ; !sIt.IsAtEnd(); ++sIt )
    {
    SImageType::IndexType index = sIt.GetIndex();
    SImageType::PointType p;
    inputImage->TransformIndexToPhysicalPoint(index, p);
    if( IsOnSphere(p, pt, radius) )
      {
      sIt.Set(maxPixel);
      }
    }

  itkUtil::WriteImage<SImageType>(inputImage, fname + "_Branded.nii.gz");
}

void
MakeBranded2DImage(SImageType::ConstPointer in,
                   landmarksConstellationDetector & myDetector,
                   const SImageType::PointType & RP,
                   const SImageType::PointType & AC,
                   const SImageType::PointType & PC,
                   const SImageType::PointType & VN4,
                   const SImageType::PointType & CM,
                   const std::string & fname)
{
  SImageType::Pointer inputImage = itkUtil::OrientImage<SImageType>(in,
                                                                    itk::SpatialOrientation::
                                                                    ITK_COORDINATE_ORIENTATION_RAI);

  RGBImageType::Pointer orientedImage = ReturnOrientedRGBImage( inputImage );

  for( unsigned int which = 0; which < 5; which++ )
    {
    SImageType::PointType pt;
    double                radius = 0.0;
    double                thickness = 2.0;

    switch( which )
      {
      case 0:
        {
        thickness = 1;
        radius = 1;
        pt = CM;
        }
        break;
      case 1:
        {
        radius = myDetector.GetModelRadius("AC");
        pt = AC;
        }
        break;
      case 2:
        {
        radius = myDetector.GetModelRadius("PC");
        pt = PC;
        }
        break;
      case 3:
        {
        radius = myDetector.GetModelRadius("RP");
        pt = RP;
        }
        break;
      case 4:
        {
        radius = myDetector.GetModelRadius("VN4");
        pt = VN4;
        }
        break;
      }

    itk::ImageRegionIterator<RGBImageType> rgbIt( orientedImage, orientedImage->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<SImageType>   sIt( inputImage, inputImage->GetLargestPossibleRegion() );
    for( ; !rgbIt.IsAtEnd() && !sIt.IsAtEnd(); ++sIt, ++rgbIt )
      {
      SImageType::IndexType index = rgbIt.GetIndex();
      SImageType::PointType p;
      orientedImage->TransformIndexToPhysicalPoint(index, p);
      if( IsOnCylinder(p, pt, pt, radius, thickness) )
        {
        RGBPixelType pixel = rgbIt.Value();

        switch( which )
          {
          case 0:
            {
            pixel.SetRed(255);
            pixel.SetGreen(255);
            pixel.SetBlue(0);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 1:
            {
            pixel.SetRed(0);
            pixel.SetGreen(255);
            pixel.SetBlue(0);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 2:
            {
            pixel.SetRed(0);
            pixel.SetGreen(0);
            pixel.SetBlue(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 3:
            {
            pixel.SetRed(255);
            pixel.SetGreen(0);
            pixel.SetBlue(0);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          case 4:
            {
            pixel.SetRed(255);
            pixel.SetGreen(0);
            pixel.SetBlue(255);
            rgbIt.Set(pixel);
            sIt.Set(255);
            }
            break;
          }
        }
      }
    }

  RGB2DImageType::Pointer TwoDImage = GenerateRGB2DImage( orientedImage );

  itkUtil::WriteImage<RGB2DImageType>(TwoDImage, fname);
}

// TODO:  Determine what the interface for WriteMRMLFile really needs to produce
// a useful file, and then limit the interface to just that.
extern void
WriteMRMLFile(std::string outputMRML,
              std::string outputLandmarksInInputSpace,
              std::string outputLandmarksInOutputSpace,
              std::string inputVolume,
              std::string outputVolume,
              std::string outputTransform,
              const LandmarksMapType & outputLandmarksInInputSpaceMap,
              const LandmarksMapType & outputLandmarksInOutputSpaceMap,
              VersorTransformType::ConstPointer versorTransform)
{
  const unsigned int LocalImageDimension = 3;
  itk::NumberToString<double>     doubleToString;

  typedef short                                      PixelType;
  typedef itk::Image<PixelType, LocalImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>            ReaderType;

  std::string mrmlFullFilename =
    itksys::SystemTools::CollapseFullPath( outputMRML.c_str() );
  std::string outputLandmarksInInputSpaceFullFilenameWithoutExtension = "";
  std::string outputLandmarksInOutputSpaceFullFilenameWithoutExtension = "";
  std::string inputVolumeFullFilename = "";
  std::string inputVolumeFilenameWithoutPath = "";
  std::string outputVolumeFullFilename = "";
  std::string outputVolumeFilenameWithoutPath = "";
  std::string outputTransformFullFilename = "";
  std::string outputTransformFilenameWithoutPath = "";

  if( outputLandmarksInInputSpace.compare("") != 0 )
    {
    outputLandmarksInInputSpace =
      itksys::SystemTools::CollapseFullPath( outputLandmarksInInputSpace.c_str() );
    outputLandmarksInInputSpaceFullFilenameWithoutExtension =
      itksys::SystemTools::GetFilenamePath( outputLandmarksInInputSpace ) + "/"
      + itksys::SystemTools::GetFilenameWithoutLastExtension( outputLandmarksInInputSpace);
    outputLandmarksInInputSpaceFullFilenameWithoutExtension =
      itksys::SystemTools::CollapseFullPath(
        outputLandmarksInInputSpaceFullFilenameWithoutExtension.c_str() );
    }

  if( outputLandmarksInOutputSpace.compare("") != 0 )
    {
    outputLandmarksInOutputSpace =
      itksys::SystemTools::CollapseFullPath( outputLandmarksInOutputSpace.c_str() );
    outputLandmarksInOutputSpaceFullFilenameWithoutExtension =
      itksys::SystemTools::GetFilenamePath( outputLandmarksInOutputSpace ) + "/"
      + itksys::SystemTools::GetFilenameWithoutLastExtension( outputLandmarksInOutputSpace );
    outputLandmarksInOutputSpaceFullFilenameWithoutExtension =
      itksys::SystemTools::CollapseFullPath( outputLandmarksInOutputSpaceFullFilenameWithoutExtension.c_str() );
    }

  if( inputVolume.compare("") != 0 )
    {
    inputVolumeFullFilename =
      itksys::SystemTools::CollapseFullPath( inputVolume.c_str() );
    inputVolumeFilenameWithoutPath =
      itksys::SystemTools::GetFilenameName( inputVolumeFullFilename );
    }

  if( outputVolume.compare("") != 0 )
    {
    outputVolumeFullFilename =
      itksys::SystemTools::CollapseFullPath( outputVolume.c_str() );
    outputVolumeFilenameWithoutPath =
      itksys::SystemTools::GetFilenameName( outputVolumeFullFilename );
    }

  if( outputTransform.compare("") != 0 )
    {
    outputTransformFullFilename =
      itksys::SystemTools::CollapseFullPath( outputTransform.c_str() );
    outputTransformFilenameWithoutPath =
      itksys::SystemTools::GetFilenameName( outputTransformFullFilename );
    }

  std::ofstream myfile( mrmlFullFilename.c_str() );
  if( !myfile.is_open() )
    {
    itkGenericExceptionMacro(<< "Cannot write mrml file!"
                             << mrmlFullFilename);
    }

  // Common mrml header
  myfile
    <<
    "<MRML  version=\"13298\" userTags=\"\">  \n  <Selection  \n  id=\"vtkMRMLSelectionNode1\"  name=\"vtkMRMLSelectionNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  activeVolumeID=\"vtkMRMLScalarVolumeNode2\"  secondaryVolumeID=\"NULL\"  activeLabelVolumeID=\"NULL\"  activeFiducialListID=\"vtkMRMLFiducialListNode1\"  activeROIListID=\"NULL\"  activeCameraID=\"NULL\"  activeViewID=\"NULL\"  activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>  \n  <Interaction  \n  id=\"vtkMRMLInteractionNode1\"  name=\"vtkMRMLInteractionNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  currentInteractionMode=\"ViewTransform\"  lastInteractionMode=\"ViewTransform\" ></Interaction>  \n  <Layout  \n  id=\"vtkMRMLLayoutNode1\"  name=\"vtkMRMLLayoutNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  currentViewArrangement=\"3\"  guiPanelVisibility=\"1\"  bottomPanelVisibility =\"0\"  guiPanelLR=\"0\"  collapseSliceControllers=\"0\"  \n  numberOfCompareViewRows=\"1\"  numberOfCompareViewColumns=\"1\"  numberOfLightboxRows=\"1\"  numberOfLightboxColumns=\"1\"  mainPanelSize=\"400\"  secondaryPanelSize=\"400\" ></Layout>  \n  <TGParameters  \n  id=\"vtkMRMLChangeTrackerNode1\"  name=\"vtkMRMLChangeTrackerNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  ROIMin=\"-1 -1 -1\"  ROIMax=\"-1 -1 -1\"  SegmentThresholdMin=\"-1\"  SegmentThresholdMax=\"-1\"  Analysis_Intensity_Flag=\"0\"  Analysis_Deformable_Flag=\"0\"  UseITK=\"1\"  RegistrationChoice=\"3\"  ROIRegistration=\"1\"  ResampleChoice=\"3\"  ResampleConst=\"0.5\" ></TGParameters>  \n  <Crosshair  \n  id=\"vtkMRMLCrosshairNode1\"  name=\"vtkMRMLCrosshairNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  crosshairMode=\"NoCrosshair\"  navigation=\"true\"  crosshairBehavior=\"Normal\"  crosshairThickness=\"Fine\"  crosshairRAS=\"0 0 0\" ></Crosshair>  \n  <Slice  \n  id=\"vtkMRMLSliceNode1\"  name=\"Green\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"374.356 322.12 0.534398\"  dimensions=\"430 370 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"-1 0 0 2.04504 0 0 1 0 0 1 0 -9.66354 0 0 0 1\"  layoutName=\"Green\"  orientation=\"Coronal\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>  \n  <SliceComposite  \n  id=\"vtkMRMLSliceCompositeNode1\"  name=\"vtkMRMLSliceCompositeNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Green\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>  \n  <Slice  \n  id=\"vtkMRMLSliceNode2\"  name=\"Red\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"364.349 313.51 1.29601\"  dimensions=\"430 370 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"-1 0 0 2.04504 0 1 0 -12.7088 0 0 1 0 0 0 0 1\"  layoutName=\"Red\"  orientation=\"Axial\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>  \n  <SliceComposite  \n  id=\"vtkMRMLSliceCompositeNode2\"  name=\"vtkMRMLSliceCompositeNode2\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Red\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>  \n  <Slice  \n  id=\"vtkMRMLSliceNode3\"  name=\"Yellow\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fieldOfView=\"374.356 322.12 1.11022\"  dimensions=\"430 370 1\"  activeSlice=\"0\"  layoutGridRows=\"1\"  layoutGridColumns=\"1\"  sliceToRAS=\"0 0 1 0 -1 0 0 -12.7088 0 1 0 -9.66354 0 0 0 1\"  layoutName=\"Yellow\"  orientation=\"Sagittal\"  jumpMode=\"1\"  sliceVisibility=\"false\"  widgetVisibility=\"false\"  useLabelOutline=\"false\"  sliceSpacingMode=\"0\"  prescribedSliceSpacing=\"1 1 1\" ></Slice>  \n  <SliceComposite  \n  id=\"vtkMRMLSliceCompositeNode3\"  name=\"vtkMRMLSliceCompositeNode3\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  backgroundVolumeID=\"vtkMRMLScalarVolumeNode1\"  foregroundVolumeID=\"\"  labelVolumeID=\"\"  compositing=\"0\"  labelOpacity=\"1\"  linkedControl=\"0\"  foregroundGrid=\"0\"  backgroundGrid=\"0\"  labelGrid=\"1\"  fiducialVisibility=\"1\"  fiducialLabelVisibility=\"1\"  sliceIntersectionVisibility=\"0\"  layoutName=\"Yellow\"  annotationMode=\"All\"  doPropagateVolumeSelection=\"1\" ></SliceComposite>  \n  <ScriptedModule  \n  id=\"vtkMRMLScriptedModuleNode1\"  name=\"vtkMRMLScriptedModuleNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\" ></ScriptedModule>  \n  <View  \n  id=\"vtkMRMLViewNode1\"  name=\"View6\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  active=\"true\"  visibility=\"true\"  fieldOfView=\"200\"  letterSize=\"0.05\"  boxVisible=\"true\"  fiducialsVisible=\"true\"  fiducialLabelsVisible=\"true\"  axisLabelsVisible=\"true\"  backgroundColor=\"0.70196 0.70196 0.90588\"  animationMode=\"Off\"  viewAxisMode=\"LookFrom\"  spinDegrees=\"2\"  spinMs=\"5\"  spinDirection=\"YawLeft\"  rotateDegrees=\"5\"  rockLength=\"200\"  rockCount=\"0\"  stereoType=\"NoStereo\"  renderMode=\"Perspective\" ></View>  \n  <Camera  \n  id=\"vtkMRMLCameraNode1\"  name=\"Camera6\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  position=\"0 500 0\"  focalPoint=\"0 0 0\"  viewUp=\"0 0 1\"  parallelProjection=\"false\"  parallelScale=\"1\"  activetag=\"vtkMRMLViewNode1\" ></Camera>  \n";

  // For fiducial landmarks in output space
  if( outputLandmarksInOutputSpace.compare("") != 0 )
    {
    myfile << "<FiducialList\n id=\"vtkMRMLFiducialListNode1\" name=\""
           << outputLandmarksInOutputSpaceFullFilenameWithoutExtension
           <<
      "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode1\" userTags=\"\" symbolScale=\"5\" symbolType=\"13\" textScale=\"4.5\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"1 0.5 0.5\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\"\n";

    LandmarksMapType::const_iterator it;
    unsigned int                     index = 0;
    for( it = outputLandmarksInOutputSpaceMap.begin(); it != outputLandmarksInOutputSpaceMap.end(); ++it )
      {
      myfile << "id " << it->first << " labeltext " << it->first << " xyz "
             << doubleToString( ( it->second )[0]) << " "
             << doubleToString( ( it->second )[1]) << " "
             << doubleToString( ( it->second )[2]);
      if( ++index < outputLandmarksInOutputSpaceMap.size() )
        {
        myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\n";
        }
      else
        {
        myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>\n";
        }
      }

    myfile
      <<
      "<FiducialListStorage\n id=\"vtkMRMLFiducialListStorageNode1\" name=\"vtkMRMLFiducialListStorageNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""
      << outputLandmarksInOutputSpace
      << "\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></FiducialListStorage>\n";
    }

  // For fiducial landmarks in input space
  if( outputLandmarksInInputSpace.compare("") != 0 )
    {
    myfile << "<FiducialList\n id=\"vtkMRMLFiducialListNode2\" name=\""
           << outputLandmarksInInputSpaceFullFilenameWithoutExtension
           <<
      "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode2\" userTags=\"\" symbolScale=\"5\" symbolType=\"13\" textScale=\"4.5\" visibility=\"1\" color=\"0.4 1 1\" selectedcolor=\"1 0.5 0.5\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\"\n";

    LandmarksMapType::const_iterator it;
    unsigned int                     index = 0;
    for( it = outputLandmarksInInputSpaceMap.begin(); it != outputLandmarksInInputSpaceMap.end(); ++it )
      {
      myfile << "id " << it->first << " labeltext " << it->first << " xyz "
             << doubleToString( ( it->second )[0]) << " "
             << doubleToString( ( it->second )[1]) << " "
             << doubleToString( ( it->second )[2]);
      if( ++index < outputLandmarksInInputSpaceMap.size() )
        {
        myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\n";
        }
      else
        {
        myfile << " orientationwxyz 0 0 0 1 selected 1 visibility 1\"></FiducialList>\n";
        }
      }

    myfile
      <<
      "<FiducialListStorage\n id=\"vtkMRMLFiducialListStorageNode2\" name=\"vtkMRMLFiducialListStorageNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""
      << outputLandmarksInInputSpace
      << "\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></FiducialListStorage>\n";
    }

  // For output volume
  if( outputVolume.compare( "" ) != 0 )
    {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( outputVolumeFullFilename );
    reader->Update();
    ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
    ImageType::SpacingType   spacing = reader->GetOutput()->GetSpacing();
    ImageType::PointType     origin = reader->GetOutput()->GetOrigin();
    myfile
      <<
      "<VolumeArchetypeStorage\n id=\"vtkMRMLVolumeArchetypeStorageNode1\" name=\"vtkMRMLVolumeArchetypeStorageNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""
      << outputVolumeFullFilename
      <<
      "\" useCompression=\"1\" readState=\"0\" writeState=\"0\" centerImage=\"0\" singleFile=\"0\" UseOrientationFromFile=\"1\"></VolumeArchetypeStorage>\n";

    myfile << "<Volume\n id=\"vtkMRMLScalarVolumeNode1\" name=\"" << outputVolumeFilenameWithoutPath
           <<
      "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLVolumeArchetypeStorageNode1\" userTags=\"\" displayNodeRef=\"vtkMRMLScalarVolumeDisplayNode1\" ijkToRASDirections=\"";
    for( unsigned int i = 0; i < LocalImageDimension; ++i )
      {
      for( unsigned int j = 0; j < LocalImageDimension; ++j )
        {
        myfile << doubleToString(direction(i, j) ) << " ";
        }
      }
    myfile
      << "\" spacing=\""
      << doubleToString(spacing[0]) << " "
      << doubleToString(spacing[1]) << " "
      << doubleToString(spacing[2])
      << "\" origin=\""
      << doubleToString(origin[0]) << " "
      << doubleToString(origin[1]) << " "
      << doubleToString(origin[2])
      <<
      "\" labelMap=\"0\"></Volume>\n<VolumeDisplay\n id=\"vtkMRMLScalarVolumeDisplayNode1\" name=\"vtkMRMLScalarVolumeDisplayNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeGrey\"  window=\"204\" level=\"153\" upperThreshold=\"32767\" lowerThreshold=\"-32768\" interpolate=\"1\" autoWindowLevel=\"1\" applyThreshold=\"0\" autoThreshold=\"0\"></VolumeDisplay>\n";
    }

  // For input volume
  if( inputVolume.compare( "" ) != 0 )
    {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputVolumeFullFilename );
    reader->Update();
    ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
    ImageType::SpacingType   spacing = reader->GetOutput()->GetSpacing();
    ImageType::PointType     origin = reader->GetOutput()->GetOrigin();
    myfile
      <<
      "<VolumeArchetypeStorage\n id=\"vtkMRMLVolumeArchetypeStorageNode2\" name=\"vtkMRMLVolumeArchetypeStorageNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""
      << inputVolumeFullFilename
      <<
      "\" useCompression=\"1\" readState=\"0\" writeState=\"0\" centerImage=\"0\" singleFile=\"0\" UseOrientationFromFile=\"1\"></VolumeArchetypeStorage>\n";

    myfile << "<Volume\n id=\"vtkMRMLScalarVolumeNode2\" name=\"" << inputVolumeFilenameWithoutPath
           <<
      "\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLVolumeArchetypeStorageNode2\" userTags=\"\" displayNodeRef=\"vtkMRMLScalarVolumeDisplayNode2\" ijkToRASDirections=\"";
    for( unsigned int i = 0; i < LocalImageDimension; ++i )
      {
      for( unsigned int j = 0; j < LocalImageDimension; ++j )
        {
        myfile << doubleToString(direction(i, j) ) << " ";
        }
      }
    myfile
      << "\" spacing=\""
      << doubleToString(spacing[0]) << " "
      << doubleToString(spacing[1]) << " "
      << doubleToString(spacing[2])
      << "\" origin=\""
      << doubleToString(origin[0]) << " "
      << doubleToString(origin[1]) << " "
      << doubleToString(origin[2])
      <<
      "\" labelMap=\"0\"></Volume>\n<VolumeDisplay\n id=\"vtkMRMLScalarVolumeDisplayNode2\" name=\"vtkMRMLScalarVolumeDisplayNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"0.5 0.5 0.5\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0\" diffuse=\"1\" selectedSpecular=\"0.5\" specular=\"0\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeGrey\"  window=\"204\" level=\"153\" upperThreshold=\"32767\" lowerThreshold=\"-32768\" interpolate=\"1\" autoWindowLevel=\"1\" applyThreshold=\"0\" autoThreshold=\"0\"></VolumeDisplay>\n";
    }

  // For output transform
  if( outputTransform.compare("") != 0 )
    {
    VersorTransformMatrixType tm = versorTransform->GetMatrix();

    myfile
      <<
      "<TransformStorage  \n  id=\"vtkMRMLTransformStorageNode1\"  name=\"vtkMRMLTransformStorageNode1\"  hideFromEditors=\"true\"  selectable=\"true\"  selected=\"false\"  fileName=\""
      << outputTransformFullFilename
      << "\"  useCompression=\"1\"  readState=\"0\"  writeState=\"0\" ></TransformStorage>  \n"
      <<
      "<LinearTransform  \n  id=\"vtkMRMLLinearTransformNode1\"  name=\""
      << outputTransformFilenameWithoutPath
      <<
      "\"  hideFromEditors=\"false\"  selectable=\"true\"  selected=\"false\"  storageNodeRef=\"vtkMRMLTransformStorageNode1\"  userTags=\"\"  matrixTransformToParent=\"";
    for( unsigned int i = 0; i <  LocalImageDimension; ++i )
      {
      for( unsigned int j = 0; j <  LocalImageDimension; ++j )
        {
        myfile << doubleToString(tm(i, j) ) << " ";
        }
      }
    myfile << "\" ></LinearTransform>  \n";
    }

  myfile << "</MRML>\n";
  myfile.close();
}

void
loadLLSModel(std::string llsModelFilename,
             std::map<std::string, std::vector<double> > & llsMeans,
             std::map<std::string, MatrixType> & llsMatrices,
             std::map<std::string, double> & searchRadii)
{
  std::ifstream myfile( llsModelFilename.c_str() );

  if( !myfile.is_open() )
    {
    itkGenericExceptionMacro(<< "Cannot open landmark model file!"
                             << llsModelFilename);
    }

  // for each landmark
  std::string line;
  while( getline(myfile, line) )
    {
    // skip newline or comments between landmarks
    if( ( line.compare(0, 1, "#") != 0 )
        && ( line.compare(0, 1, "\0") != 0 ) )
      {
      unsigned int dimension = 3; // for 3D parameters

      // read in the landmark name
      std::string name = line;

      if( !getline(myfile, line) )
        {
        itkGenericExceptionMacro(<< "Bad number of parameters info in llsModelFile!");
        }
      else
        {
        // read in mean values associated with PCA model
        unsigned int pos1 = 0;
        unsigned int pos2 = line.find(' ');
        unsigned int i = 0;
        while( pos2 < line.size() )
          {
          llsMeans[name].push_back(
            atof( line.substr(pos1, pos2 - pos1).c_str() ) );
          ++i;
          pos1 = pos2 + 1;
          pos2 = line.find(' ', pos1 + 1);
          }

        if( i != dimension )  // double check
          {
          itkGenericExceptionMacro(<< "Bad mean values in llsModelFile!");
          }
        }

      // read in search radius
      if( !getline(myfile, line) )
        {
        itkGenericExceptionMacro(<< "Bad search radius in llsModelFile!");
        }
      else
        {
        searchRadii[name] = atof( line.c_str() );
        }

      // read in the number of linear model coefficients
      unsigned int numParameters = 0;
      if( !getline(myfile, line) )
        {
        itkGenericExceptionMacro(<< "Bad number of parameters info in llsModelFile!");
        }
      else
        {
        numParameters = atoi( line.c_str() );
        }

      MatrixType coefficients; // linear model coefficients
      coefficients.set_size(dimension, numParameters);
      for( unsigned int j = 0; j < dimension; ++j )
        {
        if( !getline(myfile, line) )
          {
          itkGenericExceptionMacro(<< "Bad linear model coefficients in llsModelFile!")
          }
        else
          {
          unsigned int pos1 = 0;
          unsigned int pos2 = line.find(' ');
          unsigned int i = 0;
          while( pos2 < line.size() )
            {
            coefficients(j, i++) =
              atof( line.substr(pos1, pos2 - pos1).c_str() );
            pos1 = pos2 + 1;
            pos2 = line.find(' ', pos1 + 1);
            }

          if( i != numParameters )  // double check
            {
            itkGenericExceptionMacro(<< "Bad linear model coefficients in llsModelFile!");
            }
          }
        }
      llsMatrices[name] = coefficients;
      }
    }

  myfile.close();
}

void
writeVerificationScript(std::string outputVerificationScriptFilename,
                        std::string outputVolume,
                        std::string saveOutputLandmarksFilename)
{
  std::ofstream ScriptFile;

  ScriptFile.open( outputVerificationScriptFilename.c_str() ); // open setup
                                                               // file for
                                                               // writing
  if( !ScriptFile.is_open() )
    {
    itkGenericExceptionMacro(<< "Can't write outputVerificationScript " << outputVerificationScriptFilename);
    }

  ScriptFile << "##There is a program that reads in a T1 image, determine the AC, PC, VN4, and MPJ" << std::endl;
  ScriptFile << "##points and writes out a slicer .fcsv file with 3 landmark points in it" << std::endl;
  ScriptFile << "##There are 2 files to be read in order to confirm that the landmarks are" << std::endl;
  ScriptFile << "##being picked correctly." << std::endl;
  ScriptFile << "set volumesLogic [$::slicer3::VolumesGUI GetLogic]" << std::endl;
  ScriptFile << "set fidLogic [$::slicer3::FiducialsGUI GetLogic]" << std::endl;
  ScriptFile << std::endl;
  ScriptFile << "set fileNameT1  " << outputVolume << std::endl;
  ScriptFile << "set  landmarkName " << saveOutputLandmarksFilename << std::endl;
  ScriptFile << std::endl;
  ScriptFile << "set volumeNode [$volumesLogic AddArchetypeVolume $fileNameT1 T1brain 0]" << std::endl;
  ScriptFile << "set selectionNode [$::slicer3::ApplicationLogic GetSelectionNode]" << std::endl;
  ScriptFile << "$selectionNode SetReferenceActiveVolumeID [$volumeNode GetID]" << std::endl;
  ScriptFile << "$::slicer3::ApplicationLogic PropagateVolumeSelection" << std::endl;
  ScriptFile << std::endl;
  ScriptFile << "$fidLogic LoadFiducialList $landmarkName" << std::endl;
  ScriptFile << std::endl;
  ScriptFile.close();
}
