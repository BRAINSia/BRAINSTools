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
#include "BRAINSABCUtilities.h"
/*****************************
 * Now call the instantiations
 */
#include "BRAINSABCUtilities.hxx"
#include "LLSBiasCorrector.h"

template std::vector<FloatImageType::Pointer> DuplicateImageList<FloatImageType>(
  const std::vector<FloatImageType::Pointer> & );

template std::vector<ShortImageType::Pointer> DuplicateImageList<ShortImageType>(
  const std::vector<ShortImageType::Pointer> & );

template void NormalizeProbListInPlace<FloatImageType>(std::vector<FloatImageType::Pointer> & );

template void ZeroNegativeValuesInPlace<FloatImageType>(  std::vector<FloatImageType::Pointer> & );

MapOfFloatImageVectors
ResampleImageListToFirstKeyImage(const std::string & resamplerInterpolatorType,
                                 const MapOfFloatImageVectors & inputImageMap)
{
  muLogMacro(<< "Resampling input image map to the first key image." << std::endl);

  FloatImageType::ConstPointer KeyImageFirstRead =
    GetMapVectorFirstElement(inputImageMap).GetPointer();

  // Clear image list
  MapOfFloatImageVectors outputImageMap;

  // Resample the other images
  for(const auto & elem : inputImageMap)
    {
    auto currImageIter = elem.second.begin();
    unsigned int i(0);
    while( currImageIter != elem.second.end() )
      {
      FloatImageType::Pointer tmp =
        ResampleImageWithIdentityTransform<FloatImageType>( resamplerInterpolatorType,
                                                            0,
                                                            (*currImageIter).GetPointer(),
                                                            KeyImageFirstRead.GetPointer() );
      // Add the image
      outputImageMap[elem.first].push_back( tmp );
      ++currImageIter;
      ++i;
      }
    }
  return outputImageMap;
}

MapOfFloatImageVectors
ResampleInPlaceImageList(const std::string & resamplerInterpolatorType,
                         const MapOfFloatImageVectors & inputImageMap,
                         MapOfTransformLists & intraSubjectTransforms)
{
  muLogMacro(<< "ResampleInPlaceImageList..." << std::endl);
  /*
   * This function, first, transforms all inputImageMap to the space of the first image of the map
   * using rigid transforms (intraSubjectTransforms) and Resampling InPlace interoplation.
   * Then, it resamples all images within one modality to the voxel lattice of the fist image of that modality channel
   * using resamplerInterpolatorType and Identity transform.
   */

  MapOfFloatImageVectors resampleInPlaceImageMap;
  MapOfFloatImageVectors outputImageMap;

  PrintMapOfImageVectors(inputImageMap);

  // ResampleInPlace all images to the physical space of the first image
  //
  for(auto inputImageMapIter = inputImageMap.begin();
      inputImageMapIter != inputImageMap.end(); ++inputImageMapIter)
    {
    auto currModalIter = inputImageMapIter->second.begin();
    unsigned int i(0);
    auto xfrmIt = intraSubjectTransforms[inputImageMapIter->first].begin();
    while( currModalIter != inputImageMapIter->second.end() )
      {
      muLogMacro(<< "ResamplingInPlace input image " << inputImageMapIter->first << " #" << i
                 << " to the physical space of the first image." << std::endl);
      typedef itk::ResampleInPlaceImageFilter<FloatImageType, FloatImageType>  ResampleIPFilterType;
      typedef ResampleIPFilterType::Pointer                                    ResampleIPFilterPointer;

      typedef itk::VersorRigid3DTransform<double>   VersorRigid3DTransformType;
      const VersorRigid3DTransformType::ConstPointer tempRigidTransform =
        dynamic_cast<VersorRigid3DTransformType const *>( (*xfrmIt).GetPointer() );
      if( tempRigidTransform.IsNull() )
        {
        std::cerr << "Error in type conversion. " << __FILE__ << __LINE__ << std::endl;
        std::cerr << "ResampleInPlace is only allowed with rigid transform type." << std::endl;
        throw;
        }

      ResampleIPFilterPointer resampleIPFilter = ResampleIPFilterType::New();
      resampleIPFilter->SetInputImage( (*currModalIter) );
      resampleIPFilter->SetRigidTransform( tempRigidTransform );
      resampleIPFilter->Update();
      FloatImageType::Pointer tmp = resampleIPFilter->GetOutput();

      // Add the image
      resampleInPlaceImageMap[inputImageMapIter->first].push_back(tmp);
      ++currModalIter;
      ++xfrmIt;
      ++i;
      }
    }

  // Resample each intra subject image to the first image of its modality
  //
  muLogMacro(<< "Resampling each intra subject image to the first image of its modality using "
             << resamplerInterpolatorType << " interpolation." << std::endl);
  const FloatImageType::PixelType outsideFOVCode = vnl_huge_val( static_cast<FloatImageType::PixelType>( 1.0f ) );

  for(auto & elem : resampleInPlaceImageMap)
    {
    auto currModalIter = elem.second.begin();
    FloatImageType::Pointer currModalityKeySubjectImage = (*currModalIter).GetPointer(); // the first image of current modality

    while( currModalIter != elem.second.end() )
      {
      if( (*currModalIter).GetPointer() == currModalityKeySubjectImage.GetPointer() ) // if the current intra subject image
          // is the first image in current modality list
        {
        outputImageMap[elem.first].push_back( (*currModalIter).GetPointer() );
        }
      else // resample the current intra subject image to the first image of the current modality
        {
        FloatImageType::Pointer tmp =
          ResampleImageWithIdentityTransform<FloatImageType>( resamplerInterpolatorType,
                                                             outsideFOVCode,
                                                             (*currModalIter).GetPointer(),
                                                             currModalityKeySubjectImage.GetPointer() );

        // Zero the mask region outside FOV and also the intensities with
        // outside
        // FOV code
        typedef itk::ImageRegionIterator<FloatImageType> InternalIteratorType;
        InternalIteratorType    tmpIt( tmp, tmp->GetLargestPossibleRegion() );

        //TODO:  This code below with masking does not make sense.
        //        intraSubjectFOVIntersectionMask does not seem to do anything.
        // HACK:  We can probably remove the mask generation from here.
        // The FOV mask, regions where intensities in all channels do not
        // match FOV code
        ByteImageType::Pointer intraSubjectFOVIntersectionMask = ByteImageType::New();
        intraSubjectFOVIntersectionMask->CopyInformation( currModalityKeySubjectImage.GetPointer() );
        intraSubjectFOVIntersectionMask->SetRegions( currModalityKeySubjectImage.GetPointer()->GetLargestPossibleRegion() );
        intraSubjectFOVIntersectionMask->Allocate();
        intraSubjectFOVIntersectionMask->FillBuffer(1);

        typedef itk::ImageRegionIterator<ByteImageType> MaskIteratorType;
        MaskIteratorType maskIt( intraSubjectFOVIntersectionMask,
                                intraSubjectFOVIntersectionMask->GetLargestPossibleRegion() );
        maskIt.GoToBegin();
        tmpIt.GoToBegin();
        while( !maskIt.IsAtEnd() )
          {
          if( tmpIt.Get() == outsideFOVCode )  // Voxel came from outside
            // the original FOV during
            // registration, so
            // invalidate it.
            {
            maskIt.Set(0); // Set it as an invalid voxel in
            // intraSubjectFOVIntersectionMask
            tmpIt.Set(0);  // Set image intensity value to zero.
            }
          ++maskIt;
          ++tmpIt;
          }

        // Add the image
        outputImageMap[elem.first].push_back(tmp);
        }
      ++currModalIter;
      }
    }

  return outputImageMap;
}
