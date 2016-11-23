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
#ifndef __PhilipsDWIConverter_h
#define __PhilipsDWIConverter_h
#include "DWIConverter.h"
#include "itkExtractImageFilter.h"

/** specific converter for Philips scanners */
class PhilipsDWIConverter : public DWIConverter
{
public:
  PhilipsDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections) : DWIConverter(allHeaders,inputFileNames,
                                                                        useBMatrixGradientDirections)
    {
    }

  virtual ~PhilipsDWIConverter() {}

  virtual void LoadDicomDirectory() ITK_OVERRIDE
    {
      this->DWIConverter::LoadDicomDirectory();
      if(!this->m_MultiSliceVolume)
        {
        this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
        this->m_MeasurementFrame = this->m_Volume->GetDirection();
        this->DetermineSliceOrderIS();
        this->SetDirectionsFromSliceOrder();
        }
      // single-frame file handled specially
    }
  void ExtractDWIData() ITK_OVERRIDE
    {
      if( !this->m_MultiSliceVolume )
        {
        // assume volume interleaving
        std::cout << "Number of Slices: " << this->m_NSlice << std::endl;
        std::cout << "Number of Volumes: " << this->m_NVolume << std::endl;
        std::cout << "Number of Slices in each volume: " << this->m_SlicesPerVolume << std::endl;
        // NOTE:  Philips interleaves the directions, so the all gradient directions can be
        // determined in the first "nVolume" slices which represents the first slice from each
        // of the gradient volumes.
        for( unsigned int k = 0; k < this->m_NVolume; ++k )
          {
          std::string DiffusionDirectionality;
          bool        useSupplement49Definitions(false);
          if( this->m_Headers[k]->GetElementCSorOB(0x0018, 0x9075, DiffusionDirectionality, false) == EXIT_SUCCESS )
            {
            useSupplement49Definitions = true;
            }

          bool   B0FieldFound = false;
          double b = 0.0;
          const double zeroBValueTolerance = 0.1;  // Implausibly small value assumed to be b0 images
            {// Fill out the B-Values
            if( useSupplement49Definitions == true )
              {
              B0FieldFound = this->m_Headers[k]->GetElementFD(0x0018, 0x9087, b, false) == EXIT_SUCCESS;
              }
            else
              {
              float floatB;
              if( this->m_Headers[k]->GetElementFLorOB(0x2001, 0x1003, floatB, false) == EXIT_SUCCESS )
                {
                B0FieldFound = true;
                }
              if( B0FieldFound )
                {
                b = static_cast<double>(floatB);
                }
              std::string tag;
              this->m_Headers[k]->GetElementCSorOB(0x2001, 0x1004, tag, false );
              if( StringContains(tag, "I") && ( b > zeroBValueTolerance ) )
                {
                DiffusionDirectionality = "ISOTROPIC";
                }
              }
            }

          vnl_vector_fixed<double, 3> vect3d;
          vect3d.fill( 0.0 );
            {// Try to fill in the vect3d gradient directions
            if( useSupplement49Definitions == true )
              {
              double doubleArray[3];
              // Use alternate method to get value out of a sequence header (Some Phillips Data).
              if(this->m_Headers[k]->GetElementFD(0x0018, 0x9089, 3, doubleArray, false) == EXIT_FAILURE)
                {
                //Try alternate method.
                // std::cout << "Looking for  0018|9089 in sequence 0018,9076" << std::endl;
                // gdcm::SeqEntry *
                // DiffusionSeqEntry=this->m_Headers[k]->GetSeqEntry(0x0018,0x9076);
                itk::DCMTKSequence DiffusionSeqEntry;
                this->m_Headers[k]->GetElementSQ(0x0018, 0x9076, DiffusionSeqEntry);
                // const unsigned int
                // n=DiffusionSeqEntry->GetNumberOfSQItems();
                unsigned int n = DiffusionSeqEntry.card();
                if( n == 0 )
                  {
                  std::cout << "ERROR:  Sequence entry 0018|9076 has no items." << std::endl;
                  throw;
                  }
                DiffusionSeqEntry.GetElementFD(0x0018, 0x9089, 3, doubleArray);
                }
              vect3d[0] = doubleArray[0];
              vect3d[1] = doubleArray[1];
              vect3d[2] = doubleArray[2];
              std::cout << "===== gradient orientations:" << k << " "
                << this->m_InputFileNames[k] << " (0018,9089) " << " " << vect3d << std::endl;
              }
            else
              {
              float tmp[3];
              itk::DCMTKFileReader *hdr = this->m_Headers[k];
              if(hdr->GetElementFLorOB( 0x2005, 0x10b0, tmp[0],false) == EXIT_FAILURE ||
                 hdr->GetElementFLorOB( 0x2005, 0x10b1, tmp[1],false) == EXIT_FAILURE ||
                 hdr->GetElementFLorOB( 0x2005, 0x10b2, tmp[2],false) == EXIT_FAILURE)
                {
                // try with new philips tags.
                if(hdr->GetElementFLorOB( 0x2005, 0x12b0, tmp[0],false) == EXIT_FAILURE ||
                   hdr->GetElementFLorOB( 0x2005, 0x12b1, tmp[1],false) == EXIT_FAILURE ||
                   hdr->GetElementFLorOB( 0x2005, 0x12b2, tmp[2],false) == EXIT_FAILURE)
                  {
                  // last chance. GetELementFD throws exception on failure
                  double tmpD[3];
                  hdr->GetElementFD(0x0018,0x9089,3,tmpD);
                  tmp[0] = static_cast<float>(tmpD[0]);
                  tmp[1] = static_cast<float>(tmpD[1]);
                  tmp[2] = static_cast<float>(tmpD[2]);
                  }
                }
              vect3d[0] = static_cast<double>(tmp[0]);
              vect3d[1] = static_cast<double>(tmp[1]);
              vect3d[2] = static_cast<double>(tmp[2]);
              }
            }

          if( StringContains(DiffusionDirectionality, "ISOTROPIC") )
            {
            continue;
            }
          else if( !B0FieldFound || b < zeroBValueTolerance )
            { // Deal with b0 images
            this->m_BValues.push_back(b);
            this->m_DiffusionVectors.push_back(vect3d);
            }
          else if( StringContains(DiffusionDirectionality, "DIRECTIONAL")
            || ( DiffusionDirectionality == "NONE" ) // Some new Philips data does not specify "DIRECTIONAL"
            || ( DiffusionDirectionality == "" ) )
            { // Deal with gradient direction images
            this->m_BValues.push_back(b);
            this->m_DiffusionVectors.push_back(vect3d);
            }
          else // Have no idea why we'd be here so error out
            {
            std::cout << "ERROR: DiffusionDirectionality was "
                      << DiffusionDirectionality << "  Don't know what to do with that..." << std::endl;
            throw;
            }

          std::cout << "B-value: " << b
                    << "; diffusion direction: "
                    << this->m_DoubleConvert(vect3d[0]) << ", "
                    << this->m_DoubleConvert(vect3d[1]) << ", "
                    << this->m_DoubleConvert(vect3d[2])
                    << std::endl;
          }
        }
      else
        {
        // multi-frame file, everything is inside
        std::map<std::vector<double>, double> gradientDirectionAndBValue;
        std::map<std::string, int> sliceLocations;
        std::vector<int> ignorePhilipsSliceMultiFrame;

        this->m_BValues.clear();
        this->m_DiffusionVectors.clear();

        itk::DCMTKSequence perFrameFunctionalGroup;
        itk::DCMTKSequence innerSeq;
        double             dwbValue;

        this->m_Headers[0]->GetElementSQ(0x5200, 0x9230, perFrameFunctionalGroup);
        this->m_NSlice = perFrameFunctionalGroup.card();

        // have to determine if volume slices are interleaved
        std::string origins[2];
        for( unsigned int i = 0; i < this->m_NSlice; ++i )
          {
          itk::DCMTKItem curItem;
          perFrameFunctionalGroup.GetElementItem(i, curItem);

          // index slice locations with string origin
          itk::DCMTKSequence originSeq;
          curItem.GetElementSQ(0x0020, 0x9113, originSeq);
          std::string originString;
          originSeq.GetElementDS(0x0020, 0x0032, originString);
          ++sliceLocations[originString];
          // save origin of first 2 slices to compare and see if the
          // volume is interleaved.
          if( i < 2 )
            {
            origins[i] = originString;
            }

          itk::DCMTKSequence mrDiffusionSeq;
          curItem.GetElementSQ(0x0018, 0x9117, mrDiffusionSeq);

          std::string dirValue;
          mrDiffusionSeq.GetElementCSorOB(0x0018, 0x9075, dirValue);

          if( StringContains(dirValue, "ISO") )
            {
            ignorePhilipsSliceMultiFrame.push_back(i);
            }
          else if( StringContains(dirValue, "NONE") )
            {
            std::vector<double> v(3);
            v[0] = 0; v[1] = 0; v[2] = 0;
            unsigned int nOld = gradientDirectionAndBValue.size();
            gradientDirectionAndBValue[v] = 0;
            unsigned int nNew = gradientDirectionAndBValue.size();

            if( nOld != nNew )
              {
              vnl_vector_fixed<double, 3> vect3d;
              vect3d.fill( 0 );
              this->m_DiffusionVectors.push_back( vect3d );
              this->m_BValues.push_back( 0 );
              }
            }
          else
            {
              {
              bool preferredExtrationSucceeded = EXIT_FAILURE;
              try
                {
                preferredExtrationSucceeded = mrDiffusionSeq.GetElementDSorOB(0x0018, 0x9087, dwbValue, false);
                }
              catch(...)
                {
                preferredExtrationSucceeded = EXIT_FAILURE;
                }
              //Try alternate method.
              if( preferredExtrationSucceeded == EXIT_FAILURE )
                {
                mrDiffusionSeq.GetElementFD(0x0018, 0x9087, dwbValue);
                }
              }
            itk::DCMTKSequence volSeq;
            mrDiffusionSeq.GetElementSQ(0x0018, 0x9076, volSeq);
            double dwgVal[3];
              {
              bool preferredExtrationSucceeded = EXIT_FAILURE;
              try
                {
                preferredExtrationSucceeded = volSeq.GetElementDSorOB<double>(0x0018, 0x9089, 3, dwgVal, false);
                }
              catch(...)
                {
                preferredExtrationSucceeded = EXIT_FAILURE;
                }
              //Try alternate method.
              if( preferredExtrationSucceeded == EXIT_FAILURE )
                {
                volSeq.GetElementFD(0x0018, 0x9089, 3, dwgVal);
                }
              }
            std::vector<double> v(3);
            v[0] = dwgVal[0];
            v[1] = dwgVal[1];
            v[2] = dwgVal[2];
            unsigned int nOld = gradientDirectionAndBValue.size();
            gradientDirectionAndBValue[v] = dwbValue;
            unsigned int nNew = gradientDirectionAndBValue.size();

            if( nOld != nNew )
              {
              vnl_vector_fixed<double, 3> vect3d;
              vect3d[0] = v[0]; vect3d[1] = v[1]; vect3d[2] = v[2];
              // vect3d.normalize();
              this->m_DiffusionVectors.push_back( vect3d );

              this->m_BValues.push_back( dwbValue);
              }
            }
          }
        // update values needed for (possible) de-interleave
        this->m_SlicesPerVolume = sliceLocations.size();


        std::cout << "LPS Matrix: " << std::endl << this->m_Volume->GetDirection() << std::endl;
        std::cout << "Volume Origin: " << std::endl << this->m_Volume->GetOrigin() << std::endl;
        std::cout << "Number of slices per volume: " << this->m_SlicesPerVolume << std::endl;
        std::cout << "Slice matrix size: " << this->m_Rows << " X " << this->m_Cols << std::endl;
        std::cout << "Image resolution: " << this->m_Volume->GetSpacing() << std::endl;

        this->m_MeasurementFrame = this->m_Volume->GetDirection();

        this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
        for( unsigned int k2 = 0; k2 < this->m_BValues.size(); ++k2 )
          {
          std::cout << k2 << ": direction: "
                    << this->m_DoubleConvert(this->m_DiffusionVectors[k2][0]) << ", "
                    << this->m_DoubleConvert(this->m_DiffusionVectors[k2][1]) << ", "
                    << this->m_DoubleConvert(this->m_DiffusionVectors[k2][2])
                    << ", b-value: " << this->m_BValues[k2] << std::endl;
          }
        // de-interleave slices if the origins of the first 2 slices
        // are the same.
        if( origins[0] == origins[1] )
          {
          // interleaved image
          DeInterleaveVolume();
          }
        }
      // deal with trailing isotropic images
      unsigned long trailingVolumes = this->m_NVolume - this->m_DiffusionVectors.size();
      if(trailingVolumes > 0)
        {
        std::cout << "# of Volumes " << this->m_NVolume << " # of Diffusion Vectors "
                  << this->m_DiffusionVectors.size() << " Removing "
                  << trailingVolumes << " Isotropic volumes." << std::endl;
        typedef itk::ExtractImageFilter<VolumeType,VolumeType> ExtractImageFilterType;
        ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();

        VolumeType::RegionType desiredRegion = this->m_Volume->GetLargestPossibleRegion();
        VolumeType::SizeType desiredSize = desiredRegion.GetSize();
        desiredSize[2] -= (trailingVolumes * this->m_SlicesPerVolume);
        desiredRegion.SetSize(desiredSize);
        extractImageFilter->SetExtractionRegion(desiredRegion);
        extractImageFilter->SetInput(this->m_Volume);
        extractImageFilter->Update();
        this->m_Volume = extractImageFilter->GetOutput();
        this->m_NVolume -= trailingVolumes;
        }
    }
protected:
  /** # of trailing images to ignore */
  unsigned int        m_NTrailingImagesToIgnore;

  virtual void AddFlagsToDictionary() ITK_OVERRIDE
    {
      // relevant Philips private tags
      DcmDictEntry *PhilipsDictBValue  = new DcmDictEntry(0x2001, 0x1003, DcmVR(EVR_FL),
                                                          "B Value of diffusion weighting", 1, 1, ITK_NULLPTR, true,
                                                          "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictBValue);
      DcmDictEntry *PhilipsDictDiffusionDirection   = new DcmDictEntry(0x2001, 0x1004, DcmVR(EVR_CS),
                                                                       "Diffusion Gradient Direction", 1, 1, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirection);

      DcmDictEntry *PhilipsDictDiffusionDirectionRL = new DcmDictEntry(0x2005, 0x10b0, DcmVR(EVR_FL),
                                                                       "Diffusion Direction R/L", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionRL);
      DcmDictEntry *PhilipsDictDiffusionDirectionAP = new DcmDictEntry(0x2005, 0x10b1, DcmVR(EVR_FL),
                                                                       "Diffusion Direction A/P", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionAP);
      DcmDictEntry *PhilipsDictDiffusionDirectionFH = new DcmDictEntry(0x2005, 0x10b2, DcmVR(EVR_FL),
                                                                       "Diffusion Direction F/H", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionFH);

      // New data new uses new tags!
      DcmDictEntry *PhilipsDictDiffusionDirectionRLnew = new DcmDictEntry(0x2005, 0x12b0, DcmVR(EVR_FL),
                                                                       "Diffusion Direction R/L", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionRLnew);
      DcmDictEntry *PhilipsDictDiffusionDirectionAPnew = new DcmDictEntry(0x2005, 0x12b1, DcmVR(EVR_FL),
                                                                       "Diffusion Direction A/P", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionAPnew);
      DcmDictEntry *PhilipsDictDiffusionDirectionFHnew = new DcmDictEntry(0x2005, 0x12b2, DcmVR(EVR_FL),
                                                                       "Diffusion Direction F/H", 4, 4, ITK_NULLPTR, true,
                                                                       "dicomtonrrd");
      itk::DCMTKFileReader::AddDictEntry(PhilipsDictDiffusionDirectionFHnew);

      // relevant Philips private tags
    }

};

#endif // __PhilipsDWIConverter_h
