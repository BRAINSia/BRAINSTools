#ifndef __SiemensDWIConverter_h
#define __SiemensDWIConverter_h
#include "DWIConverter.h"
#include "StringContains.h"

/** specific converter for Siemens scanners*/
class SiemensDWIConverter : public DWIConverter
{
public:
  SiemensDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections,
                      double smallGradientThreshold) : DWIConverter(allHeaders,inputFileNames,
                                                                    useBMatrixGradientDirections),
                                                       m_SmallGradientThreshold(smallGradientThreshold),
                                                       m_MMosaic(0),
                                                       m_NMosaic(0),
                                                       m_Stride(0)
    {
    }
  virtual ~SiemensDWIConverter() {}

  /** Siemens datasets are either in the
   *  normal sequential volume arrangement or
   *  in mosaic datasets, where each slice contains
   *  a collection of 2D slices arranged in a single
   *  mosaic slice.
   */
  virtual void LoadDicomDirectory()
    {
      this->DWIConverter::LoadDicomDirectory();
      std::string ImageType;
      this->m_Headers[0]->GetElementCS(0x0008, 0x0008, ImageType);
      if(StringContains(ImageType,"MOSAIC"))
        {
        this->m_NVolume = this->m_NSlice;
        this->m_Stride = 1; // Stride used in extracting the bval/gvec.
        std::cout << "Siemens SliceMosaic......" << std::endl;

        this->m_SliceOrderIS = false;

        // for siemens mosaic image, figure out mosaic slice order from 0029|1010
        // copy information stored in 0029,1010 into a string for parsing
        std::string tag;
        this->m_Headers[0]->GetElementOB(0x0029, 0x1010, tag);
        // parse SliceNormalVector from 0029,1010 tag
        std::vector<double> valueArray(0);
        int                 nItems = ExtractSiemensDiffusionInformation(tag, "SliceNormalVector", valueArray);
        if( nItems != 3 )  // did not find enough information
          {
          std::cout << "Warning: Cannot find complete information on SliceNormalVector in 0029|1010" << std::endl;
          std::cout << "         Slice order may be wrong." << std::endl;
          }
        else if( valueArray[2] > 0 )
          {
          m_SliceOrderIS = true;
          }

        // parse NumberOfImagesInMosaic from 0029,1010 tag
        valueArray.resize(0);
        nItems = ExtractSiemensDiffusionInformation(tag, "NumberOfImagesInMosaic", valueArray);
        if( nItems == 0 )  // did not find enough information
          {
          std::cout << "Warning: Cannot find complete information on NumberOfImagesInMosaic in 0029|1010" << std::endl;
          std::cout << "         Resulting image may contain empty slices." << std::endl;
          }
        else
          {
          this->m_SlicesPerVolume = static_cast<int>(valueArray[0]);
          this->m_MMosaic = static_cast<int>(ceil(sqrt(valueArray[0]) ) );
          this->m_NMosaic = this->m_MMosaic;
          }
        std::cout << "Mosaic in " << this->m_MMosaic << " X " << this->m_NMosaic
                  << " blocks (total number of blocks = " << valueArray[0] << ")." << std::endl;
        this->DeMosaic();
        }
      else
        {
        // expect normal 'array of volumes' organization
        // the superclass' LoadDicomDirectory will detect interleaved
        // slices organization and de-interleave the gradient volume.
        this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
        this->m_Stride = this->m_SlicesPerVolume; // Stride used in extracting the bval/gvec.
        this->DetermineSliceOrderIS();
        }
      this->SetDirectionsFromSliceOrder();
    }
  /** find the bvalues and gradient vectors */
  void ExtractDWIData()
    {
      for( unsigned int k = 0; k < this->m_NSlice; k += this->m_Stride )
        {
        // in Siemens, this entry is a 'CSA Header' which is blob
        // of mixed text & binary data.  Pretty annoying but there you
        // have it.
        std::string diffusionInfoString;;
        this->m_Headers[k]->GetElementOB( 0x0029, 0x1010, diffusionInfoString );

        // parse B_value from 0029,1010 tag
        std::vector<double> valueArray(0);

        int nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);
        if( nItems != 1 )
          {
          vnl_vector_fixed<double, 3> vect3d;
          // B_Value is missing -- the punt position is to count this
          // volume as having a B_value & Gradient Direction of zero
          std::cout << "Warning: Cannot find complete information on B_value in 0029|1010" << std::endl;
          this->m_BValues.push_back( 0.0 );
          vect3d.fill( 0.0 );
          this->m_DiffusionVectors.push_back(vect3d);
          continue;
          }

        // we got a 'valid' B-value
        // If we're trusting the gradient directions in the header,
        // then all we need to do here is save the bValue.
        if( !this->m_UseBMatrixGradientDirections )
          {
          valueArray.resize(0);
          ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);

          this->m_BValues.push_back( valueArray[0] );
          }
        else
          {
          // JTM - Patch from UNC: fill the nhdr header with the gradient directions and
          // bvalues computed out of the BMatrix
          valueArray.resize(0);
          nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_matrix", valueArray);
          vnl_matrix_fixed<double, 3, 3> bMatrix;

          if( nItems == 6 )
            {
            std::cout << "=============================================" << std::endl;
            std::cout << "BMatrix calculations..." << std::endl;
            // UNC comments: We get the value of the b-value tag in the header.
            // We won't use it as is, but just to locate the B0 images.
            // This check must be added, otherwise the bmatrix of the B0 is not
            // read properly (it's not an actual field in the DICOM header of the B0).
            std::vector<double> bval_tmp(0);
            bool                b0_image = false;

            // UNC comments: Get the bvalue
            nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", bval_tmp);
            if( bval_tmp[0] == 0 )
              {
              b0_image = true;
              }

            // UNC comments: The principal eigenvector of the bmatrix is to be extracted as
            // it's the gradient direction and trace of the matrix is the b-value

            // UNC comments: Fill out the 3x3 bmatrix with the 6 components read from the
            // DICOM header.
            bMatrix[0][0] = valueArray[0];
            bMatrix[0][1] = valueArray[1];
            bMatrix[0][2] = valueArray[2];
            bMatrix[1][1] = valueArray[3];
            bMatrix[1][2] = valueArray[4];
            bMatrix[2][2] = valueArray[5];
            bMatrix[1][0] = bMatrix[0][1];
            bMatrix[2][0] = bMatrix[0][2];
            bMatrix[2][1] = bMatrix[1][2];

            // UNC comments: Computing the decomposition
            vnl_svd<double> svd(bMatrix);

            // UNC comments: Extracting the principal eigenvector i.e. the gradient direction
            vnl_vector_fixed<double, 3> vect3d;
            vect3d[0] = svd.U(0, 0);
            vect3d[1] = svd.U(1, 0);
            vect3d[2] = svd.U(2, 0);

            std::cout << "BMatrix: " << std::endl;
            std::cout << bMatrix[0][0] << std::endl;
            std::cout << bMatrix[0][1] << "\t" << bMatrix[1][1] << std::endl;
            std::cout << bMatrix[0][2] << "\t" << bMatrix[1][2] << "\t" << bMatrix[2][2] << std::endl;

            // UNC comments: The b-value si the trace of the bmatrix
            const double bvalue = bMatrix[0][0] + bMatrix[1][1] + bMatrix[2][2];
            std::cout << bvalue << std::endl;
            // UNC comments: Even if the bmatrix is null, the svd decomposition set the 1st eigenvector
            // to (1,0,0). So we force the gradient direction to 0 if the bvalue is null
            if( (b0_image == true) || (bvalue == 0) )
              {
              std::cout << "B0 image detected: gradient direction and bvalue forced to 0" << std::endl;
              vect3d[0] = 0;
              vect3d[1] = 0;
              vect3d[2] = 0;
              std::cout << "Gradient coordinates: " << this->m_DoubleConvert(vect3d[0])
                        << " " << this->m_DoubleConvert(vect3d[1])
                        << " " << this->m_DoubleConvert(vect3d[2]) << std::endl;
              this->m_BValues.push_back(0);
              }
            else
              {
              std::cout << "Gradient coordinates: " << this->m_DoubleConvert(vect3d[0])
                        << " " << this->m_DoubleConvert(vect3d[1])
                        << " " << this->m_DoubleConvert(vect3d[2]) << std::endl;
              this->m_BValues.push_back(bvalue);
              }
            this->m_DiffusionVectors.push_back(vect3d);
            }
          else
            {
            // silently returning zero gradient vectors is a problem,
            // but it is also necessary for some fiels.
            valueArray.resize(0);
            ExtractSiemensDiffusionInformation(diffusionInfoString, "B_value", valueArray);
            vnl_vector_fixed<double, 3> vect3d;
            this->m_BValues.push_back( valueArray[0] );
            vect3d[0] = 0;
            vect3d[1] = 0;
            vect3d[2] = 0;
            this->m_DiffusionVectors.push_back(vect3d);
            }
          }
        }

      if( this->m_UseBMatrixGradientDirections == false )
        {
        for( unsigned int k = 0; k < this->m_NSlice; k += this->m_Stride )
          {
          std::cout << "=======================================" << std::endl << std::endl;
          std::string diffusionInfoString;
          this->m_Headers[k]->GetElementOB(0x0029, 0x1010, diffusionInfoString );

          std::vector<double>         valueArray;
          vnl_vector_fixed<double, 3> vect3d;

          // parse DiffusionGradientDirection from 0029,1010 tag
          valueArray.resize(0);
          int nItems =
            ExtractSiemensDiffusionInformation(diffusionInfoString, "DiffusionGradientDirection", valueArray);
          if( nItems != 3 )  // did not find enough information
            {
            std::cout << "Warning: Cannot find complete information on DiffusionGradientDirection in 0029|1010"
                      << std::endl;
            vect3d.fill( 0 );
            this->m_DiffusionVectors.push_back(vect3d);
            }
          else
            {
            std::cout << "Number of Directions : " << nItems << std::endl;
            std::cout << "   Directions 0: " << valueArray[0] << std::endl;
            std::cout << "   Directions 1: " << valueArray[1] << std::endl;
            std::cout << "   Directions 2: " << valueArray[2] << std::endl;
            double DiffusionVector_magnitude;
            vect3d[0] = valueArray[0];
            vect3d[1] = valueArray[1];
            vect3d[2] = valueArray[2];

            DiffusionVector_magnitude = sqrt((vect3d[0] * vect3d[0]) +
                                             (vect3d[1] * vect3d[1]) +
                                             (vect3d[2] * vect3d[2]) );

            std::cout << "DiffusionVector_magnitude " << DiffusionVector_magnitude << std::endl;
            if( DiffusionVector_magnitude <= this->m_SmallGradientThreshold )
              {
              std::cout << "ERROR: Gradient vector with unreasonably small magnitude exists." << std::endl;
              std::cout << "Gradient #" << k << " with magnitude " << DiffusionVector_magnitude << std::endl;
              std::cout << "Please set useBMatrixGradientDirections to calculate gradient directions "
                        << "from the scanner B Matrix to alleviate this problem." << std::endl;
              throw;
              }

            // vect3d.normalize();
            this->m_DiffusionVectors.push_back(vect3d);
            int p = this->m_BValues.size();
            std::cout << "Image#: " << k
                      << " BV: " << this->m_BValues[p - 1] << " GD: "
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][0]) << ","
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][1]) << ","
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][2])
                      << std::endl;
            }
          }
        }
      //
      // test gradients. It is OK for one or more guide images to have
      // zero gradients, but all gradients == 0 is an error. It means
      // that the gradient data is missing.
      DiffusionVecVectorType::iterator nonZ =
        std::find_if(this->m_DiffusionVectors.begin(),
                     this->m_DiffusionVectors.end(),
                     SiemensDWIConverter::IsZeroMag);
      if(nonZ == this->m_DiffusionVectors.end())
        {
        itkGenericExceptionMacro(<< "Dataset has no diffusion vectors");
        }
    }
private:
  static bool IsZeroMag(DiffusionVectorType vec)
    {
      return vec.magnitude() != 0.0;
    }
protected:
  /** turn a mosaic image back into a sequential volume image */
  void DeMosaic()
    {
      // de-mosaic
      this->m_Rows /= this->m_MMosaic;
      this->m_Cols /= this->m_NMosaic;

      // center the volume since the image position patient given in the
      // dicom header was useless
      this->m_Origin[0] = -(this->m_Rows * (this->m_NRRDSpaceDirection[0][0])
                         + this->m_Cols * (this->m_NRRDSpaceDirection[0][1])
                         + this->m_SlicesPerVolume * (this->m_NRRDSpaceDirection[0][2]) ) / 2.0;
      this->m_Origin[1] = -(this->m_Rows * (this->m_NRRDSpaceDirection[1][0])
                         + this->m_Cols * (this->m_NRRDSpaceDirection[1][1])
                         + this->m_SlicesPerVolume * (this->m_NRRDSpaceDirection[1][2]) ) / 2.0;
      this->m_Origin[2] = -(this->m_Rows * (this->m_NRRDSpaceDirection[2][0])
                         + this->m_Cols * (this->m_NRRDSpaceDirection[2][1])
                         + this->m_SlicesPerVolume * (this->m_NRRDSpaceDirection[2][2]) ) / 2.0;

      VolumeType::Pointer img = this->m_Volume;

      VolumeType::RegionType region = img->GetLargestPossibleRegion();
      VolumeType::SizeType   size = region.GetSize();

      VolumeType::SizeType dmSize = size;
      unsigned int         original_slice_number = dmSize[2] * m_SlicesPerVolume;
      dmSize[0] /= this->m_MMosaic;
      dmSize[1] /= this->m_NMosaic;
      dmSize[2] = this->m_NVolume * this->m_SlicesPerVolume;

      region.SetSize( dmSize );
      this->m_Volume = VolumeType::New();
      this->m_Volume->CopyInformation( img );
      this->m_Volume->SetRegions( region );
      this->m_Volume->Allocate();

      VolumeType::RegionType dmRegion = this->m_Volume->GetLargestPossibleRegion();
      dmRegion.SetSize(2, 1);
      region.SetSize(0, dmSize[0]);
      region.SetSize(1, dmSize[1]);
      region.SetSize(2, 1);

      for( unsigned int k = 0; k < original_slice_number; ++k )
        {
        unsigned int new_k = k /* - bad_slice_counter */;

        dmRegion.SetIndex(2, new_k);
        itk::ImageRegionIteratorWithIndex<VolumeType> dmIt( this->m_Volume, dmRegion );

        // figure out the mosaic region for this slice
        int sliceIndex = k;

        // int nBlockPerSlice = this->m_Mosaic*this->m_NMosaic;
        int slcMosaic = sliceIndex / (m_SlicesPerVolume);
        sliceIndex -= slcMosaic * m_SlicesPerVolume;
        int colMosaic = sliceIndex / this->m_MMosaic;
        int rawMosaic = sliceIndex - this->m_MMosaic * colMosaic;
        region.SetIndex( 0, rawMosaic * dmSize[0] );
        region.SetIndex( 1, colMosaic * dmSize[1] );
        region.SetIndex( 2, slcMosaic );

        itk::ImageRegionConstIteratorWithIndex<VolumeType> imIt( img, region );
        for( dmIt.GoToBegin(), imIt.GoToBegin(); !dmIt.IsAtEnd(); ++dmIt, ++imIt )
          {
          dmIt.Set( imIt.Get() );
          }
        }
    }

  unsigned int ConvertFromCharPtr(const char *s)
    {
      unsigned int rval = 0;

      // assume little-endian
      for( unsigned i = 0; i < sizeof(unsigned int); ++i )
        {
        rval += ( (unsigned int)s[i]) << (i * 8);
        }
      return rval;
    }

  /** pull data out of Siemens scans.
   *
   *  Siemens sticks most of the DTI information into a single
   *  OB-format entry.  This is actually rigidly structured, but
   *  this function depends on the needed data living at fixed offset
   *  from the beginning of the name of each tag, and ignores the
   *  internal structure documented in the Siemens Dicom Compliance
   *  document.
   */
  unsigned int
  ExtractSiemensDiffusionInformation(const std::string & tagString,
                                     const std::string & nameString,
                                     std::vector<double>& valueArray)
    {
      std::cerr << "TagString: " << tagString << std::endl;
      ::size_t atPosition = tagString.find( nameString );

      if( atPosition == std::string::npos )
        {
        return 0;
        }
      while( true )  // skip nameString inside a quotation
        {
        std::string nextChar = tagString.substr( atPosition + nameString.size(), 1 );

        if( nextChar.c_str()[0] == 0 )
          {
          break;
          }
        else
          {
          atPosition = tagString.find( nameString, atPosition + 2 );
          }
        }

      if( atPosition == std::string::npos )
        {
        return 0;
        }
      std::string  infoAsString = tagString.substr( atPosition, tagString.size() - atPosition + 1 );
      const char * infoAsCharPtr = infoAsString.c_str();

      unsigned int vm = ConvertFromCharPtr(infoAsCharPtr + 64);
      {
      std::string vr = infoAsString.substr( 68, 2 );

      // std::cout << "\tName String: " << nameString << std::endl;
      // std::cout << "\tVR: " << vr << std::endl;
      // std::cout << "\tVM: " << vm << std::endl;
      // std::cout << "Local String: " << infoAsString.substr(0,80) << std::endl;

      /* This hack is required for some Siemens VB15 Data */
      if( ( nameString == "DiffusionGradientDirection" ) && (vr != "FD") )
        {
        bool loop = true;
        while( loop )
          {
          atPosition = tagString.find( nameString, atPosition + 26 );
          if( atPosition == std::string::npos )
            {
            // std::cout << "\tFailed to find DiffusionGradientDirection Tag - returning" << vm << std::endl;
            return 0;
            }
          infoAsString = tagString.substr( atPosition, tagString.size() - atPosition + 1 );
          infoAsCharPtr = infoAsString.c_str();
          // std::cout << "\tOffset to new position" << std::endl;
          // std::cout << "\tNew Local String: " << infoAsString.substr(0,80) << std::endl;
          vm = ConvertFromCharPtr(infoAsCharPtr + 64);
          vr = infoAsString.substr( 68, 2 );
          if( vr == "FD" )
            {
            loop = false;
            }
          // std::cout << "\tVR: " << vr << std::endl;
          // std::cout << "\tVM: " << vm << std::endl;
          }
        }
      else
        {
        // std::cout << "\tUsing initial position" << std::endl;
        }
      // std::cout << "\tArray Length: " << vm << std::endl;
      }

      unsigned int offset = 84;
      for( unsigned int k = 0; k < vm; ++k )
        {
        const int itemLength = ConvertFromCharPtr(infoAsCharPtr + offset + 4);
        const int strideSize = static_cast<int>(ceil(static_cast<double>(itemLength) / 4) * 4);
        if( infoAsString.length() < offset + 16 + itemLength )
          {
          // data not available or incomplete
          return 0;
          }
        const std::string valueString = infoAsString.substr( offset + 16, itemLength );
        valueArray.push_back( atof(valueString.c_str() ) );
        offset += 16 + strideSize;
        }
      return vm;
    }
  virtual void AddFlagsToDictionary()
    {
      // relevant Siemens private tags
      DcmDictEntry *SiemensMosiacParameters = new DcmDictEntry(0x0051, 0x100b, DcmVR(EVR_IS),
                                                               "Mosiac Matrix Size", 1, 1, 0, true,
                                                               "dicomtonrrd");
      DcmDictEntry *SiemensDictNMosiac = new DcmDictEntry(0x0019, 0x100a, DcmVR(EVR_US),
                                                          "Number of Images In Mosaic", 1, 1, 0, true,
                                                          "dicomtonrrd");
      DcmDictEntry *SiemensDictBValue = new DcmDictEntry(0x0019, 0x100c, DcmVR(EVR_IS),
                                                         "B Value of diffusion weighting", 1, 1, 0, true,
                                                         "dicomtonrrd");
      DcmDictEntry *SiemensDictDiffusionDirection = new DcmDictEntry(0x0019, 0x100e, DcmVR(EVR_FD),
                                                                     "Diffusion Gradient Direction", 3, 3, 0, true,
                                                                     "dicomtonrrd");
      DcmDictEntry *SiemensDictDiffusionMatrix = new DcmDictEntry(0x0019, 0x1027, DcmVR(EVR_FD),
                                                                  "Diffusion Matrix", 6, 6, 0, true,
                                                                  "dicomtonrrd");
      DcmDictEntry *SiemensDictShadowInfo = new DcmDictEntry(0x0029, 0x1010, DcmVR(EVR_OB),
                                                             "Siemens DWI Info", 1, 1, 0, true,
                                                             "dicomtonrrd");

      // relevant Siemens private tags
      itk::DCMTKFileReader::AddDictEntry(SiemensMosiacParameters);
      itk::DCMTKFileReader::AddDictEntry(SiemensDictNMosiac);
      itk::DCMTKFileReader::AddDictEntry(SiemensDictBValue);
      itk::DCMTKFileReader::AddDictEntry(SiemensDictDiffusionDirection);
      itk::DCMTKFileReader::AddDictEntry(SiemensDictDiffusionMatrix);
      itk::DCMTKFileReader::AddDictEntry(SiemensDictShadowInfo);
    }
private:
  double      m_SmallGradientThreshold;
  unsigned int m_MMosaic;
  unsigned int m_NMosaic;
  unsigned int m_Stride;
};

#endif // __SiemensDWIConverter_h
