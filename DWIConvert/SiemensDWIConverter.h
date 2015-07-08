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

  template <typename T>
  T CSAExtractFromString(const char *ptr)
    {
      T rval = *(reinterpret_cast<const T *>(ptr));
      itk::ByteSwapper<T>::SwapFromSystemToLittleEndian(&rval);
      return rval;
    }

  class CSAItem : public std::vector<std::string >
  {
  public:
    typedef std::vector<std::string> SuperClass;

    itk::uint32_t vm;
    std::string vr;

    CSAItem(unsigned int length) : SuperClass(length),
                                   vm(0)
      {
      }
    CSAItem() : SuperClass()
      {
      }
    CSAItem( const CSAItem & other ) : SuperClass(other.size())
      {
        *this = other;
      }
    CSAItem & operator=(const CSAItem &other)
      {
        this->resize(0);
        for(CSAItem::const_iterator it = other.begin();
            it != other.end();
            ++it)
          {
          this->push_back(*it);
          }
        this->vm = other.vm;
        this->vr = other.vr;
        return *this;
      }

    template <typename T>
    std::vector<T> AsVector() const
      {
        std::vector<T> rval;
        for(unsigned i = 0; i < this->size(); ++i)
          {
          if(! (*this)[i].empty())
            {
            T val = 0;
            std::stringstream convert((*this)[i]);
            convert >> val;
            rval.push_back(val);
            }
          }
        return rval;
      }
    void DebugPrint() const
      {
        std::cerr << "  VM = " << this->vm << " VR = " << this->vr << std::endl
                  << "    ";
        bool firstTime(false);
        for(CSAItem::const_iterator it = this->begin();
            it != this->end(); ++it)
          {
          if(firstTime)
            {
            firstTime = false;
            }
          else
            {
            std::cerr << " ";
            }
          std::cerr << *it;
          }
        std::cerr << std::endl;
      }
  };

  class CSAHeader : public std::map<std::string,CSAItem>
  {
  public:
    void DebugPrint() const
      {
        for(CSAHeader::const_iterator it = this->begin();
            it != this->end(); ++it)
          {
          std::cerr << it->first << std::endl;
          it->second.DebugPrint();
          }
      }
  };

  void  DecodeCSAHeader(CSAHeader &header, const std::string &infoString)
    {
      //
      // the reference used to write this code is here:
      // http://nipy.sourceforge.net/nibabel/dicom/siemens_csa.html
      const char *info = infoString.c_str();

      const bool isCSA2 = info[0] == 'S' && info[1] == 'V'
        && info[2] == '1' && info[3] == '0';
      unsigned int offset;

      if(isCSA2)
        {
        offset = 8; // past SV10 + unused 4 bytes
        }
      else
        {
        offset = 0;
        }
      const itk::uint32_t numberOfTags =
        this->CSAExtractFromString<itk::uint32_t>(info+offset);
      offset += sizeof(itk::uint32_t); // skip numberOfTags;
      offset += sizeof(itk::uint32_t); // skip unused2

      for(unsigned i = 0; i < numberOfTags; ++i)
        {
        // tag name is 64 bytes null terminated.
        std::string tagName = info + offset;
        offset += 64;                        // skip tag name
        itk::uint32_t vm = this->CSAExtractFromString<itk::uint32_t>(info+offset);
        offset += sizeof(itk::uint32_t);

        CSAItem current(vm);
        current.vm = vm;

        // vr = 3 bytes of string + 1 for pad
        char vr[4];
        for(unsigned j = 0; j < 3; ++j)
          {
          vr[j] = info[offset + j];
          }
        vr[3] = '\0';
        current.vr = vr;
        offset += 4; // after VR
        offset += 4; // skip syngodt

        const itk::int32_t nItems =
          this->CSAExtractFromString<itk::int32_t>(info + offset);
        offset += 4;
        offset += 4; // skip xx

        for(int j = 0; j < nItems; ++j)
          {
          // 4 items in XX, first being item length
          const itk::int32_t  itemLength =
            this->CSAExtractFromString<itk::int32_t>(info + offset);
          offset += 16;
          std::string valueString;
          valueString = info + offset;
          offset += itemLength;
          while((offset % 4) != 0)
            {
            ++offset;
            }
          if(j < static_cast<int>(vm))
            {
            current[j] = valueString;
            }
          }
        header[tagName] = current;
        }
    }

  /** Siemens datasets are either in the
   *  normal sequential volume arrangement or
   *  in mosaic datasets, where each slice contains
   *  a collection of 2D slices arranged in a single
   *  mosaic slice.
   */
  virtual void LoadDicomDirectory() ITK_OVERRIDE
    {
      this->DWIConverter::LoadDicomDirectory();
      std::string ImageType;
      this->m_MeasurementFrame.SetIdentity();
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
        this->DetermineSliceOrderIS();
        this->SetDirectionsFromSliceOrder();
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
        this->SetDirectionsFromSliceOrder();
        }
    }
  /** find the bvalues and gradient vectors */
  void ExtractDWIData() ITK_OVERRIDE
    {
      for( unsigned int k = 0; k < this->m_NSlice; k += this->m_Stride )
        {
        // in Siemens, this entry is a 'CSA Header' which is blob
        // of mixed text & binary data.  Pretty annoying but there you
        // have it.
        std::string diffusionInfoString;;
        this->m_Headers[k]->GetElementOB( 0x0029, 0x1010, diffusionInfoString );

        CSAHeader csaHeader;
        this->DecodeCSAHeader(csaHeader,diffusionInfoString);
        // csaHeader.DebugPrint();

        // parse B_value from 0029,1010 tag
        std::vector<double> valueArray(0);
        CSAHeader::const_iterator csaIt;
        if( !this->m_UseBMatrixGradientDirections )
          {
          int nItems(0);
          if((csaIt = csaHeader.find("B_value")) == csaHeader.end())
            {
            nItems = 0;
            }
          else
            {
            // we got a 'valid' B-value
            // If we're trusting the gradient directions in the header,
            // then all we need to do here is save the bValue.
            valueArray = csaIt->second.AsVector<double>();
            nItems = valueArray.size();
            }
          vnl_vector_fixed<double, 3> vect3d;
          vect3d.fill( 0.0 );

          double currentBValue = 0.0;
          if( nItems == 1 )
            {
            currentBValue = valueArray[0] ;
            }
          else
            {
            // B_Value is missing -- the punt position is to count this
            // volume as having a B_value & Gradient Direction of zero
            std::cout << "Warning: Cannot find complete information on B_value in 0029|1010" << std::endl;
            }

          std::cout << "---------- BVALUE: " << currentBValue << std::endl;
          if ( currentBValue < 1.0 )
          {
            this->m_BValues.push_back( 0.0 );
            this->m_DiffusionVectors.push_back(vect3d);
            continue;
          }
          this->m_BValues.push_back( currentBValue );
          if((csaIt = csaHeader.find("DiffusionGradientDirection")) == csaHeader.end())
            {
            nItems = 0;
            }
          else
            {
            valueArray = csaIt->second.AsVector<double>();
            nItems = valueArray.size();
            }

          // parse DiffusionGradientDirection from 0029,1010 tag
          valueArray.resize(0);
          nItems = ExtractSiemensDiffusionInformation(diffusionInfoString, "DiffusionGradientDirection", valueArray);
          if( nItems != 3 )  // did not find enough information
            {
            std::cout << "Warning: Cannot find complete information on DiffusionGradientDirection in 0029|1010"
                      << std::endl;
            vect3d.fill( 0 );
            this->m_DiffusionVectors.push_back(vect3d);
            }
          else
            {
            double DiffusionVector_magnitude = sqrt((valueArray[0] * valueArray[0]) +
              (valueArray[1] * valueArray[1]) +
              (valueArray[2] * valueArray[2]) );
            if ( DiffusionVector_magnitude > 1.1 )
              {
              //Gradient vectors are supposed to be unit vectors!
              // If coded as [ 1.0001 1.0001 1.0001 ]  then it is really a B0 image.
              // This is ugly hack but works around a persistent dicom coding problem
              // on some scanners
              std::fill(valueArray.begin(), valueArray.end(), 0.0);
              DiffusionVector_magnitude=0.0;
              }
            else if( DiffusionVector_magnitude <= this->m_SmallGradientThreshold )
              {
              std::cout << "ERROR: Gradient vector with unreasonably small magnitude exists." << std::endl;
              std::cout << "Gradient #" << k << " with magnitude " << DiffusionVector_magnitude << std::endl;
              std::cout << "Please set useBMatrixGradientDirections to calculate gradient directions "
                        << "from the scanner B Matrix to alleviate this problem." << std::endl;
              throw;
              }

            std::cout << "Number of Directions : " << nItems << std::endl;
            std::cout << "   Directions 0: " << valueArray[0] << std::endl;
            std::cout << "   Directions 1: " << valueArray[1] << std::endl;
            std::cout << "   Directions 2: " << valueArray[2] << std::endl;
            std::cout << "DiffusionVector_magnitude " << DiffusionVector_magnitude << std::endl;

            vect3d[0] = valueArray[0];
            vect3d[1] = valueArray[1];
            vect3d[2] = valueArray[2];
            // vect3d.normalize();
            this->m_DiffusionVectors.push_back(vect3d);
            const int p = this->m_BValues.size();
            std::cout << "Image#: " << k
                      << " BV: " << this->m_BValues[p - 1] << " GD: "
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][0]) << ","
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][1]) << ","
                      << this->m_DoubleConvert(this->m_DiffusionVectors[k / this->m_Stride][2])
                      << std::endl;
            }
          }
        else
          {
          int nItems;
          // UNC comments: We get the value of the b-value tag in the header.
          // We won't use it as is, but just to locate the B0 images.
          // This check must be added, otherwise the bmatrix of the B0 is not
          // read properly (it's not an actual field in the DICOM header of the B0).
          std::vector<double> bval_tmp(0);
          bool                b0_image = false;
          // UNC comments: Get the bvalue
          if((csaIt = csaHeader.find("B_value")) == csaHeader.end())
            {
            nItems = 0;
            }
          else
            {
            // we got a 'valid' B-value
            // If we're trusting the gradient directions in the header,
            // then all we need to do here is save the bValue.
            bval_tmp  = csaIt->second.AsVector<double>();
            }

          if( bval_tmp[0] == 0 )
            {
            b0_image = true;
            }

          // JTM - Patch from UNC: fill the nhdr header with the gradient directions and
          // bvalues computed out of the BMatrix
          if((csaIt = csaHeader.find("B_matrix")) == csaHeader.end())
            {
            nItems = 0;
            }
          else
            {
            valueArray = csaIt->second.AsVector<double>();
            nItems = valueArray.size();
            }
          vnl_matrix_fixed<double, 3, 3> bMatrix;

          if( nItems == 6 && !b0_image)
            {
            std::cout << "=============================================" << std::endl;
            std::cout << "BMatrix calculations..." << std::endl;

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
            this->m_BValues.push_back( bval_tmp[0] );
            vnl_vector_fixed<double, 3> vect3d;
            vect3d[0] = 0;
            vect3d[1] = 0;
            vect3d[2] = 0;
            this->m_DiffusionVectors.push_back(vect3d);
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
        std::cerr << this->m_InputFileNames[0] << " has no non-zero diffusion vectors" << std::endl;
        // itkGenericExceptionMacro(<< "" );
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
      this->m_Cols /= this->m_MMosaic;
      this->m_Rows /= this->m_NMosaic;

      // center the volume since the image position patient given in the
      // dicom header was useless
      this->m_Origin[0] = -(this->m_Cols * (this->m_NRRDSpaceDirection[0][0])
                         + this->m_Rows * (this->m_NRRDSpaceDirection[0][1])
                         + this->m_SlicesPerVolume * (this->m_NRRDSpaceDirection[0][2]) ) / 2.0;
      this->m_Origin[1] = -(this->m_Cols * (this->m_NRRDSpaceDirection[1][0])
                         + this->m_Rows * (this->m_NRRDSpaceDirection[1][1])
                         + this->m_SlicesPerVolume * (this->m_NRRDSpaceDirection[1][2]) ) / 2.0;
      this->m_Origin[2] = -(this->m_Cols * (this->m_NRRDSpaceDirection[2][0])
                         + this->m_Rows * (this->m_NRRDSpaceDirection[2][1])
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
        const size_t infoAsStringLength=infoAsString.length();

        if( infoAsStringLength < ( offset + 16 + itemLength ) )
          {
          // data not available or incomplete
          return 0;
          }
        const std::string valueString = infoAsString.substr( offset + 16, itemLength );
        const double componentValue =  atof(valueString.c_str() );
        valueArray.push_back( componentValue );
        offset += 16 + strideSize;
        }
      return vm;
    }
  virtual void AddFlagsToDictionary() ITK_OVERRIDE
    {
      // relevant Siemens private tags
      DcmDictEntry *SiemensMosiacParameters = new DcmDictEntry(0x0051, 0x100b, DcmVR(EVR_IS),
                                                               "Mosiac Matrix Size", 1, 1, ITK_NULLPTR, true,
                                                               "dicomtonrrd");
      DcmDictEntry *SiemensDictNMosiac = new DcmDictEntry(0x0019, 0x100a, DcmVR(EVR_US),
                                                          "Number of Images In Mosaic", 1, 1, ITK_NULLPTR, true,
                                                          "dicomtonrrd");
      DcmDictEntry *SiemensDictBValue = new DcmDictEntry(0x0019, 0x100c, DcmVR(EVR_IS),
                                                         "B Value of diffusion weighting", 1, 1, ITK_NULLPTR, true,
                                                         "dicomtonrrd");
      DcmDictEntry *SiemensDictDiffusionDirection = new DcmDictEntry(0x0019, 0x100e, DcmVR(EVR_FD),
                                                                     "Diffusion Gradient Direction", 3, 3, ITK_NULLPTR, true,
                                                                     "dicomtonrrd");
      DcmDictEntry *SiemensDictDiffusionMatrix = new DcmDictEntry(0x0019, 0x1027, DcmVR(EVR_FD),
                                                                  "Diffusion Matrix", 6, 6, ITK_NULLPTR, true,
                                                                  "dicomtonrrd");
      DcmDictEntry *SiemensDictShadowInfo = new DcmDictEntry(0x0029, 0x1010, DcmVR(EVR_OB),
                                                             "Siemens DWI Info", 1, 1, ITK_NULLPTR, true,
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
