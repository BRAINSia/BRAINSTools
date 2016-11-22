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
#ifndef __DWIConverter_h
#define __DWIConverter_h
#include <vector>
#include <algorithm>
#include "itkMatrix.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkDCMTKFileReader.h"
#include "itkDCMTKImageIO.h"
#include "StringContains.h"
#include <algorithm>
#include "DWIMetaDataDictionaryValidator.h"

/** the DWIConverter is a base class for all scanner-specific
 *  converters.  It handles the tasks that are required for all
 *  scanners. In particular it loads the DICOM directory, and fills
 *  out various data fields needed by the DWIConvert program in order
 *  to write out NRRD and other files.
 */
class DWIConverter
{
public:
  typedef short                               PixelValueType;
  typedef itk::Image<PixelValueType, 3>       VolumeType;
  typedef VolumeType::SpacingType             SpacingType;
  typedef itk::ImageSeriesReader<VolumeType>  ReaderType;
  typedef ReaderType::FileNamesContainer      FileNamesContainer;
  typedef itk::ImageFileReader<VolumeType>    SingleFileReaderType;
  typedef itk::DCMTKSeriesFileNames           InputNamesGeneratorType;
  typedef std::vector<itk::DCMTKFileReader *> DCMTKFileVector;
  typedef itk::Matrix<double, 3, 3>           RotationMatrixType;
  typedef itk::Vector<double, 3>              PointType;
  DWIConverter(const DCMTKFileVector &allHeaders,
               const FileNamesContainer &inputFileNames,
               const bool useBMatrixGradientDirections) : m_Headers(allHeaders),
                                                    m_InputFileNames(inputFileNames),
                                                    m_Rows(0),
                                                    m_Cols(0),
                                                    m_SlicesPerVolume(0),
                                                    m_XRes(0.0),
                                                    m_YRes(0.0),
                                                    m_SliceSpacing(0.0),
                                                    m_MultiSliceVolume(false),
                                                    m_SliceOrderIS(true),
                                                    m_NSlice(0),
                                                    m_NVolume(0),
                                                    m_UseBMatrixGradientDirections(useBMatrixGradientDirections),
                                                    m_IsInterleaved(false)
    {
      { // Set the m_Origin only 1 time, and then use it every time
      // origin
      double origin[3];
      m_Headers[0]->GetOrigin(origin);
      this->m_Origin[0] = origin[0];
      this->m_Origin[1] = origin[1];
      this->m_Origin[2] = origin[2];
      }

      this->m_NRRDSpaceDefinition = "left-posterior-superior";;
      this->m_MeasurementFrame.SetIdentity();
      this->m_LPSDirCos.SetIdentity();
      this->m_SpacingMatrix.SetIdentity();
    }

  virtual ~DWIConverter() {}

  virtual void LoadDicomDirectory()
    {
      //
      // add vendor-specific flags to dictionary
      this->AddFlagsToDictionary();
      //
      // load the volume, either single or multivolume.
      m_NSlice = this->m_InputFileNames.size();
      itk::DCMTKImageIO::Pointer dcmtkIO = itk::DCMTKImageIO::New();
      if( this->m_InputFileNames.size() > 1 )
        {
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetImageIO( dcmtkIO );
        reader->SetFileNames( this->m_InputFileNames );
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cerr << "Exception thrown while reading DICOM volume"
                    << std::endl;
          std::cerr << excp << std::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = false;
        }
      else
        {
        SingleFileReaderType::Pointer reader =
          SingleFileReaderType::New();
        reader->SetImageIO( dcmtkIO );
        reader->SetFileName( this->m_InputFileNames[0] );
        m_NSlice = this->m_InputFileNames.size();
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::cerr << "Exception thrown while reading the series" << std::endl;
          std::cerr << excp << std::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = true;
        }

      // figure out image dimensions
      m_Headers[0]->GetElementUS(0x0028, 0x0010, this->m_Rows);
      m_Headers[0]->GetElementUS(0x0028, 0x0011, this->m_Cols);

      // spacing
      double spacing[3];
      m_Headers[0]->GetSpacing(spacing);
      m_YRes = spacing[1];
      m_XRes = spacing[0];
      m_SliceSpacing = spacing[2];

      // a map of ints keyed by the slice location string
      // reported in the dicom file.  The number of slices per
      // volume is the same as the number of unique slice locations
      std::map<std::string, int> sliceLocations;
      //
      // check for interleave
      if( !this->m_MultiSliceVolume )
        {
        // Make a hash of the sliceLocations in order to get the correct
        // count.  This is more reliable since SliceLocation may not be available.
        std::vector<int>         sliceLocationIndicator;
        std::vector<std::string> sliceLocationStrings;

        sliceLocationIndicator.resize( this->m_NSlice );
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          std::string originString;
          this->m_Headers[k]->GetElementDS(0x0020, 0x0032, originString );
          sliceLocationStrings.push_back( originString );
          sliceLocations[originString]++;
          }

        // this seems like a crazy way to figure out if slices are
        // interleaved, but it works. Perhaps replace with comparing
        // the reported location between the first two slices?
        // Would be less clever-looking and devious, but would require
        // less computation.
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          std::map<std::string, int>::iterator it = sliceLocations.find( sliceLocationStrings[k] );
          sliceLocationIndicator[k] = distance( sliceLocations.begin(), it );
          }

        // sanity check on # of volumes versus # of dicom files
        if(this->m_Headers.size() % sliceLocations.size() != 0)
          {
          itkGenericExceptionMacro(<< "Missing DICOM Slice files: Number of slice files ("
                            << this->m_Headers.size() << ") not evenly divisible by"
                            << " the number of slice locations ");
          }

        this->m_SlicesPerVolume = sliceLocations.size();
        std::cout << "=================== this->m_SlicesPerVolume:" << this->m_SlicesPerVolume << std::endl;


        // if the this->m_SlicesPerVolume == 1, de-interleaving won't do
        // anything so there's no point in doing it.
        if( this->m_NSlice >= 2 && this->m_SlicesPerVolume > 1 )
          {
          if( sliceLocationIndicator[0] != sliceLocationIndicator[1] )
            {
            std::cout << "Dicom images are ordered in a volume interleaving way." << std::endl;
            }
          else
            {
            std::cout << "Dicom images are ordered in a slice interleaving way." << std::endl;
            this->m_IsInterleaved = true;
            // reorder slices into a volume interleaving manner
            DeInterleaveVolume();
            }
          }
        }

    // check ImageOrientationPatient and figure out slice direction in
    // L-P-I (right-handed) system.
    // In Dicom, the coordinate frame is L-P by default. Look at
    // http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301
    double dirCosArray[6];
    // 0020,0037 -- Image Orientation (Patient)
    this->m_Headers[0]->GetDirCosArray(dirCosArray);
    double *dirCosArrayP = dirCosArray;
    for( unsigned i = 0; i < 2; ++i )
      {
      for( unsigned j = 0; j < 3; ++j, ++dirCosArrayP )
        {
        this->m_LPSDirCos[j][i] = *dirCosArrayP;
        }
      }

    // Cross product, this gives I-axis direction
    this->m_LPSDirCos[0][2] = this->m_LPSDirCos[1][0] * this->m_LPSDirCos[2][1] -
      this->m_LPSDirCos[2][0] * this->m_LPSDirCos[1][1];
    this->m_LPSDirCos[1][2] = this->m_LPSDirCos[2][0] * this->m_LPSDirCos[0][1] -
      this->m_LPSDirCos[0][0] * this->m_LPSDirCos[2][1];
    this->m_LPSDirCos[2][2] = this->m_LPSDirCos[0][0] * this->m_LPSDirCos[1][1] -
      this->m_LPSDirCos[1][0] * this->m_LPSDirCos[0][1];

    std::cout << "ImageOrientationPatient (0020:0037): ";
    std::cout << "LPS Orientation Matrix" << std::endl;
    std::cout << this->m_LPSDirCos << std::endl;

    this->m_SpacingMatrix.Fill(0.0);
    this->m_SpacingMatrix[0][0] = this->m_XRes;
    this->m_SpacingMatrix[1][1] = this->m_YRes;
    this->m_SpacingMatrix[2][2] = this->m_SliceSpacing;
    std::cout << "this->m_SpacingMatrix" << std::endl;
    std::cout << this->m_SpacingMatrix << std::endl;

    std::cout << "NRRDSpaceDirection" << std::endl;
    std::cout << this->GetNRRDSpaceDirection() << std::endl;

    }
  /** extract dwi data -- vendor specific so must happen in subclass
   *  implementing this method.
   */
  virtual void ExtractDWIData() = 0;

  /** access methods for image data */
  const DWIMetaDataDictionaryValidator::GradientTableType &GetDiffusionVectors() const { return this->m_DiffusionVectors; }

  const DWIMetaDataDictionaryValidator::GradientTableType computeScaledDiffusionVectors() const
  {
    //TODO: This can be further simplified.  Move private computeScaled(...) into this logic
    const DWIMetaDataDictionaryValidator::GradientTableType& UnitNormDiffusionVectors = this->GetDiffusionVectors();
    const std::vector<double>& bValues = this->GetBValues();
    const double maxBvalue = this->GetMaxBValue();
    const DWIMetaDataDictionaryValidator::GradientTableType BvalueScaledDiffusionVectors =
      this->computeScaledDiffusionVectors(UnitNormDiffusionVectors, bValues, maxBvalue);
    return BvalueScaledDiffusionVectors;
  }

  const std::vector<double> &GetBValues() const { return this->m_BValues; }
  double GetMaxBValue() const { return ComputeMaxBvalue( this->m_BValues); }
  void SetMeasurementFrameIdentity()
  {
    this->m_MeasurementFrame.SetIdentity();
  }

  VolumeType::Pointer GetDiffusionVolume() { return this->m_Volume; }

  SpacingType GetSpacing()
    {
      SpacingType spacing;
      spacing[0] = this->m_XRes;
      spacing[1] = this->m_YRes;
      spacing[2] = this->m_SliceSpacing;
      return spacing;
    }

  VolumeType::PointType GetOrigin() const
    {
      VolumeType::PointType rval;
      rval[0] = this->m_Origin[0];
      rval[1] = this->m_Origin[1];
      rval[2] = this->m_Origin[2];
      return rval;
    }

  RotationMatrixType   GetLPSDirCos() const { return this->m_LPSDirCos; }

  RotationMatrixType GetMeasurementFrame() const { return this->m_MeasurementFrame; }

  RotationMatrixType GetNRRDSpaceDirection() const { return  this->m_LPSDirCos * this->m_SpacingMatrix; }

  unsigned int GetNVolume() const { return this->m_NVolume; }

  std::string GetNRRDSpaceDefinition() const { return this->m_NRRDSpaceDefinition; }

  unsigned short GetRows() const { return m_Rows; }

  unsigned short GetCols() const { return m_Cols; }

  unsigned int GetSlicesPerVolume() const { return m_SlicesPerVolume; }

protected:
  /* determine if slice order is inferior to superior */
  void DetermineSliceOrderIS()
    {
      double image0Origin[3];
      image0Origin[0]=m_Origin[0];
      image0Origin[1]=m_Origin[1];
      image0Origin[2]=m_Origin[2];
      std::cout << "Slice 0: " << image0Origin[0] << " "
                << image0Origin[1] << " " << image0Origin[2] << std::endl;

      // assume volume interleaving, i.e. the second dicom file stores
      // the second slice in the same volume as the first dicom file
      double image1Origin[3];

      unsigned long nextSlice = 0;
      if (this->m_Headers.size() > 1)
        {
        // assuming multiple files is invalid for single-file volume: http://www.na-mic.org/Bug/view.php?id=4105
        nextSlice = this->m_IsInterleaved ? this->m_NVolume : 1;
        }

      this->m_Headers[nextSlice]->GetOrigin(image1Origin);
      std::cout << "Slice " << nextSlice << ": " << image1Origin[0] << " " << image1Origin[1]
                << " " << image1Origin[2] << std::endl;

      image1Origin[0] -= image0Origin[0];
      image1Origin[1] -= image0Origin[1];
      image1Origin[2] -= image0Origin[2];
      const DWIConverter::RotationMatrixType & NRRDSpaceDirection = this->GetNRRDSpaceDirection();
      double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2])
        + image1Origin[1] * (NRRDSpaceDirection[1][2])
        + image1Origin[2] * (NRRDSpaceDirection[2][2]);
      if( x1 < 0 )
        {
        this->m_SliceOrderIS = false;
        }
    }
  /** the SliceOrderIS flag can be computed (as above) but if it's
   *  invariant, the derived classes can just set the flag. This method
   *  fixes up the m_LPDirCos after the flag is set.
   */
  void SetDirectionsFromSliceOrder()
    {
      if(this->m_SliceOrderIS)
        {
        std::cout << "Slice order is IS" << std::endl;
        }
      else
        {
        std::cout << "Slice order is SI" << std::endl;
        this->m_LPSDirCos[0][2] = -this->m_LPSDirCos[0][2];
        this->m_LPSDirCos[1][2] = -this->m_LPSDirCos[1][2];
        this->m_LPSDirCos[2][2] = -this->m_LPSDirCos[2][2];
        }
    }

  /* given a sequence of dicom files where all the slices for location
   * 0 are folled by all the slices for location 1, etc. This method
   * transforms it into a sequence of volumes
   */
  void
  DeInterleaveVolume()
    {
      size_t NVolumes = this->m_NSlice / this->m_SlicesPerVolume;

      VolumeType::RegionType R = this->m_Volume->GetLargestPossibleRegion();

      R.SetSize(2, 1);
      std::vector<VolumeType::PixelType> v(this->m_NSlice);
      std::vector<VolumeType::PixelType> w(this->m_NSlice);

      itk::ImageRegionIteratorWithIndex<VolumeType> I(this->m_Volume, R );
      // permute the slices by extracting the 1D array of voxels for
      // a particular {x,y} position, then re-ordering the voxels such
      // that all the voxels for a particular volume are adjacent
      for( I.GoToBegin(); !I.IsAtEnd(); ++I )
        {
        VolumeType::IndexType idx = I.GetIndex();
        // extract all values in one "column"
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          idx[2] = k;
          v[k] = this->m_Volume->GetPixel( idx );
          }
        // permute
        for( unsigned int k = 0; k < NVolumes; ++k )
          {
          for( unsigned int m = 0; m < this->m_SlicesPerVolume; ++m )
            {
            w[(k * this->m_SlicesPerVolume) + m] = v[(m * NVolumes) + k];
            }
          }
        // put things back in order
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          idx[2] = k;
          this->m_Volume->SetPixel( idx, w[k] );
          }
        }
    }

  double ComputeMaxBvalue(const std::vector<double> &bValues) const
  {
    double maxBvalue(0.0);
    for( unsigned int k = 0; k < bValues.size(); ++k )
    {
      if( bValues[k] > maxBvalue )
      {
        maxBvalue = bValues[k];
      }
    }
    return maxBvalue;
  }

  DWIMetaDataDictionaryValidator::GradientTableType
  computeScaledDiffusionVectors( const DWIMetaDataDictionaryValidator::GradientTableType &UnitNormDiffusionVectors,
    const std::vector<double> &bValues,
    const double maxBvalue) const
  {
    DWIMetaDataDictionaryValidator::GradientTableType BvalueScaledDiffusionVectors;
    for( unsigned int k = 0; k < UnitNormDiffusionVectors.size(); ++k )
    {
      vnl_vector_fixed<double,3> vec(3);
      float scaleFactor = 0;
      if( maxBvalue > 0 )
      {
        scaleFactor = sqrt( bValues[k] / maxBvalue );
      }
      std::cout << "Scale Factor for Multiple BValues: " << k << " -- sqrt( " << bValues[k] << " / " << maxBvalue << " ) = "
                << scaleFactor << std::endl;
      for( unsigned ind = 0; ind < 3; ++ind )
      {
        vec[ind] = UnitNormDiffusionVectors[k][ind] * scaleFactor;
      }
      BvalueScaledDiffusionVectors.push_back(vec);
    }
    return BvalueScaledDiffusionVectors;
  }

  /** add vendor-specific flags; */
  virtual void AddFlagsToDictionary() = 0;
  /** one file reader per DICOM file in dataset */
  const DCMTKFileVector     m_Headers;
  /** the names of all the filenames, needed to use
   *  itk::ImageSeriesReader
   */
  const FileNamesContainer  m_InputFileNames;

  /** dimensions */
  unsigned short      m_Rows;
  unsigned short      m_Cols;
  unsigned int        m_SlicesPerVolume;

  /** spacing */
  double              m_XRes;
  double              m_YRes;
  double              m_SliceSpacing;

  /** image origin */
  PointType            m_Origin;

  /** measurement from for gradients if different than patient
   *  reference frame.
   */
  RotationMatrixType   m_MeasurementFrame;
  /** potentially the measurement frame */
  RotationMatrixType   m_LPSDirCos;
  /** matrix with just spacing information, used a couple places */
  RotationMatrixType   m_SpacingMatrix;
  /** the current dataset is represented in a single file */
  bool                m_MultiSliceVolume;
  /** slice order is inferior/superior? */
  bool                m_SliceOrderIS;
  /** the image read from the DICOM dataset */
  VolumeType::Pointer m_Volume;
  /** number of total slices */
  unsigned int        m_NSlice;
  /** number of gradient volumes */
  unsigned int        m_NVolume;

  /** list of B Values for each volume */
  std::vector<double>  m_BValues;
  /** list of gradient vectors */
  DWIMetaDataDictionaryValidator::GradientTableType  m_DiffusionVectors;
  /** double conversion instance, for optimal printing of numbers as
   *  text
   */
  itk::NumberToString<double> m_DoubleConvert;
  /** use the BMatrix to compute gradients in Siemens data instead of
   *  the reported graients. which are in many cases bogus.
   */
  const bool                  m_UseBMatrixGradientDirections;
  /** track if images is interleaved */
  bool                        m_IsInterleaved;
  /** again this is always the same (so far) but someone thought
   *  it might be important to change it.
   */
  std::string                 m_NRRDSpaceDefinition;
};

#endif // __DWIConverter_h
