//
// Created by Johnson, Hans J on 11/24/16.
//
#ifndef BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
#define BRAINSTOOLS_DWIDICOMCONVERTERBASE_H

#include <vector>
#include <iostream>
#include "DWIConverter.h"
#include "itkDCMTKSeriesFileNames.h"
#include "itkMacro.h"
#include "itkDCMTKImageIO.h"
#include "itkImage.h"
#include "itkDCMTKFileReader.h"
#include "itkNumberToString.h"
#include "DWIConvertUtils.h"

class DWIDICOMConverterBase : public DWIConverter {
 public:

  typedef itk::DCMTKSeriesFileNames           InputNamesGeneratorType;
  typedef std::vector<itk::DCMTKFileReader *> DCMTKFileVector;


  DWIDICOMConverterBase(const DCMTKFileVector &allHeaders,
               const FileNamesContainer &inputFileNames,
               const bool useBMatrixGradientDirections,
               const bool FSLFileFormatHorizontalBy3Rows) :
                                                    DWIConverter(inputFileNames, FSLFileFormatHorizontalBy3Rows),
                                                    m_UseBMatrixGradientDirections(useBMatrixGradientDirections),
                                                    m_Headers(allHeaders)

    {

    }


  /**
   * @brief Return common fields.  Does nothing for FSL
   * @return empty map
   */
  virtual CommonDicomFieldMapType GetCommonDicomFieldsMap() const ITK_OVERRIDE;

  virtual void LoadFromDisk() ITK_OVERRIDE {
    this->LoadDicomDirectory();
  }

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
          std::__1::cerr << "Exception thrown while reading DICOM volume"
                         << std::__1::endl;
          std::__1::cerr << excp << std::__1::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = false;
        }
      else
        {
        itk::ImageFileReader<Volume3DUnwrappedType>::Pointer reader =
          itk::ImageFileReader<Volume3DUnwrappedType>::New();
        reader->SetImageIO( dcmtkIO );
        reader->SetFileName( this->m_InputFileNames[0] );
        m_NSlice = this->m_InputFileNames.size();
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::__1::cerr << "Exception thrown while reading the series" << std::__1::endl;
          std::__1::cerr << excp << std::__1::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = true;
        }
      {
      // origin
      double origin[3];
      m_Headers[0]->GetOrigin(origin);
      Volume3DUnwrappedType::PointType imOrigin;
      imOrigin[0] = origin[0];
      imOrigin[1] = origin[1];
      imOrigin[2] = origin[2];
      this->m_Volume->SetOrigin(imOrigin);
      }
      // spacing
      {

      double spacing[3];
      m_Headers[0]->GetSpacing(spacing);
      SpacingType imSpacing;
      imSpacing[0] = spacing[0];
      imSpacing[1] = spacing[1];
      imSpacing[2] = spacing[2];
      m_Volume->SetSpacing(imSpacing);
      }

      // a map of ints keyed by the slice location string
      // reported in the dicom file.  The number of slices per
      // volume is the same as the number of unique slice locations
      std::__1::map<std::__1::string, int> sliceLocations;
      //
      // check for interleave
      if( !this->m_MultiSliceVolume )
        {
        // Make a hash of the sliceLocations in order to get the correct
        // count.  This is more reliable since SliceLocation may not be available.
        std::__1::vector<int>         sliceLocationIndicator;
        std::__1::vector<std::__1::string> sliceLocationStrings;

        sliceLocationIndicator.resize( this->m_NSlice );
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          std::__1::string originString;
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
        std::__1::cout << "=================== this->m_SlicesPerVolume:" << this->m_SlicesPerVolume << std::__1::endl;


        // if the this->m_SlicesPerVolume == 1, de-interleaving won't do
        // anything so there's no point in doing it.
        if( this->m_NSlice >= 2 && this->m_SlicesPerVolume > 1 )
          {
          if( sliceLocationIndicator[0] != sliceLocationIndicator[1] )
            {
            std::__1::cout << "Dicom images are ordered in a volume interleaving way." << std::__1::endl;
            }
          else
            {
            std::__1::cout << "Dicom images are ordered in a slice interleaving way." << std::__1::endl;
            this->m_IsInterleaved = true;
            // reorder slices into a volume interleaving manner
            DeInterleaveVolume();
            }
          }
        }

    {
    Volume3DUnwrappedType::DirectionType LPSDirCos;
    LPSDirCos.SetIdentity();

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
        LPSDirCos[j][i] = *dirCosArrayP;
        }
      }

    // Cross product, this gives I-axis direction
    LPSDirCos[0][2] = LPSDirCos[1][0] * LPSDirCos[2][1] - LPSDirCos[2][0] * LPSDirCos[1][1];
    LPSDirCos[1][2] = LPSDirCos[2][0] * LPSDirCos[0][1] - LPSDirCos[0][0] * LPSDirCos[2][1];
    LPSDirCos[2][2] = LPSDirCos[0][0] * LPSDirCos[1][1] - LPSDirCos[1][0] * LPSDirCos[0][1];

    this->m_Volume->SetDirection(LPSDirCos);
    }
    std::__1::cout << "ImageOrientationPatient (0020:0037): ";
    std::__1::cout << "LPS Orientation Matrix" << std::__1::endl;
    std::__1::cout << this->m_Volume->GetDirection() << std::__1::endl;


    //TODO: Remove __1
    std::__1::cout << "this->m_SpacingMatrix" << std::__1::endl;
    std::__1::cout << this->GetSpacingMatrix() << std::__1::endl;

    std::__1::cout << "NRRDSpaceDirection" << std::__1::endl;
    std::__1::cout << this->GetNRRDSpaceDirection() << std::__1::endl;
      //TODO: Add metadata to the DWI images
      {
        //<element tag="0008,0060" vr="CS" vm="1" len="2" name="Modality">MR</element>
        this->_addToStringDictionary("0008","0060","Modality",DCM_CS);
        //<element tag="0008,0070" vr="LO" vm="1" len="18" name="Manufacturer">GE MEDICAL SYSTEMS</element>
        this->_addToStringDictionary("0008","0070","Manufacturer",DCM_LO);
        //<element tag="0008,1090" vr="LO" vm="1" len="10" name="ManufacturerModelName">SIGNA HDx</element>
        this->_addToStringDictionary("0008","1090","ManufacturerModelName",DCM_LO);
        //<element tag="0018,0087" vr="DS" vm="1" len="2" name="MagneticFieldStrength">3</element>
        this->_addToStringDictionary("0018","0087","MagneticFieldStrength",DCM_DS);
        //<element tag="0018,1020" vr="LO" vm="3" len="42" name="SoftwareVersions">14\LX\MR Software release:14.0_M5A_0828.b</element>
        this->_addToStringDictionary("0018","1020","SoftwareVersions",DCM_LO);
        //<element tag="0018,0022" vr="CS" vm="2" len="12" name="ScanOptions">EPI_GEMS\PFF</element>
        this->_addToStringDictionary("0018","0022","ScanOptions",DCM_CS);
        //<element tag="0018,0023" vr="CS" vm="1" len="2" name="MRAcquisitionType">2D</element>
        this->_addToStringDictionary("0018","0023","MRAcquisitionType",DCM_CS);
        //<element tag="0018,0080" vr="DS" vm="1" len="6" name="RepetitionTime">12000</element>
        this->_addToStringDictionary("0018","0080","RepetitionTime",DCM_DS);
        //<element tag="0018,0081" vr="DS" vm="1" len="4" name="EchoTime">74.7</element>
        this->_addToStringDictionary("0018","0081","EchoTime",DCM_DS);
        //<element tag="0018,0083" vr="DS" vm="1" len="2" name="NumberOfAverages">1</element>
        this->_addToStringDictionary("0018","0083","NumberOfAverages",DCM_DS);
        //<element tag="0018,1314" vr="DS" vm="1" len="2" name="FlipAngle">90</element>
        this->_addToStringDictionary("0018","1314","FlipAngle",DCM_DS);
      }
    }


protected:
  enum VRType{
    DCM_CS,
    DCM_LO,
    DCM_SH,
    DCM_DS,
  };

  /**
   * @brief Create a dictionary of dicom extracted information about the scans
   *   using dcmtk dcm2xml tool, identify the information desired to be kept
   *   <element tag="0018,1314" vr="DS" vm="1" len="2" name="FlipAngle">90</element>
   * @param dcm_primary_name "0018" in example above
   * @param dcm_seconary_name "1314" in example above
   * @param dcm_human_readable_name "FlipAngle" in example above
   * @param vr "DCM_DS" for enumeration as indicated by vr in example above
   */
  void _addToStringDictionary(const std::string dcm_primary_name,
    const std::string dcm_seconary_name,
    const std::string dcm_human_readable_name, const enum VRType vr)
  {
    int dcm_primary_code;
    {
      std::istringstream iss(dcm_primary_name);
      iss >> std::hex >> dcm_primary_code;
    }
    int dcm_secondary_code;
    {
      std::istringstream iss(dcm_seconary_name);
      iss >> std::hex >> dcm_secondary_code;
    }
    std::string stringValue="UNKNOWN";
    const bool throwException=false;
    switch(vr)
    {
    case DCM_CS:
      this->m_Headers[0]->GetElementCS(dcm_primary_code, dcm_secondary_code, stringValue, throwException );
      break;
    case DCM_LO:
      this->m_Headers[0]->GetElementLO(dcm_primary_code, dcm_secondary_code, stringValue, throwException );
      break;
    case DCM_SH:
      this->m_Headers[0]->GetElementLO(dcm_primary_code, dcm_secondary_code, stringValue, throwException );
      break;
    case DCM_DS:
      this->m_Headers[0]->GetElementDS(dcm_primary_code, dcm_secondary_code, stringValue, throwException );
      break;
    default:
      stringValue="INVALIDDR";
    }
    // in NRRD key name DICOM_0008_0060_Modality:=MR
    std::string map_name("DICOM_");
    map_name+=dcm_primary_name+"_"+dcm_seconary_name+"_"+dcm_human_readable_name;
    this->m_CommonDicomFieldsMap[map_name] = stringValue;
  }

  /* determine if slice order is inferior to superior */
  void DetermineSliceOrderIS()
    {
      double image0Origin[3];
      image0Origin[0]=this->m_Volume->GetOrigin()[0];
      image0Origin[1]=this->m_Volume->GetOrigin()[1];
      image0Origin[2]=this->m_Volume->GetOrigin()[2];
      std::__1::cout << "Slice 0: " << image0Origin[0] << " "
                     << image0Origin[1] << " " << image0Origin[2] << std::__1::endl;

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
      std::__1::cout << "Slice " << nextSlice << ": " << image1Origin[0] << " " << image1Origin[1]
                     << " " << image1Origin[2] << std::__1::endl;

      image1Origin[0] -= image0Origin[0];
      image1Origin[1] -= image0Origin[1];
      image1Origin[2] -= image0Origin[2];
      const RotationMatrixType & NRRDSpaceDirection = this->GetNRRDSpaceDirection();
      double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2])
        + image1Origin[1] * (NRRDSpaceDirection[1][2])
        + image1Origin[2] * (NRRDSpaceDirection[2][2]);
      if( x1 < 0 )
        {
        this->m_SliceOrderIS = false;
        }
    }

  /** force use of the BMatrix to compute gradients in Siemens data instead of
   *  the reported gradients. which are in many cases bogus.
   */
  const bool m_UseBMatrixGradientDirections;
  /** one file reader per DICOM file in dataset */
  const DCMTKFileVector    m_Headers;
};

#endif //BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
