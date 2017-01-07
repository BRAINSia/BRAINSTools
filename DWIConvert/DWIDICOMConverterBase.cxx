//
// Created by Johnson, Hans J on 11/24/16.
//

#include "DWIDICOMConverterBase.h"

/**
 * @brief Return common fields.  Does nothing for FSL
 * @return empty map
 */
DWIDICOMConverterBase::CommonDicomFieldMapType
DWIDICOMConverterBase::GetCommonDicomFieldsMap() const
{
return this->m_CommonDicomFieldsMap;
}

DWIDICOMConverterBase::DWIDICOMConverterBase(const DCMTKFileVector &allHeaders,
                                             const FileNamesContainer &inputFileNames,
                                             const bool useBMatrixGradientDirections,
                                             const bool FSLFileFormatHorizontalBy3Rows) :
        DWIConverter(inputFileNames, FSLFileFormatHorizontalBy3Rows),
        m_UseBMatrixGradientDirections(useBMatrixGradientDirections),
        m_Headers(allHeaders),
        m_MultiSliceVolume(false),
        m_SliceOrderIS(true),
        m_IsInterleaved(false),
        m_SlicesPerVolume(0),
        m_NSlice(0),
        m_NVolume(0)
{

}

void DWIDICOMConverterBase::LoadDicomDirectory()
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
      m_3DUnwrappedVolume = reader->GetOutput();
    m_MultiSliceVolume = false;
  }
  else
  {
    itk::ImageFileReader<Volume3DType>::Pointer reader =
            itk::ImageFileReader<Volume3DType>::New();
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
      m_3DUnwrappedVolume = reader->GetOutput();
    m_MultiSliceVolume = true;
  }
  {
    // origin
    double origin[3];
    m_Headers[0]->GetOrigin(origin);
    Volume3DType::PointType imOrigin;
    imOrigin[0] = origin[0];
    imOrigin[1] = origin[1];
    imOrigin[2] = origin[2];
    this->m_3DUnwrappedVolume->SetOrigin(imOrigin);
  }
  // spacing
  {

    double spacing[3];
    m_Headers[0]->GetSpacing(spacing);
    SpacingType imSpacing;
    imSpacing[0] = spacing[0];
    imSpacing[1] = spacing[1];
    imSpacing[2] = spacing[2];
      m_3DUnwrappedVolume->SetSpacing(imSpacing);
  }

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

  {
    Volume3DType::DirectionType LPSDirCos;
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

    this->m_3DUnwrappedVolume->SetDirection(LPSDirCos);
  }
  std::cout << "ImageOrientationPatient (0020:0037): ";
  std::cout << "LPS Orientation Matrix" << std::endl;
  std::cout << this->m_3DUnwrappedVolume->GetDirection() << std::endl;

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

void DWIDICOMConverterBase::LoadFromDisk()
{
  this->LoadDicomDirectory();


}


/**
 * @brief Create a dictionary of dicom extracted information about the scans
 *   using dcmtk dcm2xml tool, identify the information desired to be kept
 *   <element tag="0018,1314" vr="DS" vm="1" len="2" name="FlipAngle">90</element>
 * @param dcm_primary_name "0018" in example above
 * @param dcm_seconary_name "1314" in example above
 * @param dcm_human_readable_name "FlipAngle" in example above
 * @param vr "DCM_DS" for enumeration as indicated by vr in example above
 */
void DWIDICOMConverterBase::_addToStringDictionary(const std::string dcm_primary_name,
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

/** the SliceOrderIS flag can be computed (as above) but if it's
 *  invariant, the derived classes can just set the flag. This method
 *  fixes up the VolumeDirectionCos after the flag is set.
 */
void DWIDICOMConverterBase::SetDirectionsFromSliceOrder()
{
  if(this->m_SliceOrderIS)
  {
    std::cout << "Slice order is IS" << std::endl;
  }
  else
  {
    std::cout << "Slice order is SI" << std::endl;
    Volume3DType::DirectionType LPSDirCos = this->m_3DUnwrappedVolume->GetDirection();
    LPSDirCos[0][2] *= -1;
    LPSDirCos[1][2] *= -1;
    LPSDirCos[2][2] *= -1;
    this->m_3DUnwrappedVolume->SetDirection(LPSDirCos);
    //Need to update the measurement frame too!
    this->m_MeasurementFrame[0][2] *= -1;
    this->m_MeasurementFrame[1][2] *= -1;
    this->m_MeasurementFrame[2][2] *= -1;
  }
}

/* given a sequence of dicom files where all the slices for location
 * 0 are folled by all the slices for location 1, etc. This method
 * transforms it into a sequence of volumes
 */
void DWIDICOMConverterBase::DeInterleaveVolume()
{
  size_t NVolumes = this->m_NSlice / this->m_SlicesPerVolume;

  Volume3DType::RegionType R = this->m_3DUnwrappedVolume->GetLargestPossibleRegion();

  R.SetSize(2, 1);
  std::vector<Volume3DType::PixelType> v(this->m_NSlice);
  std::vector<Volume3DType::PixelType> w(this->m_NSlice);

  itk::ImageRegionIteratorWithIndex<Volume3DType> I(this->m_3DUnwrappedVolume, R );
  // permute the slices by extracting the 1D array of voxels for
  // a particular {x,y} position, then re-ordering the voxels such
  // that all the voxels for a particular volume are adjacent
  for( I.GoToBegin(); !I.IsAtEnd(); ++I )
  {
    Volume3DType::IndexType idx = I.GetIndex();
    // extract all values in one "column"
    for( unsigned int k = 0; k < this->m_NSlice; ++k )
    {
      idx[2] = k;
      v[k] = this->m_3DUnwrappedVolume->GetPixel( idx );
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
      this->m_3DUnwrappedVolume->SetPixel( idx, w[k] );
    }
  }
}
/* determine if slice order is inferior to superior */
void DWIDICOMConverterBase::DetermineSliceOrderIS()
{
  double image0Origin[3];
  image0Origin[0]=this->m_3DUnwrappedVolume->GetOrigin()[0];
  image0Origin[1]=this->m_3DUnwrappedVolume->GetOrigin()[1];
  image0Origin[2]=this->m_3DUnwrappedVolume->GetOrigin()[2];
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
  const RotationMatrixType & NRRDSpaceDirection = GetNRRDSpaceDirection<Volume3DType>(this->m_3DUnwrappedVolume);
  double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2])
              + image1Origin[1] * (NRRDSpaceDirection[1][2])
              + image1Origin[2] * (NRRDSpaceDirection[2][2]);
  if( x1 < 0 )
  {
    this->m_SliceOrderIS = false;
  }
}


Volume4DType::Pointer DWIDICOMConverterBase::ThreeDUnwrappedToFourDImage(Volume3DType::Pointer img) const
{
  const int nVolumes = this->GetNVolume();

  Volume3DType::SizeType size3D(img->GetLargestPossibleRegion().GetSize());
  Volume3DType::DirectionType direction3D(img->GetDirection());
  Volume3DType::SpacingType spacing3D(img->GetSpacing());
  Volume3DType::PointType origin3D(img->GetOrigin());

  Volume4DType::RegionType region4D;
  {
    Volume4DType::SizeType size4D;
    size4D[0] = size3D[0];
    size4D[1] = size3D[1];
    size4D[2] = size3D[2]/nVolumes;
    size4D[3] = nVolumes;
    Volume4DType::IndexType index4D;
    index4D.Fill(0);
    region4D.SetIndex(index4D);
    region4D.SetSize(size4D);

    if ((size4D[2]*nVolumes)!=size3D[2]) {
      itkGenericExceptionMacro(
              << "#of slices in volume not evenly divisible by"
              << " the number of volumes: slices = " << size3D[2]
              << " volumes = " << nVolumes << " left-over slices = "
              << size3D[2] % nVolumes << std::endl);
    }
  }
  Volume4DType::DirectionType direction4D;
  direction4D.SetIdentity();
  Volume4DType::SpacingType   spacing4D;
  spacing4D.Fill(1.0);
  Volume4DType::PointType     origin4D;
  origin4D.Fill(0.0);
  for( unsigned i = 0; i < 3; ++i )
  {
    for( unsigned j = 0; j < 3; ++j )
    {
      direction4D[i][j] = direction3D[i][j];
    }
    spacing4D[i] = spacing3D[i];
    origin4D[i] = origin3D[i];
  }


  Volume4DType::Pointer img4D = Volume4DType::New();
  img4D->SetRegions(region4D);
  img4D->SetDirection(direction4D);
  img4D->SetSpacing(spacing4D);
  img4D->SetOrigin(origin4D);
  img4D->Allocate();

  {
    img4D->SetMetaDataDictionary(img->GetMetaDataDictionary());
    //Set the qform and sfrom codes for the MetaDataDictionary.
    itk::MetaDataDictionary & thisDic = img4D->GetMetaDataDictionary();
    itk::EncapsulateMetaData< std::string >( thisDic, "qform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
    itk::EncapsulateMetaData< std::string >( thisDic, "sform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
  }

  const size_t bytecount = img4D->GetLargestPossibleRegion().GetNumberOfPixels()
                           * sizeof(PixelValueType);

  memcpy(img4D->GetBufferPointer(), img->GetBufferPointer(), bytecount);
  return img4D;
}


Volume3DType::Pointer DWIDICOMConverterBase::FourDToThreeDUnwrappedImage(Volume4DType::Pointer img4D) const
{
  Volume4DType::SizeType      size4D(img4D->GetLargestPossibleRegion().GetSize() );
  Volume4DType::DirectionType direction4D(img4D->GetDirection() );
  Volume4DType::SpacingType   spacing4D(img4D->GetSpacing() );
  Volume4DType::PointType     origin4D(img4D->GetOrigin() );

  Volume3DType::RegionType region3D;
  {
    Volume3DType::SizeType size3D;
    size3D[0] = size4D[0];
    size3D[1] = size4D[1];
    const int nVolumes = img4D->GetLargestPossibleRegion().GetSize()[3];
    size3D[2] = size4D[2] * nVolumes;

    Volume3DType::IndexType index3D;
    index3D.Fill(0);
    region3D.SetIndex(index3D);
    region3D.SetSize(size3D);

    if( (size4D[2] * nVolumes) != size3D[2] )
    {
      itkGenericExceptionMacro(
              << "#of slices in volume not evenly divisible by"
              << " the number of volumes: slices = " << size3D[2]
              << " volumes = " << nVolumes << " left-over slices = "
              << size3D[2] % nVolumes << std::endl);
    }
  }
  Volume3DType::DirectionType direction3D;
  direction3D.SetIdentity();
  Volume3DType::SpacingType   spacing3D;
  spacing3D.Fill(1.0);
  Volume3DType::PointType     origin3D;
  origin3D.Fill(0.0);
  for( unsigned i = 0; i < 3; ++i )
  {
    for( unsigned j = 0; j < 3; ++j )
    {
      direction3D[i][j] = direction4D[i][j];
    }
    spacing3D[i] = spacing4D[i];
    origin3D[i] = origin4D[i];
  }

  Volume3DType::Pointer img = Volume3DType::New();
  img->SetRegions(region3D);
  img->SetDirection(direction3D);
  img->SetSpacing(spacing3D);
  img->SetOrigin(origin3D);
  img->Allocate();

  {
    img->SetMetaDataDictionary(img4D->GetMetaDataDictionary());
    //Set the qform and sfrom codes for the MetaDataDictionary.
    itk::MetaDataDictionary & thisDic = img->GetMetaDataDictionary();
    itk::EncapsulateMetaData< std::string >( thisDic, "qform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
    itk::EncapsulateMetaData< std::string >( thisDic, "sform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
  }

  const size_t bytecount = img->GetLargestPossibleRegion().GetNumberOfPixels()
                           * sizeof(PixelValueType);

  memcpy(img->GetBufferPointer(), img4D->GetBufferPointer(), bytecount);
  return img;
}

unsigned int DWIDICOMConverterBase::GetSlicesPerVolume() const { return m_SlicesPerVolume; }
unsigned int DWIDICOMConverterBase::GetNVolume() const { return this->m_NVolume; }