//
// Created by Johnson, Hans J on 11/24/16.
//

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>


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

DWIDICOMConverterBase::DWIDICOMConverterBase(DCMTKFileVector            allHeaders,
                                             const FileNamesContainer & inputFileNames,
                                             const bool                 useBMatrixGradientDirections)
  : DWIConverter(inputFileNames)
  , m_UseBMatrixGradientDirections(useBMatrixGradientDirections)
  , m_Headers(std::move(allHeaders))
  , m_MultiSliceVolume(false)
  , m_SliceOrderIS(true)
{}

namespace
{
constexpr double samePositionTolerance = 0.01; // mm along the slice normal

std::vector<double>
SliceDistances(const DWIDICOMConverterBase::DCMTKFileVector & headers)
{
  double dircos[6] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
  headers[0]->GetDirCosArray(dircos);
  const double        normal[3] = { dircos[1] * dircos[5] - dircos[2] * dircos[4],
                                    dircos[2] * dircos[3] - dircos[0] * dircos[5],
                                    dircos[0] * dircos[4] - dircos[1] * dircos[3] };
  std::vector<double> distances;
  distances.reserve(headers.size());
  for (const auto & header : headers)
  {
    double origin[3] = { 0.0, 0.0, 0.0 };
    header->GetOrigin(origin);
    distances.push_back(normal[0] * origin[0] + normal[1] * origin[1] + normal[2] * origin[2]);
  }
  return distances;
}
} // namespace

void
DWIDICOMConverterBase::OrderSlicesVolumeMajor(DCMTKFileVector & headers)
{
  const size_t nFiles = headers.size();
  if (nFiles < 2)
  {
    return;
  }
  const std::vector<double> distances = SliceDistances(headers);
  std::vector<itk::int32_t> instanceNumbers(nFiles, 0);
  for (size_t i = 0; i < nFiles; ++i)
  {
    headers[i]->GetElementIS(0x0020, 0x0013, instanceNumbers[i], false);
  }
  const auto acquiredBefore = [&](size_t a, size_t b) {
    if (instanceNumbers[a] != instanceNumbers[b])
    {
      return instanceNumbers[a] < instanceNumbers[b];
    }
    return headers[a]->GetFileName() < headers[b]->GetFileName();
  };

  std::vector<size_t> byDistance(nFiles);
  std::iota(byDistance.begin(), byDistance.end(), 0);
  std::stable_sort(
    byDistance.begin(), byDistance.end(), [&](size_t a, size_t b) { return distances[a] < distances[b]; });
  std::vector<std::vector<size_t>> positions;
  for (size_t k = 0; k < nFiles; ++k)
  {
    const size_t i = byDistance[k];
    if (positions.empty() || distances[i] - distances[positions.back().front()] > samePositionTolerance)
    {
      positions.emplace_back();
    }
    positions.back().push_back(i);
  }
  const size_t nVolumes = nFiles / positions.size();
  for (auto & position : positions)
  {
    if (position.size() != nVolumes || nVolumes * positions.size() != nFiles)
    {
      itkGenericExceptionMacro(<< "Missing DICOM slice files: " << nFiles << " files at " << positions.size()
                               << " slice positions do not repeat the same number of times");
    }
    std::sort(position.begin(), position.end(), acquiredBefore);
  }

  std::stable_sort(
    positions.begin(), positions.end(), [&](const std::vector<size_t> & a, const std::vector<size_t> & b) {
      return acquiredBefore(a.front(), b.front());
    });
  bool monotonic = true;
  for (size_t s = 2; s < positions.size(); ++s)
  {
    const double previousStep = distances[positions[s - 1].front()] - distances[positions[s - 2].front()];
    const double step = distances[positions[s].front()] - distances[positions[s - 1].front()];
    monotonic = monotonic && (previousStep * step > 0.0);
  }
  if (!monotonic)
  {
    std::sort(positions.begin(), positions.end(), [&](const std::vector<size_t> & a, const std::vector<size_t> & b) {
      return distances[a.front()] < distances[b.front()];
    });
  }

  DCMTKFileVector ordered;
  ordered.reserve(nFiles);
  for (size_t v = 0; v < nVolumes; ++v)
  {
    for (const auto & position : positions)
    {
      ordered.push_back(std::move(headers[position[v]]));
    }
  }
  headers = std::move(ordered);
}

void
DWIDICOMConverterBase::ValidateVolumeEncodings(
  const std::function<void(unsigned int, double &, vnl_vector_fixed<double, 3> &)> & readSliceEncoding) const
{
  constexpr double bTolerance = 1.0;
  constexpr double gradientTolerance = 1e-4;
  for (unsigned int first = 0; first + this->m_SlicesPerVolume <= this->m_NSlice; first += this->m_SlicesPerVolume)
  {
    double                      firstB = 0.0;
    vnl_vector_fixed<double, 3> firstGradient(0.0);
    readSliceEncoding(first, firstB, firstGradient);
    for (unsigned int k = first + 1; k < first + this->m_SlicesPerVolume; ++k)
    {
      double                      b = 0.0;
      vnl_vector_fixed<double, 3> gradient(0.0);
      readSliceEncoding(k, b, gradient);
      const bool bothB0 = b < bTolerance && firstB < bTolerance;
      if (std::abs(b - firstB) > bTolerance || (!bothB0 && (gradient - firstGradient).inf_norm() > gradientTolerance))
      {
        itkGenericExceptionMacro(<< "Volume " << first / this->m_SlicesPerVolume
                                 << " mixes diffusion encodings: " << this->m_Headers[first]->GetFileName()
                                 << " has b=" << firstB << " gradient " << firstGradient << " but "
                                 << this->m_Headers[k]->GetFileName() << " has b=" << b << " gradient " << gradient);
      }
    }
  }
}

void
DWIDICOMConverterBase::LoadDicomDirectory()
{
  //
  // add vendor-specific flags to dictionary
  this->AddFlagsToDictionary();
  //
  // load the volume, either single or multivolume.
  m_NSlice = this->m_InputFileNames.size();
  auto dcmtkIO = itk::DCMTKImageIO::New();
  if (this->m_InputFileNames.size() > 1)
  {
    auto reader = ReaderType::New();
    reader->SetImageIO(dcmtkIO);
    reader->SetFileNames(this->m_InputFileNames);
    try
    {
      reader->Update();
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Exception thrown while reading DICOM volume" << std::endl;
      std::cerr << excp << std::endl;
      throw;
    }
    m_Volume = reader->GetOutput();
    m_MultiSliceVolume = false;
  }
  else
  {
    auto reader = itk::ImageFileReader<Volume3DUnwrappedType>::New();
    reader->SetImageIO(dcmtkIO);
    reader->SetFileName(this->m_InputFileNames[0]);
    m_NSlice = this->m_InputFileNames.size();
    try
    {
      reader->Update();
    }
    catch (const itk::ExceptionObject & excp)
    {
      std::cerr << "Exception thrown while reading the series" << std::endl;
      std::cerr << excp << std::endl;
      throw;
    }
    m_Volume = reader->GetOutput();
    m_MultiSliceVolume = true;
  }
  {
    // origin
    double origin[3];
    m_Headers[0]->GetOrigin(origin);
    Volume3DUnwrappedType::PointType imOrigin(origin);
    this->m_Volume->SetOrigin(imOrigin);
  }
  // spacing
  {

    double spacing[3];
    // m_Headers[0]->GetSpacing(spacing);
    getDicomSpacing(spacing);
    SpacingType imSpacing(spacing);
    m_Volume->SetSpacing(imSpacing);
  }
  m_thickness = readThicknessFromDicom();

  // Unsigned pixels (#457) are stored in the short volume with their bits intact and written as unsigned short.
  {
    unsigned short pixelRep = 1; // absent means signed
    m_Headers[0]->GetElementUS(0x0028, 0x0103, pixelRep, false);
    double rescaleIntercept = 0.0;
    m_Headers[0]->GetElementDS<double>(0x0028, 0x1052, 1, &rescaleIntercept, false);
    this->m_UnsignedPixelData = pixelRep == 0 && rescaleIntercept >= 0.0;
  }

  if (!this->m_MultiSliceVolume)
  {
    // Headers arrive volume by volume (OrderSlicesVolumeMajor); a volume ends where the first position repeats.
    const std::vector<double> distances = SliceDistances(this->m_Headers);
    unsigned int              slicesPerVolume = 1;
    while (slicesPerVolume < this->m_NSlice &&
           std::abs(distances[slicesPerVolume] - distances[0]) > samePositionTolerance)
    {
      ++slicesPerVolume;
    }
    if (this->m_NSlice % slicesPerVolume != 0)
    {
      itkGenericExceptionMacro(<< "Missing DICOM Slice files: Number of slice files (" << this->m_NSlice
                               << ") not evenly divisible by the number of slice locations (" << slicesPerVolume
                               << ")");
    }
    this->m_SlicesPerVolume = slicesPerVolume;
    std::cout << "=================== this->m_SlicesPerVolume:" << this->m_SlicesPerVolume << std::endl;
  }

  {
    Volume3DUnwrappedType::DirectionType LPSDirCos = Volume3DUnwrappedType::DirectionType::GetIdentity();

    // check ImageOrientationPatient and figure out slice direction in
    // L-P-I (right-handed) system.
    // In Dicom, the coordinate frame is L-P by default. Look at
    // http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301
    double dirCosArray[6];
    // 0020,0037 -- Image Orientation (Patient)
    this->m_Headers[0]->GetDirCosArray(dirCosArray);
    double * dirCosArrayP = dirCosArray;
    for (unsigned i = 0; i < 2; ++i)
    {
      for (unsigned j = 0; j < 3; ++j, ++dirCosArrayP)
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
  std::cout << "ImageOrientationPatient (0020:0037): ";
  std::cout << "LPS Orientation Matrix" << std::endl;
  std::cout << this->m_Volume->GetDirection() << std::endl;


  // INFO: Remove __1
  std::cout << "this->m_SpacingMatrix" << std::endl;
  std::cout << this->GetSpacingMatrix() << std::endl;

  std::cout << "NRRDSpaceDirection" << std::endl;
  std::cout << this->GetNRRDSpaceDirection() << std::endl;
  // INFO: Add metadata to the DWI images
  {
    //<element tag="0008,0060" vr="CS" vm="1" len="2" name="Modality">MR</element>
    this->_addToStringDictionary("0008", "0060", "Modality", DCM_CS);
    //<element tag="0008,0070" vr="LO" vm="1" len="18" name="Manufacturer">GE MEDICAL SYSTEMS</element>
    this->_addToStringDictionary("0008", "0070", "Manufacturer", DCM_LO);
    //<element tag="0008,1090" vr="LO" vm="1" len="10" name="ManufacturerModelName">SIGNA HDx</element>
    this->_addToStringDictionary("0008", "1090", "ManufacturerModelName", DCM_LO);
    //<element tag="0018,0087" vr="DS" vm="1" len="2" name="MagneticFieldStrength">3</element>
    this->_addToStringDictionary("0018", "0087", "MagneticFieldStrength", DCM_DS);
    //<element tag="0018,1020" vr="LO" vm="3" len="42" name="SoftwareVersions">14\LX\MR Software
    // release:14.0_M5A_0828.b</element>
    this->_addToStringDictionary("0018", "1020", "SoftwareVersions", DCM_LO);
    //<element tag="0018,0022" vr="CS" vm="2" len="12" name="ScanOptions">EPI_GEMS\PFF</element>
    this->_addToStringDictionary("0018", "0022", "ScanOptions", DCM_CS);
    //<element tag="0018,0023" vr="CS" vm="1" len="2" name="MRAcquisitionType">2D</element>
    this->_addToStringDictionary("0018", "0023", "MRAcquisitionType", DCM_CS);
    //<element tag="0018,0080" vr="DS" vm="1" len="6" name="RepetitionTime">12000</element>
    this->_addToStringDictionary("0018", "0080", "RepetitionTime", DCM_DS);
    //<element tag="0018,0081" vr="DS" vm="1" len="4" name="EchoTime">74.7</element>
    this->_addToStringDictionary("0018", "0081", "EchoTime", DCM_DS);
    //<element tag="0018,0083" vr="DS" vm="1" len="2" name="NumberOfAverages">1</element>
    this->_addToStringDictionary("0018", "0083", "NumberOfAverages", DCM_DS);
    //<element tag="0018,1314" vr="DS" vm="1" len="2" name="FlipAngle">90</element>
    this->_addToStringDictionary("0018", "1314", "FlipAngle", DCM_DS);
  }
}

void
DWIDICOMConverterBase::LoadFromDisk()
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
void
DWIDICOMConverterBase::_addToStringDictionary(const std::string & dcm_primary_name,
                                              const std::string & dcm_seconary_name,
                                              const std::string & dcm_human_readable_name,
                                              const enum VRType   vr)
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
  std::string stringValue = "UNKNOWN";
  const bool  throwException = false;
  switch (vr)
  {
    case DCM_CS:
      this->m_Headers[0]->GetElementCS(dcm_primary_code, dcm_secondary_code, stringValue, throwException);
      break;
    case DCM_LO:
    case DCM_SH:
      this->m_Headers[0]->GetElementLO(dcm_primary_code, dcm_secondary_code, stringValue, throwException);
      break;
    case DCM_DS:
      this->m_Headers[0]->GetElementDS(dcm_primary_code, dcm_secondary_code, stringValue, throwException);
      break;
    default:
      stringValue = "INVALIDDR";
  }
  // in NRRD key name DICOM_0008_0060_Modality:=MR
  std::string map_name("DICOM_");
  map_name += dcm_primary_name + "_" + dcm_seconary_name + "_" + dcm_human_readable_name;
  this->m_CommonDicomFieldsMap[map_name] = stringValue;
}

/** the SliceOrderIS flag can be computed (as above) but if it's
 *  invariant, the derived classes can just set the flag. This method
 *  fixes up the VolumeDirectionCos after the flag is set.
 */
void
DWIDICOMConverterBase::SetDirectionsFromSliceOrder()
{
  if (this->m_SliceOrderIS)
  {
    std::cout << "Slice order is IS" << std::endl;
  }
  else
  {
    std::cout << "Slice order is SI" << std::endl;
    Volume3DUnwrappedType::DirectionType LPSDirCos = this->m_Volume->GetDirection();
    LPSDirCos[0][2] *= -1;
    LPSDirCos[1][2] *= -1;
    LPSDirCos[2][2] *= -1;
    this->m_Volume->SetDirection(LPSDirCos);
    // Need to update the measurement frame too!
    this->m_MeasurementFrame[0][2] *= -1;
    this->m_MeasurementFrame[1][2] *= -1;
    this->m_MeasurementFrame[2][2] *= -1;
  }
}

/* given a sequence of dicom files where all the slices for location
 * 0 are folled by all the slices for location 1, etc. This method
 * transforms it into a sequence of volumes
 */
void
DWIDICOMConverterBase::DeInterleaveVolume()
{
  size_t NVolumes = this->m_NSlice / this->m_SlicesPerVolume;

  Volume3DUnwrappedType::RegionType R = this->m_Volume->GetLargestPossibleRegion();

  R.SetSize(2, 1);
  std::vector<Volume3DUnwrappedType::PixelType> v(this->m_NSlice);
  std::vector<Volume3DUnwrappedType::PixelType> w(this->m_NSlice);

  itk::ImageRegionConstIteratorWithIndex<Volume3DUnwrappedType> I(this->m_Volume, R);
  // permute the slices by extracting the 1D array of voxels for
  // a particular {x,y} position, then re-ordering the voxels such
  // that all the voxels for a particular volume are adjacent
  for (I.GoToBegin(); !I.IsAtEnd(); ++I)
  {
    Volume3DUnwrappedType::IndexType idx = I.GetIndex();
    // extract all values in one "column"
    for (unsigned int k = 0; k < this->m_NSlice; ++k)
    {
      idx[2] = k;
      v[k] = this->m_Volume->GetPixel(idx);
    }
    // permute
    for (unsigned int k = 0; k < NVolumes; ++k)
    {
      for (unsigned int m = 0; m < this->m_SlicesPerVolume; ++m)
      {
        w[(k * this->m_SlicesPerVolume) + m] = v[(m * NVolumes) + k];
      }
    }
    // put things back in order
    for (unsigned int k = 0; k < this->m_NSlice; ++k)
    {
      idx[2] = k;
      this->m_Volume->SetPixel(idx, w[k]);
    }
  }
}
/* determine if slice order is inferior to superior */
void
DWIDICOMConverterBase::DetermineSliceOrderIS()
{
  double image0Origin[3] = { this->m_Volume->GetOrigin()[0],
                             this->m_Volume->GetOrigin()[1],
                             this->m_Volume->GetOrigin()[2] };
  std::cout << "Slice 0: " << image0Origin[0] << " " << image0Origin[1] << " " << image0Origin[2] << std::endl;

  // assume volume interleaving, i.e. the second dicom file stores
  // the second slice in the same volume as the first dicom file
  double image1Origin[3];

  unsigned long nextSlice = 0;
  if (this->m_Headers.size() > 1)
  {
    // assuming multiple files is invalid for single-file volume: http://www.na-mic.org/Bug/view.php?id=4105
    nextSlice = 1;
  }

  this->m_Headers[nextSlice]->GetOrigin(image1Origin);
  std::cout << "Slice " << nextSlice << ": " << image1Origin[0] << " " << image1Origin[1] << " " << image1Origin[2]
            << std::endl;

  image1Origin[0] -= image0Origin[0];
  image1Origin[1] -= image0Origin[1];
  image1Origin[2] -= image0Origin[2];
  const RotationMatrixType & NRRDSpaceDirection = this->GetNRRDSpaceDirection();
  {
    double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2]) + image1Origin[1] * (NRRDSpaceDirection[1][2]) +
                image1Origin[2] * (NRRDSpaceDirection[2][2]);
    if (x1 < 0)
    {
      this->m_SliceOrderIS = false;
    }
  }
}

/*
 * According Dicom standard:(DICOM PS3.6 2016b - Data Dictionary)
 * (0018, 0050) indicates slice thickness.
 *  which is also consistent with Dcom2iix software.
 *  In the 348-bytes header of NifTi, there is no place to store thickness information.
 *  And 348 bytes which are already occupied compactly leave no space for private information extension.
 *  In another words, you will loss thickness information from Dicom to FSL, or from Nrrd to FSL.
 * */
double
DWIDICOMConverterBase::readThicknessFromDicom() const
{
  double thickness = 0.0;
  m_Headers[0]->GetElementDS<double>(0x0018, 0x0050, 1, &thickness);
  return thickness;
}

// This getDicomSpacing method is abstracted from itk::itkDCMTKFileReader::GetSpacing of version Feb 27, 2017
// as currently ITK::GetSpacing can not correctly handle zSpace case
int
DWIDICOMConverterBase::getDicomSpacing(double * const spacing) const
{
  double _spacing[3] = { 0.0, 0.0, 0.0 };
  //
  // There are several tags that can have spacing, and we're going
  // from most to least desirable, starting with PixelSpacing, which
  // is guaranteed to be in patient space.
  // Imager Pixel spacing is inter-pixel spacing at the sensor front plane
  // Pixel spacing

  // first, shared function groups sequence, then
  // per-frame groups sequence
  int rval = m_Headers[0]->GetElementDS<double>(0x0028, 0x0030, 2, _spacing, false);
  if (rval != EXIT_SUCCESS)
  {
    // imager pixel spacing
    rval = m_Headers[0]->GetElementDS<double>(0x0018, 0x1164, 2, &_spacing[0], false);
    if (rval != EXIT_SUCCESS)
    {
      // Nominal Scanned PixelSpacing
      rval = m_Headers[0]->GetElementDS<double>(0x0018, 0x2010, 2, &_spacing[0], false);
    }
  }

  if (rval == EXIT_SUCCESS)
  {
    // slice thickness
    spacing[0] = _spacing[1];
    spacing[1] = _spacing[0];
    /*
     * According Dicom standard (DICOM PS3.6 2016b - Data Dictionary)
     * (0028, 0030) indicates physical X,Y spacing inside a slice;
     * (0018, 0088) indicates physical Z spacing between slices;
     *  which above are also consistent with Dcom2iix software.
     *  when we can not get (0018, 0088),we will revert to previous
     *  behavior and use (0018, 0050) thickness as a proxy to spacing.
     * */
    if (m_Headers[0]->GetElementDS<double>(0x0018, 0x0088, 1, &_spacing[2], false) == EXIT_SUCCESS ||
        m_Headers[0]->GetElementDS<double>(0x0018, 0x0050, 1, &_spacing[2], false) == EXIT_SUCCESS)
    {
      spacing[2] = _spacing[2];
    }
    else
    {
      // punt, thicknes of 1
      spacing[2] = 1.0;
    }
    return rval;
  }
  // this is for multiframe images -- preferentially use the shared
  // functional group, and then the per-frame functional group
  unsigned short candidateSequences[2] = {
    0x9229, // check for Shared Functional Group Sequence first
    0x9230, // check the Per-frame Functional Groups Sequence
  };
  for (unsigned short candidateSequence : candidateSequences)
  {
    brains::DCMTKSequence spacingSequence;
    rval = m_Headers[0]->GetElementSQ(0x5200, candidateSequence, spacingSequence, false);
    if (rval == EXIT_SUCCESS)
    {
      brains::DCMTKItem item;
      rval = spacingSequence.GetElementItem(0, item, false);
      if (rval == EXIT_SUCCESS)
      {
        brains::DCMTKSequence subSequence;
        // Pixel Measures Sequence
        rval = item.GetElementSQ(0x0028, 0x9110, subSequence, false);
        if (rval == EXIT_SUCCESS)
        {
          /*
           * According Dicom standard (DICOM PS3.6 2016b - Data Dictionary)
           * (0028, 0030) indicates physical X,Y spacing inside a slice;
           * (0018, 0088) indicates physical Z spacing between slices;
           *  which above are also consistent with Dcom2iix software.
           *  when we can not get (0018, 0088),we will revert to previous
           *  behavior and use (0018, 0050) thickness as a proxy to spacing.
           * */
          if (subSequence.GetElementDS<double>(0x0028, 0x0030, 2, _spacing, false) == EXIT_SUCCESS)
          {
            spacing[0] = _spacing[1];
            spacing[1] = _spacing[0];
            if (subSequence.GetElementDS<double>(0x0018, 0x0088, 1, &_spacing[2], false) == EXIT_SUCCESS ||
                subSequence.GetElementDS<double>(0x0018, 0x0050, 1, &_spacing[2], false) == EXIT_SUCCESS)
            {
              spacing[2] = _spacing[2];
            }
            else
            {
              // punt, zSpace of 1
              spacing[2] = 1.0;
            }
            break;
          }
        }
      }
    }
  }
  return rval;
}

// ---------------------------------------------------------------------------
// Fix for issue #294: Supplement 49 standard gradient tag reading.
//
// DICOM Supplement 49 (ftp://medical.nema.org/medical/dicom/final/sup49_ft.pdf)
// defines [0018,9087] (DiffusionBValue) and [0018,9089] (DiffusionGradientOrientation)
// inside the SharedFunctionalGroupsSequence (5200,9229).  Hitachi scanners have
// always used these tags; other vendors (e.g. Toshiba, some generic implementations)
// also emit them.  This method is a shared fallback so that any
// DWIDICOMConverterBase subclass can call it when vendor-specific tags are absent.
// ---------------------------------------------------------------------------
bool
DWIDICOMConverterBase::TryExtractSupp49DWIData()
{
  // Require at least one header and at least one slice per volume.
  if (this->m_Headers.empty() || this->m_SlicesPerVolume == 0)
  {
    return false;
  }

  std::vector<double>                      bValues;
  std::vector<vnl_vector_fixed<double, 3>> diffusionVectors;

  for (unsigned int k = 0; k < this->m_NSlice; k += this->m_SlicesPerVolume)
  {
    brains::DCMTKSequence sharedFGS;
    // (5200,9229) SharedFunctionalGroupsSequence
    if (this->m_Headers[k]->GetElementSQ(0x5200, 0x9229, sharedFGS) != EXIT_SUCCESS)
    {
      return false; // tags not present; caller should try vendor-specific path
    }

    double b = 0.0;
    if (sharedFGS.GetElementFD(0x0018, 0x9087, b) != EXIT_SUCCESS)
    {
      return false;
    }
    bValues.push_back(b);

    double doubleArray[3] = { 0.0, 0.0, 0.0 };
    if (sharedFGS.GetElementFD(0x0018, 0x9089, 3, doubleArray) != EXIT_SUCCESS)
    {
      return false;
    }
    vnl_vector_fixed<double, 3> vect3d(doubleArray);
    diffusionVectors.push_back(vect3d);
  }

  // Only commit results if we successfully read all volumes.
  this->m_BValues = bValues;
  this->m_DiffusionVectors = diffusionVectors;
  return true;
}
