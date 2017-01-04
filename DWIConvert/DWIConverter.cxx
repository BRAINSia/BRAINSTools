//
// Created by Johnson, Hans J on 12/3/16.
//

#include "DWIConverter.h"
#include "itkFlipImageFilter.h"
DWIConverter::DWIConverter( const FileNamesContainer &inputFileNames, const bool FSLFileFormatHorizontalBy3Rows )
        :
        m_InputFileNames(inputFileNames),
        m_FSLFileFormatHorizontalBy3Rows(FSLFileFormatHorizontalBy3Rows),
        m_SlicesPerVolume(0),
        m_NSlice(0),
        m_NVolume(0),
        m_NRRDSpaceDefinition("left-posterior-superior")
{
  this->m_MeasurementFrame.SetIdentity();
}

DWIConverter::~DWIConverter() {}

DWIConverter::RotationMatrixType DWIConverter::GetSpacingMatrix() const
{
  RotationMatrixType SpacingMatrix;
  SpacingMatrix.Fill(0.0);
  SpacingMatrix[0][0] = this->m_Volume->GetSpacing()[0];
  SpacingMatrix[1][1] = this->m_Volume->GetSpacing()[1];
  SpacingMatrix[2][2] = this->m_Volume->GetSpacing()[2];
  return SpacingMatrix;
}

const DWIMetaDataDictionaryValidator::GradientTableType&  DWIConverter::GetDiffusionVectors() const { return this->m_DiffusionVectors; }


void DWIConverter::ConvertToSingleBValueScaledDiffusionVectors()
{
  const double maxBvalue = this->GetMaxBValue();
  {
    DWIMetaDataDictionaryValidator::GradientTableType BvalueScaledDiffusionVectors(0);
    BvalueScaledDiffusionVectors.reserve(m_DiffusionVectors.size());
    for( unsigned int k = 0; k < m_DiffusionVectors.size(); ++k )
    {
      vnl_vector_fixed<double,3> vec(3);
      float scaleFactor = 0;
      if( maxBvalue > 0 )
      {
        scaleFactor = sqrt( this->m_BValues[k] / maxBvalue );
      }
      std::cout << "Scale Factor for Multiple BValues: " << k
                << " -- sqrt( " << this->m_BValues[k] << " / " << maxBvalue << " ) = "
                << scaleFactor << std::endl;
      for( unsigned ind = 0; ind < 3; ++ind )
      {
        vec[ind] = this->m_DiffusionVectors[k][ind] * scaleFactor;
      }
      BvalueScaledDiffusionVectors.push_back(vec);
      this->m_BValues[k] = maxBvalue;
    }
    this->m_DiffusionVectors = BvalueScaledDiffusionVectors;
  }
}

void DWIConverter::ConvertToMutipleBValuesUnitScaledBVectors()
{
  const double maxBvalue = this->GetMaxBValue();
  {
    for( unsigned int k = 0; k < m_DiffusionVectors.size(); ++k )
    {
      double mag = m_DiffusionVectors[k].magnitude();
      if( std::abs( mag*mag - 1.0 ) < 0.01 ) //if less than 1% differnece
      {
        mag = 1.0;  //If less than 1% difference, then assume 100%
        //This is to avoid numerical instability with computing magnitudes of gradients
      }
      m_DiffusionVectors[k].normalize();
      this->m_BValues[k] = itk::Math::Round<double>( maxBvalue*mag*mag );
    }
  }
}

Volume4DType::Pointer DWIConverter::OrientForFSLConventions( const bool toFSL)
{
  static const double FSLDesiredDirectionFlipsWRTLPS[4] = {1,-1,1,1};
  static const double DicomDesiredDirectionFlipsWRTLPS[4] = {1,1,1,1};
  this->ConvertBVectorsToIdentityMeasurementFrame();
  this->ConvertToMutipleBValuesUnitScaledBVectors();


  Volume4DType::Pointer image4D = ThreeDToFourDImage(this->GetDiffusionVolume());
  Volume4DType::DirectionType direction=image4D->GetDirection();
  direction.GetVnlMatrix().get_row(0).magnitude();
  //LPS to RAI as FSL desires images to be formatted for viewing purposes.
  // This conversion makes FSLView display the images in
  // a way that is most easily interpretable.
  typedef itk::FlipImageFilter<Volume4DType> FlipperType;
  FlipperType::Pointer myFlipper = FlipperType::New();
  myFlipper->SetInput( image4D ) ;
  FlipperType::FlipAxesArrayType arrayAxisFlip;
  for(size_t i=0; i< Volume4DType::ImageDimension; ++i)
  {
    if( toFSL )
    {
    arrayAxisFlip[i] = ( FSLDesiredDirectionFlipsWRTLPS[i]*direction(i,i) < -0.5 ); // i.e. a negative magnitude greater than 0.5
    }
    else
    {
    arrayAxisFlip[i] = ( DicomDesiredDirectionFlipsWRTLPS[i]*direction(i,i) < -0.5 ); // i.e. a negative magnitude greater than 0.5
    }
    //This is necesssary to ensure that the BVEC file is consistent with FSL orientation assumptions
    for(size_t g =0 ; g < this->m_DiffusionVectors.size(); ++g)
    {
      this->m_DiffusionVectors[g][i]  *= ( arrayAxisFlip[i] ? -1 : 1 );
    }
  }
  /* Debugging information for identifying orientation!
  std::cout << "arrayAxisFlip" << std::endl;
  for(int i =0; i < 11; ++i)
  {
    std::cout << arrayAxisFlip << std::endl;
    std::cout << m_IsInterleaved << std::endl;
    std::cout << m_SliceOrderIS << std::endl;
  }
   */
  //
  // FSL wants the second and third dimensions flipped with regards to LPS orientation
  // FSL wants the second and third dimeinsions flipped with regards to LPS orientation
  myFlipper->SetFlipAxes(arrayAxisFlip);
  myFlipper->FlipAboutOriginOff();  //Flip the image and direction cosignes
  // this is similar to a transform of [1 0 0; 0 -1 0; 0 0 -1]
  myFlipper->Update();
  Volume4DType::Pointer temp = myFlipper->GetOutput();
  temp->SetMetaDataDictionary( image4D->GetMetaDataDictionary());
  this->m_Volume = FourDToThreeDImage(temp);
  return temp;
}

const std::vector<double>& DWIConverter::GetBValues() const { return this->m_BValues; }
void  DWIConverter::SetBValues( const std::vector<double> & inBValues ) { this->m_BValues = inBValues; }
double DWIConverter::GetMaxBValue() const { return ComputeMaxBvalue( this->m_BValues); }

DWIConverter::Volume3DUnwrappedType::Pointer DWIConverter::GetDiffusionVolume() const { return this->m_Volume; }

DWIConverter::SpacingType DWIConverter::GetSpacing() const
{
  return this->m_Volume->GetSpacing();
}

DWIConverter::Volume3DUnwrappedType::PointType DWIConverter::GetOrigin() const
{
  return this->m_Volume->GetOrigin();
}

void DWIConverter::SetOrigin(DWIConverter::Volume3DUnwrappedType::PointType origin)
{
  return this->m_Volume->SetOrigin(origin);
}


DWIConverter::RotationMatrixType   DWIConverter::GetLPSDirCos() const { return this->m_Volume->GetDirection(); }

DWIConverter::RotationMatrixType DWIConverter::GetMeasurementFrame() const { return this->m_MeasurementFrame; }

DWIConverter::RotationMatrixType DWIConverter::GetNRRDSpaceDirection() const { return  this->m_Volume->GetDirection() * this->GetSpacingMatrix(); }

unsigned int DWIConverter::GetSlicesPerVolume() const { return m_SlicesPerVolume; }
unsigned int DWIConverter::GetNVolume() const { return this->m_NVolume; }
std::string DWIConverter::GetNRRDSpaceDefinition() const { return this->m_NRRDSpaceDefinition; }

unsigned short DWIConverter::GetRows() const { return this->m_Volume->GetLargestPossibleRegion().GetSize()[0]; }

unsigned short DWIConverter::GetCols() const { return this->m_Volume->GetLargestPossibleRegion().GetSize()[1]; }

void DWIConverter::ReadGradientInformation(const std::string& inputBValues, const std::string &inputBVectors, const std::string &inputVolumeNameTemplate)
{// override gradients embedded in file with an external FSL Formatted files
  std::string _inputBValues = inputBValues;
  std::string baseDirectory = itksys::SystemTools::GetParentDirectory(inputVolumeNameTemplate);
  if( CheckArg<std::string>("B Values", inputBValues, "") == EXIT_FAILURE )
  {
    std::vector<std::string> pathElements;
    pathElements.push_back(baseDirectory);
    pathElements.push_back("/");
    pathElements.push_back( itksys::SystemTools::GetFilenameWithoutExtension (inputVolumeNameTemplate) + ".bval");
    _inputBValues = itksys::SystemTools::JoinPath(pathElements);
    std::cout << "   From template " << inputVolumeNameTemplate << std::endl;
    std::cout << "   defaulting to: " << _inputBValues << std::endl;
  }
  std::string _inputBVectors = inputBVectors;
  if( CheckArg<std::string>("B Vectors", inputBVectors, "") == EXIT_FAILURE )
  {
    std::vector<std::string> pathElements;
    pathElements.push_back(baseDirectory);
    pathElements.push_back("/");
    pathElements.push_back( itksys::SystemTools::GetFilenameWithoutExtension(inputVolumeNameTemplate) + ".bvec" );
    _inputBVectors = itksys::SystemTools::JoinPath(pathElements);
    std::cout << "   From template " << inputVolumeNameTemplate << std::endl;
    std::cout << "   defaulting to: " << _inputBVectors << std::endl;
  }

  unsigned int                      bValCount = 0;
  std::vector<double>               BVals;

  if( ReadBVals(BVals, bValCount, _inputBValues) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "ERROR reading Bvals " << _inputBValues);
  }
  DWIMetaDataDictionaryValidator::GradientTableType BVecs;
  unsigned int                      bVecCount = 0;
  if( ReadBVecs(BVecs, bVecCount, _inputBVectors,true) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "ERROR reading Bvals " << _inputBVectors);
  }
  if( bValCount != bVecCount )
  {
    itkGenericExceptionMacro( << "Mismatch between count of B Vectors ("
                                      << bVecCount << ") and B Values ("
                                      << bValCount << ")" << std::endl);
  }

  size_t numGradients = BVals.size();
  if( numGradients != this->GetNVolume() )
  {
    itkGenericExceptionMacro( << "number of Gradients doesn't match number of volumes:"
                                      << numGradients << " != " << this->GetNVolume() << std::endl);
  }
  //We need to zero out BVecs for Zero BVals
  for(size_t i = 0; i < bValCount; ++i)
  {
    if ( BVals[i] < 1 )
    {
      BVecs[i][0] = 0;
      BVecs[i][1] = 0;
      BVecs[i][2] = 0;
    }
    else
    {
      BVecs[i].normalize();  //Ensure that they are unit vectors as required FSL.
    }
  }
  this->m_DiffusionVectors = BVecs;
  this->m_BValues = BVals;
  return;
}

void DWIConverter::ConvertBVectorsToIdentityMeasurementFrame()
{
  DWIMetaDataDictionaryValidator::GradientTableType gradientVectors;
  const vnl_matrix_fixed<double, 3, 3> InverseMeasurementFrame = this->GetMeasurementFrame().GetInverse();
  // grab the diffusion vectors.
  for (unsigned int k = 0; k<this->m_DiffusionVectors.size(); ++k) {
    DWIMetaDataDictionaryValidator::GradientDirectionType vec;
    // For scanners, the measurement frame for the gradient directions is the same as the
    // Excerpt from http://teem.sourceforge.net/nrrd/format.html definition of "measurement frame:"
    // There is also the possibility that a measurement frame
    // should be recorded for an image even though it is storing
    // only scalar values (e.g., a sequence of diffusion-weighted MR
    // images has a measurement frame for the coefficients of
    // the diffusion-sensitizing gradient directions, and
    // the measurement frame field is the logical store
    // this information).
    // It was noticed on oblique Philips DTI scans that the prescribed protocol directions were
    // rotated by the ImageOrientationPatient amount and recorded in the DICOM header.
    // In order to compare two different scans to determine if the same protocol was prosribed,
    // it is necessary to multiply each of the recorded diffusion gradient directions by
    // the inverse of the LPSDirCos.
    vnl_vector_fixed<double, 3> RotatedScaledDiffusionVectors =
            InverseMeasurementFrame*(this->m_DiffusionVectors[k]);
    for (unsigned ind = 0; ind<3; ++ind) {
      vec[ind] = RotatedScaledDiffusionVectors[ind];
    }
    gradientVectors.push_back(vec);
  }
  this->m_DiffusionVectors = gradientVectors;
  this->m_MeasurementFrame.SetIdentity();
}

std::string  DWIConverter::MakeFileComment(
        const std::string& version,
        bool useBMatrixGradientDirections,
        bool useIdentityMeaseurementFrame,
        double smallGradientThreshold,
        const std::string inputFileType) const
{
  std::stringstream commentSection;
  {
    commentSection << "#" << std::endl << "#" << std::endl;
    commentSection << "# This file was created by DWIConvert version " << version << std::endl
                   << "# https://github.com/BRAINSia/BRAINSTools" << std::endl
                   << "# part of the BRAINSTools package." << std::endl
                   << "# Command line options:" << std::endl
                   << "# --inputFileType " << inputFileType << std::endl;
    if( std::abs( smallGradientThreshold- 0.2 ) > 1e-4 )
    {
      commentSection << "# --smallGradientThreshold " << smallGradientThreshold << std::endl;
    }
    if (useIdentityMeaseurementFrame) {
      commentSection << "# --useIdentityMeasurementFrame" << std::endl;
    }
    if (useBMatrixGradientDirections) {
      commentSection << "# --useBMatrixGradientDirections" << std::endl;
    }
  }
  return commentSection.str();
}

void DWIConverter::ManualWriteNRRDFile(
        const std::string& outputVolumeHeaderName,
        const std::string commentstring) const
{
  const size_t extensionPos = outputVolumeHeaderName.find(".nhdr");
  const bool nrrdSingleFileFormat = ( extensionPos != std::string::npos ) ? false : true;

  std::string outputVolumeDataName;
  if( extensionPos != std::string::npos )
  {
    outputVolumeDataName = outputVolumeHeaderName.substr(0, extensionPos);
    outputVolumeDataName += ".raw";
  }

  itk::NumberToString<double> DoubleConvert;
  std::ofstream header;
  // std::string headerFileName = outputDir + "/" + outputFileName;

  const double maxBvalue = this->GetMaxBValue();
  header.open(outputVolumeHeaderName.c_str(), std::ios_base::out | std::ios_base::binary);
  header << "NRRD0005" << std::endl
         << std::setprecision(17) << std::scientific;

  header << commentstring;

  // stamp with DWIConvert branding

  if (!nrrdSingleFileFormat) {
    header << "content: exists(" << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << ",0)"
           << std::endl;
  }
  header << "type: short" << std::endl;
  header << "dimension: 4" << std::endl;
  header << "space: " << this->GetNRRDSpaceDefinition() << "" << std::endl;

  const DWIConverter::RotationMatrixType& NRRDSpaceDirection = this->GetNRRDSpaceDirection();
  header << "sizes: " << this->GetCols()
         << " " << this->GetRows()
         << " " << this->GetSlicesPerVolume()
         << " " << this->GetNVolume() << std::endl;
  header << "thicknesses:  NaN  NaN " << DoubleConvert(this->GetSpacing()[2]) << " NaN" << std::endl;
  // need to check
  header << "space directions: "
         << "("
         << DoubleConvert(NRRDSpaceDirection[0][0]) << ","
         << DoubleConvert(NRRDSpaceDirection[1][0]) << ","
         << DoubleConvert(NRRDSpaceDirection[2][0])
         << ") "
         << "("
         << DoubleConvert(NRRDSpaceDirection[0][1]) << ","
         << DoubleConvert(NRRDSpaceDirection[1][1]) << ","
         << DoubleConvert(NRRDSpaceDirection[2][1]) << ") "
         << "("
         << DoubleConvert(NRRDSpaceDirection[0][2]) << ","
         << DoubleConvert(NRRDSpaceDirection[1][2]) << ","
         << DoubleConvert(NRRDSpaceDirection[2][2])
         << ") none" << std::endl;
  header << "centerings: cell cell cell ???" << std::endl;
  header << "kinds: space space space list" << std::endl;

  header << "endian: little" << std::endl;
  header << "encoding: raw" << std::endl;
  header << "space units: \"mm\" \"mm\" \"mm\"" << std::endl;

  const DWIConverter::Volume3DUnwrappedType::PointType ImageOrigin = this->GetOrigin();
  header << "space origin: "
         << "(" << DoubleConvert(ImageOrigin[0])
         << "," << DoubleConvert(ImageOrigin[1])
         << "," << DoubleConvert(ImageOrigin[2]) << ") " << std::endl;
  if (!nrrdSingleFileFormat) {
    header << "data file: " << itksys::SystemTools::GetFilenameName(outputVolumeDataName) << std::endl;
  }

  {
    DWIConverter::RotationMatrixType MeasurementFrame = this->GetMeasurementFrame();
    header << "measurement frame: "
           << "(" << DoubleConvert(MeasurementFrame[0][0]) << ","
           << DoubleConvert(MeasurementFrame[1][0]) << ","
           << DoubleConvert(MeasurementFrame[2][0]) << ") "
           << "(" << DoubleConvert(MeasurementFrame[0][1]) << ","
           << DoubleConvert(MeasurementFrame[1][1]) << ","
           << DoubleConvert(MeasurementFrame[2][1]) << ") "
           << "(" << DoubleConvert(MeasurementFrame[0][2]) << ","
           << DoubleConvert(MeasurementFrame[1][2]) << ","
           << DoubleConvert(MeasurementFrame[2][2]) << ")"
           << std::endl;
  }

  for(std::map<std::string,std::string>::const_iterator it=this->m_CommonDicomFieldsMap.begin();
      it != this->m_CommonDicomFieldsMap.end(); ++it)
  {
    header << it->first << ":=" << it->second << std::endl;
  }

  header << "modality:=DWMRI" << std::endl;
  // this is the norminal BValue, i.e. the largest one.
  header << "DWMRI_b-value:=" << DoubleConvert(maxBvalue) << std::endl;

  //  the following three lines are for older NRRD format, where
  //  baseline images are always in the begining.
  //  header << "DWMRI_gradient_0000:=0  0  0" << std::endl;
  //  header << "DWMRI_NEX_0000:=" << nBaseline << std::endl;
  //  need to check
  {
    const DWIMetaDataDictionaryValidator::GradientTableType & gradientVectors = this->m_DiffusionVectors;
    unsigned int gradientVecIndex = 0;
    for (unsigned int k = 0; k<gradientVectors.size(); ++k) {
      header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << ":="
             << DoubleConvert(gradientVectors[gradientVecIndex][0]) << "   "
             << DoubleConvert(gradientVectors[gradientVecIndex][1]) << "   "
             << DoubleConvert(gradientVectors[gradientVecIndex][2])
             << std::endl;
      ++gradientVecIndex;
    }
  }
  // write data in the same file is .nrrd was chosen
  header << std::endl;;
  if (nrrdSingleFileFormat) {
    unsigned long nVoxels = this->GetDiffusionVolume()->GetBufferedRegion().GetNumberOfPixels();
    header.write(reinterpret_cast<char*>(this->GetDiffusionVolume()->GetBufferPointer()),
                 nVoxels*sizeof(short));
  }
  else {
    // if we're writing out NRRD, and the split header/data NRRD
    // format is used, write out the image as a raw volume.
    itk::ImageFileWriter<DWIConverter::Volume3DUnwrappedType>::Pointer
            rawWriter = itk::ImageFileWriter<DWIConverter::Volume3DUnwrappedType>::New();
    itk::RawImageIO<PixelValueType,3>::Pointer rawIO
            = itk::RawImageIO<PixelValueType,3>::New();
    rawWriter->SetImageIO(rawIO);
    rawIO->SetByteOrderToLittleEndian();
    rawWriter->SetFileName(outputVolumeDataName.c_str());
    rawWriter->SetInput(this->GetDiffusionVolume());
    try {
      rawWriter->Update();
    }
    catch (itk::ExceptionObject& excp) {
      std::cerr << "Exception thrown while writing the series to"
                << outputVolumeDataName << " " << excp << std::endl;
      std::cerr << excp << std::endl;
    }
  }
  header.close();
  return;
}

Volume4DType::Pointer DWIConverter::ThreeDToFourDImage(Volume3DUnwrappedType::Pointer img) const
{
  const int nVolumes = this->GetNVolume();

  DWIConverter::Volume3DUnwrappedType::SizeType size3D(img->GetLargestPossibleRegion().GetSize());
  DWIConverter::Volume3DUnwrappedType::DirectionType direction3D(img->GetDirection());
  DWIConverter::Volume3DUnwrappedType::SpacingType spacing3D(img->GetSpacing());
  DWIConverter::Volume3DUnwrappedType::PointType origin3D(img->GetOrigin());

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

DWIConverter::Volume3DUnwrappedType::Pointer DWIConverter::FourDToThreeDImage(Volume4DType::Pointer img4D) const
{
  Volume4DType::SizeType      size4D(img4D->GetLargestPossibleRegion().GetSize() );
  Volume4DType::DirectionType direction4D(img4D->GetDirection() );
  Volume4DType::SpacingType   spacing4D(img4D->GetSpacing() );
  Volume4DType::PointType     origin4D(img4D->GetOrigin() );

  Volume3DUnwrappedType::RegionType region3D;
  {
    Volume3DUnwrappedType::SizeType size3D;
    size3D[0] = size4D[0];
    size3D[1] = size4D[1];
    const int nVolumes = img4D->GetLargestPossibleRegion().GetSize()[3];
    size3D[2] = size4D[2] * nVolumes;

    Volume3DUnwrappedType::IndexType index3D;
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
  Volume3DUnwrappedType::DirectionType direction3D;
  direction3D.SetIdentity();
  Volume3DUnwrappedType::SpacingType   spacing3D;
  spacing3D.Fill(1.0);
  Volume3DUnwrappedType::PointType     origin3D;
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

  Volume3DUnwrappedType::Pointer img = Volume3DUnwrappedType::New();
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

void DWIConverter::WriteFSLFormattedFileSet(const std::string& outputVolumeHeaderName,
                                            const std::string outputBValues, const std::string outputBVectors, Volume4DType::Pointer img4D) const
{
  const double trace = this->m_MeasurementFrame[0][0] * this->m_MeasurementFrame[1][1] *
                       this->m_MeasurementFrame[2][2];
  if( std::abs( trace - 1.0 ) > 1e-4 )
  {
    itkGenericExceptionMacro( << "ERROR:  Only identity measurement frame allow for writing FSL formatted files "
                                      << std::endl);
  }

  {
    //Set the qform and sfrom codes for the MetaDataDictionary.
    itk::MetaDataDictionary & thisDic = img4D->GetMetaDataDictionary();
    itk::EncapsulateMetaData< std::string >( thisDic, "qform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
    itk::EncapsulateMetaData< std::string >( thisDic, "sform_code_name", "NIFTI_XFORM_SCANNER_ANAT" );
  }
  itk::ImageFileWriter<Volume4DType>::Pointer imgWriter = itk::ImageFileWriter<Volume4DType>::New();
  imgWriter->SetInput( img4D );
  imgWriter->SetFileName( outputVolumeHeaderName.c_str() );
  try
  {
    imgWriter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while writing "
              << outputVolumeHeaderName << std::endl;
    std::cerr << excp << std::endl;
    throw;
  }
  // FSL output of gradients & BValues
  std::string outputFSLBValFilename;
  //
  const size_t extensionPos = this->has_valid_nifti_extension(outputVolumeHeaderName);
  if( outputBValues == "" )
  {
    outputFSLBValFilename = outputVolumeHeaderName.substr(0, extensionPos);
    outputFSLBValFilename += ".bval";
  }
  else
  {
    outputFSLBValFilename = outputBValues;
  }
  std::string outputFSLBVecFilename;
  if( outputBVectors == "" )
  {
    outputFSLBVecFilename = outputVolumeHeaderName.substr(0, extensionPos);
    outputFSLBVecFilename += ".bvec";
  }
  else
  {
    outputFSLBVecFilename = outputBVectors;
  }

  // write out in FSL format
  if( WriteBValues<double>(this->m_BValues, outputFSLBValFilename) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "Failed to write FSL BVal File: " << outputFSLBValFilename << std::endl;);
  }
  if( WriteBVectors( this->m_DiffusionVectors , outputFSLBVecFilename) != EXIT_SUCCESS )
  {
    itkGenericExceptionMacro(<< "Failed to write FSL BVec File: " << outputFSLBVecFilename << std::endl;);
  }
}

void DWIConverter::SetAllowLossyConversion(const bool newValue) { this->m_allowLossyConversion = newValue; }

double DWIConverter::ComputeMaxBvalue(const std::vector<double> &bValues) const
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

size_t DWIConverter::has_valid_nifti_extension( std::string outputVolumeHeaderName ) const
{
  const size_t NUMEXT=2;
  const char * const extList [NUMEXT] = {".nii.gz", ".nii"};
  for(size_t i = 0 ; i < NUMEXT; ++i)
  {
    const size_t extensionPos = outputVolumeHeaderName.find(extList[i]);
    if( extensionPos != std::string::npos )
    {
      return extensionPos;
    }
  }
  {
    std::cerr << "FSL Format output chosen, "
              << "but output Volume not a recognized "
              << "NIfTI filename " << outputVolumeHeaderName
              << std::endl;
    exit(1);
  }
  return std::string::npos;
}

DWIConverter::Volume3DUnwrappedType::Pointer DWIConverter::getVolumePointer(){
  return m_Volume;
}