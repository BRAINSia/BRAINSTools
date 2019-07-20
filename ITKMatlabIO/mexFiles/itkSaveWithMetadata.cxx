#include "itkImage.h"
#include "itksys/SystemTools.hxx"
#include "itkMetaDataObject.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "nrrdCommon.h"
#include "itkNumberToString.h"
#include "myMexPrintf.h" //add by HuiXie Nov 15th, 2016


// add these 2 function to replace the use of vnl_vector_fixed, by HuiXie Nov 15th, 2016
double
getVecMagnitude( const std::vector< double > vec )
{
  double    sum = 0.0;
  int const size = vec.size();
  for ( int i = 0; i < size; ++i )
  {
    sum += vec[i] * vec[i];
  }
  return sqrt( sum );
}

void
normalizeVec( std::vector< double > & vec )
{
  double    magnitude = getVecMagnitude( vec );
  int const size = vec.size();
  for ( int i = 0; i < size; ++i )
  {
    vec[i] = vec[i] / magnitude;
    if ( vec[i] < 1e-10 )
      vec[i] = 0;
  }
}

void
printVecValue( const std::vector< double > vec )
{
  int const size = vec.size();
  for ( int i = 0; i < size; ++i )
  {
    std::cerr << "i= " << i << "  vec[i]= " << vec[i] << std::endl;
  }
  std::cout << std::endl;
}

inline itk::ImageIOBase::IOComponentType
typeMtoITK( const mxClassID mtype )
{

  switch ( mtype )
  {
    case mxINT8_CLASS:
      return itk::ImageIOBase::CHAR;
    case mxUINT8_CLASS:
      return itk::ImageIOBase::UCHAR;
    case mxINT16_CLASS:
      return itk::ImageIOBase::SHORT;
    case mxUINT16_CLASS:
      return itk::ImageIOBase::USHORT;
    case mxINT32_CLASS:
      return itk::ImageIOBase::INT;
    case mxUINT32_CLASS:
      return itk::ImageIOBase::UINT;
    case mxINT64_CLASS:
      return itk::ImageIOBase::LONG;
    case mxUINT64_CLASS:
      return itk::ImageIOBase::ULONG;
    case mxSINGLE_CLASS:
      return itk::ImageIOBase::FLOAT;
    case mxDOUBLE_CLASS:
      return itk::ImageIOBase::DOUBLE;
    default:
      break;
  } // end switch
  return itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;
} // end typeMtoITK

/** need a specialization that only calls
 * SetNumberOfComponentsPerPixel for VectorImages
 */
using DWIImage = itk::VectorImage< double, 3 >;

// TODO: Remove this code

template < typename TImage >
void
SetNumberOfComponentsPerPixel( typename TImage::Pointer im, unsigned numComponents )
{
  im->SetNumberOfComponentsPerPixel( numComponents );
}

#if 0
template<>
void SetNumberOfComponentsPerPixel<DWIImage>(DWIImage::Pointer im, unsigned numComponents)
{
  im->SetNumberOfComponentsPerPixel(numComponents);
}
#endif

template < typename TImage >
void
WriteITKImageFromMatlabStructure( const MatlabStructManager & msm, const char * filename )
{
  using ImageType = TImage;
  using WriterType = itk::ImageFileWriter< ImageType >;

  const unsigned int           numDims = msm.GetNumberOfDimensions( "data" );
  const mwSize * const         mSize = msm.GetDimensions( "data" );
  typename ImageType::SizeType itkSize;
  size_t                       numPixels = 1;
  {
    unsigned int axIdx;
    for ( axIdx = 0; axIdx < ImageType::ImageDimension; axIdx++ )
    {
      itkSize[axIdx] = mSize[axIdx];
      numPixels *= itkSize[axIdx];
    }

    //???? this code has problem, it should never be executed. by Hui Xie
    if ( axIdx < numDims )
    {
      numPixels *= mSize[axIdx];
    }
  }

  /** spaceorigin **/
  typename ImageType::PointType itkOrigin;
  itkOrigin.Fill( 0.0 );
  const double * spaceorigin_temp = (double *)mxGetData( msm.GetField( "spaceorigin" ) );

  const unsigned int spatialDims = numDims; // instead of 3
  for ( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++ )
  {
    if ( sdIdx < spatialDims )
    {
      itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
    }
    else
    {
      itkOrigin[sdIdx] = 0.0;
    }
  }

  typename ImageType::SpacingType   itkSpacing;
  typename ImageType::DirectionType itkDirection;

  // fill out extra dimension for (usually) 4 D -- add zero to 4th
  // entry in direction vector, except for the last vector, which
  // needs to be 0 0 0 1
  itkDirection.SetIdentity();
  /** spacedirections **/
  const double * spacedirections_temp = (double *)mxGetData( msm.GetField( "spacedirections" ) );
  for ( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx )
  {
    // vnl_vector_fixed<double, spatialDims> vec;
    std::vector< double > vec;
    vec.reserve( spatialDims );
    for ( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
    {
      const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
      // vec[sdIdx] = spacedirections_temp[sdir_offset];
      vec.push_back( spacedirections_temp[sdir_offset] );
    }
    // itkSpacing[axIdx] = vec.magnitude();
    itkSpacing[axIdx] = getVecMagnitude( vec );
    // vec.normalize();
    normalizeVec( vec );
    for ( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
    {
      itkDirection[sdIdx][axIdx] = vec[sdIdx];
    }
  }

  double      flipFactors[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int * nrrdSpaceDefinition = (int *)mxGetData( msm.GetField( "space" ) );
  if ( nrrdSpaceDefinition )
  {
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch ( *nrrdSpaceDefinition )
    {
      case nrrdSpaceLeftPosteriorSuperior:
        // already in ITK space format
        break;
      case nrrdSpaceRightAnteriorSuperior:
        for ( unsigned int i = 0; i < ImageType::ImageDimension; ++i )
        {
          itkDirection[i][0] *= -1;
          itkDirection[i][1] *= -1;
        }
        flipFactors[0] = -1.0;
        flipFactors[1] = -1.0;
        break;
      case nrrdSpaceLeftAnteriorSuperior:
        for ( unsigned int i = 0; i < ImageType::ImageDimension; ++i )
        {
          itkDirection[i][0] *= -1;
        }
        flipFactors[0] = -1.0;
        break;
      default:
        myMexPrintf( "nrrdSpaceDefinition failed into any preset case in WriteITKImageFromMatlabStructure function." );
        break;
    }
  }

  typename ImageType::Pointer im = ImageType::New();
  im->SetDirection( itkDirection );
  im->SetOrigin( itkOrigin );
  im->SetSpacing( itkSpacing );
  im->SetRegions( itkSize );
  size_t ComponentsPerPixel = 1;

  // if we have a vector image
  if ( numDims > ImageType::ImageDimension )
  {
    ComponentsPerPixel = mSize[ImageType::ImageDimension];
    SetNumberOfComponentsPerPixel< ImageType >( im, ComponentsPerPixel );
  }
  im->Allocate();
  const mxArray * const dataMx = msm.GetField( "data" );
  // Note: Matlab returns the internal value types, not a vector of those types.
  const typename TImage::InternalPixelType * voxels =
    static_cast< typename TImage::InternalPixelType * >( mxGetData( dataMx ) );
  {
    std::copy( voxels, voxels + numPixels * ComponentsPerPixel, im->GetBufferPointer() );
  }

  itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();

  // space direction
  std::string spaceDirKey( "NRRD_" );
  spaceDirKey += airEnumStr( nrrdSpace, nrrdField_space );
  itk::EncapsulateMetaData< std::string >(
    thisDic, spaceDirKey, std::string( airEnumStr( nrrdSpace, nrrdSpaceLeftPosteriorSuperior ) ) );

  // fill out metadata
  // Measurement Frame
  mxArray * mxMeasurementFrame = msm.GetField( "measurementframe" );
  /** measurementframe **/
  if ( msm.isDWIdata() && mxMeasurementFrame )
  {
    const double * const mxMeasurementFrame_temp = (double *)mxGetData( mxMeasurementFrame );
    const mwSize * const measurementFrameSize = msm.GetDimensions( "data" );
    if ( mxMeasurementFrame != 0 )
    {
      using MeasurementMatType = std::vector< std::vector< double > >;
      MeasurementMatType measurementFrame;
      std::string        measurementFrameFieldName = "NRRD_";
      measurementFrameFieldName += airEnumStr( nrrdField, nrrdField_measurement_frame );
      for ( unsigned int i = 0, count = 0; i < ImageType::ImageDimension; ++i )
      {
        // if the saved measurement frame is N-d and weire at N+1-D...
        std::vector< double > tmp;
        if ( i >= measurementFrameSize[0] )
        {
          for ( unsigned int j = 0; j < ImageType::ImageDimension; ++j )
          {
            if ( j == i )
            {
              tmp.push_back( 1.0 );
            }
            else
            {
              tmp.push_back( 0.0 );
            } // end else
          }   // end for
        }
        else
        {
          for ( unsigned int j = 0; j < ImageType::ImageDimension; ++j )
          {
            if ( j < measurementFrameSize[0] )
            {
              tmp.push_back( mxMeasurementFrame_temp[count] * flipFactors[j] );
              ++count;
            }
            else if ( j == i )
            {
              tmp.push_back( 1.0 );
            }
            else
            {
              tmp.push_back( 0.0 );
            }
          }
        }
        measurementFrame.push_back( tmp );
      }
      itk::EncapsulateMetaData< MeasurementMatType >( thisDic, measurementFrameFieldName, measurementFrame );
    }
  }

  // centerings.
  const mxArray * const centerings_temp_MxArray = msm.GetField( "centerings" );
  if ( centerings_temp_MxArray )
  {
    const int * centerings_temp = (int *)mxGetData( centerings_temp_MxArray );
    for ( unsigned int i = 0; i < numDims; ++i )
    {
      std::stringstream attName;
      attName << "NRRD_centerings[" << i << "]";
      std::string val = airEnumStr( nrrdCenter, centerings_temp[i] );
      itk::EncapsulateMetaData< std::string >( thisDic, attName.str(), val );
    }
  }

  // kinds
  const mxArray * const kinds_temp_MxArray = msm.GetField( "kinds" );
  const int * const     kinds_temp = (int *)mxGetData( kinds_temp_MxArray );
  if ( kinds_temp_MxArray )
  {
    for ( unsigned int axIdx = 0; axIdx < numDims; ++axIdx )
    {
      std::stringstream attName;
      attName << "NRRD_kinds[" << axIdx << "]";
      std::string val = airEnumStr( nrrdKind, kinds_temp[axIdx] );
      itk::EncapsulateMetaData< std::string >( thisDic, attName.str(), val );
    }
  }

  // space units is not necessary
  //  std::string spaceUnitString;
  //  for (unsigned int i = 0; i < numDims; ++i)
  //  {
  //    if (2 == kinds_temp[i]) //space
  //      spaceUnitString.append("\"mm\" ");
  //  }
  //
  //  itk::EncapsulateMetaData<std::string>(thisDic, std::string("space units"), spaceUnitString);

  if ( msm.isDWIdata() )
  {
    // One of the directions is the gradient list instead of a spaceDim
    { // Get the DWI specific fields.
      const mxArray * const gradientdirections = msm.GetField( "gradientdirections" );
      if ( gradientdirections != 0 )
      {
        std::stringstream ss;
        ss.precision( 10 );
        const unsigned numGradients = mxGetNumberOfElements( gradientdirections ) / 3;
        if ( numGradients > 0 )
        {
          const double * const gradients = static_cast< double * >( mxGetData( gradientdirections ) );
          for ( unsigned i = 0; i < numGradients; ++i )
          {
            ss.str( std::string() );
            ss << gradients[i] << " " << gradients[i + numGradients] << " " << gradients[i + ( 2 * numGradients )];
            std::stringstream gradName;
            gradName << GRADIENT_PREFIX << std::setw( 4 ) << std::setfill( '0' ) << i;
            itk::EncapsulateMetaData< std::string >( thisDic, gradName.str(), ss.str() );
          }
          itk::EncapsulateMetaData< std::string >( thisDic, std::string( "modality" ), std::string( "DWMRI" ) );
          const double * const bvalue = static_cast< double * >( mxGetData( msm.GetField( "bvalue" ) ) );
          ss << *bvalue;
          itk::EncapsulateMetaData< std::string >( thisDic, std::string( "DWMRI_b-value" ), ss.str() );
        }
      }
      // OK
    }
  }

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename );
  writer->SetInput( im );
  try
  {
    writer->Write();
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::string msg( "itkSaveWithMetaData: can't read\n" );
    msg += filename;
    msg += " ";
    msg += excp.what();
    char errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
    snprintf( errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: at line %d", msg.c_str(), __LINE__ );
    mexErrMsgTxt( errBuff );
  }
} // end WriteITKImageFromMatlabStructure

// comment by Hui Xie on Nov 15th, 2016 as unused code
/*
template<typename TImage>
void
WriteDWIFile(const MatlabStructManager &msm, const char *filename) {
    using ImageType = TImage;
    using WriterType = itk::ImageFileWriter<ImageType>;

    const unsigned int numDims = msm.GetNumberOfDimensions("data");
    const mwSize *const mSize = msm.GetDimensions("data");
    myMexPrintf(" size =\n");

    typename ImageType::SizeType itkSize;
    size_t numPixels = 1;
    {
        unsigned int axIdx;
        for (axIdx = 0; axIdx < ImageType::ImageDimension; axIdx++) {
            itkSize[axIdx] = mSize[axIdx];
            numPixels *= itkSize[axIdx];
            myMexPrintf("%lf, %d, %d\n", itkSize[axIdx], axIdx, ImageType::ImageDimension);
        }
        if (axIdx < numDims) {
            numPixels *= mSize[axIdx];
        }
    }

    */
/** spaceorigin **/      /*
    typename ImageType::PointType itkOrigin;
    itkOrigin.Fill(0.0);
    const double *spaceorigin_temp = (double *) mxGetData(msm.GetField("spaceorigin"));
    // TODO:  Make work for 2D, but currently only works for 3D and 3D vectors.
    constexpr unsigned int spatialDims = 3; // msm.GetNumberOfDimensions("data");
    if (spatialDims != 3) {
        myMexPrintf("ERROR: ONLY 3D images (3D+gradients OK) supported.\n");
        return;
    }
    myMexPrintf("Space Origin =\n");
    for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++) {
        if (sdIdx < spatialDims) {
            itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
        } else {
            itkOrigin[sdIdx] = 0.0;
        }
        myMexPrintf("%lf, %d, %d\n", itkOrigin[sdIdx], sdIdx, ImageType::ImageDimension);
    }
     
    typename ImageType::SpacingType itkSpacing;
    typename ImageType::DirectionType itkDirection;
    //
    // fill out extra dimension for (usually) 4 D -- add zero to 4th
    // entry in direction vector, except for the last vector, which
    // needs to be 0 0 0 1
    itkDirection.SetIdentity();
    */
/** spacedirections **/  /*
const double *spacedirections_temp = (double *) mxGetData(msm.GetField("spacedirections"));
for (unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx) {
    vnl_vector_fixed<double, spatialDims> vec;
    for (unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx) {
        const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
        vec[sdIdx] = spacedirections_temp[sdir_offset];
    }
    itkSpacing[axIdx] = vec.magnitude();
    vec.normalize();
    for (unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx) {
        itkDirection[sdIdx][axIdx] = vec[sdIdx];
    }
}
myMexPrintf("directions =\n");
for (unsigned int axIdx = 0; axIdx < ImageType::ImageDimension; ++axIdx) {
    for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx) {
        myMexPrintf("%lf, %d, %d\n", itkDirection[sdIdx][axIdx], sdIdx, ImageType::ImageDimension);
    }
}
 
double flipFactors[4] = {1.0, 1.0, 1.0, 1.0};
const int *nrrdSpaceDefinition = (int *) mxGetData(msm.GetField("space"));
if (nrrdSpaceDefinition) {
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch (*nrrdSpaceDefinition) {
        case nrrdSpaceLeftPosteriorSuperior:
            // already in ITK space format
            break;
        case nrrdSpaceRightAnteriorSuperior:
            for (unsigned int i = 0; i < ImageType::ImageDimension; ++i) {
                itkDirection[i][0] *= -1;
                itkDirection[i][1] *= -1;
            }
            flipFactors[0] = -1.0;
            flipFactors[1] = -1.0;
            break;
        case nrrdSpaceLeftAnteriorSuperior:
            for (unsigned int i = 0; i < ImageType::ImageDimension; ++i) {
                itkDirection[i][0] *= -1;
            }
            flipFactors[0] = -1.0;
            break;
    }
}
 
typename ImageType::Pointer im = ImageType::New();
im->SetDirection(itkDirection);
im->SetOrigin(itkOrigin);
im->SetSpacing(itkSpacing);
im->SetRegions(itkSize);
size_t ComponentsPerPixel = 1;
//
// if we have a vector image
if (numDims > ImageType::ImageDimension) {
    ComponentsPerPixel = mSize[ImageType::ImageDimension];
    SetNumberOfComponentsPerPixel<ImageType>(im, ComponentsPerPixel);
}
im->Allocate();
const mxArray *const dataMx = msm.GetField("data");
//Note: Matlab returns the internal value types, not a vector of those types.
const typename TImage::InternalPixelType *voxels =
        static_cast<typename TImage::InternalPixelType *>( mxGetData(dataMx));
if (numDims > ImageType::ImageDimension) {
    itk::ImageRegionIterator<TImage> it(im, im->GetLargestPossibleRegion());
    while (!it.IsAtEnd()) {
        itk::VariableLengthVector<typename TImage::InternalPixelType> vlv(voxels, ComponentsPerPixel, false);
        voxels += ComponentsPerPixel;
        it.Set(vlv);
    }
}
 
itk::MetaDataDictionary &thisDic = im->GetMetaDataDictionary();
//
// space direction
std::string spaceDirKey("NRRD_space");
spaceDirKey += airEnumStr(nrrdSpace, nrrdField_space);
 
myMexPrintf("%s\n", spaceDirKey.c_str());
itk::EncapsulateMetaData<std::string>(thisDic, spaceDirKey,
                                      std::string(airEnumStr(nrrdSpace, nrrdSpaceLeftPosteriorSuperior)));
 
//
// fill out metadata
// Measurement Frame
mxArray *mxMeasurementFrame = msm.GetField("measurementframe");
*/
/** measurementframe **/ /*
if (mxMeasurementFrame) {
   const double *const mxMeasurementFrame_temp = (double *) mxGetData(mxMeasurementFrame);
   const mwSize *const measurementFrameSize = msm.GetDimensions("data");
   if (mxMeasurementFrame != 0) {
       using MeasurementMatType = std::vector<std::vector<double> >;
       MeasurementMatType measurementFrame;
       std::string measurementFrameFieldName = "NRRD_";
       measurementFrameFieldName += airEnumStr(nrrdField, nrrdField_measurement_frame);
       for (unsigned int i = 0, count = 0; i < ImageType::ImageDimension; ++i) {
           // if the saved measurement frame is N-d and weire at N+1-D...
           std::vector<double> tmp;
           if (i >= measurementFrameSize[0]) {
               for (unsigned int j = 0; j < ImageType::ImageDimension; ++j) {
                   if (j == i) {
                       tmp.push_back(1.0);
                   } else {
                       tmp.push_back(0.0);
                   }
               }
           } else {
               for (unsigned int j = 0; j < ImageType::ImageDimension; ++j) {
                   if (j < measurementFrameSize[0]) {
                       tmp.push_back(mxMeasurementFrame_temp[count] * flipFactors[j]);
                       ++count;
                   } else if (j == i) {
                       tmp.push_back(1.0);
                   } else {
                       tmp.push_back(0.0);
                   }
               }
           }
           measurementFrame.push_back(tmp);
       }
       itk::EncapsulateMetaData<MeasurementMatType>(thisDic, measurementFrameFieldName, measurementFrame);
   }
}
//
// centerings.
const mxArray *const centerings_temp_MxArray = msm.GetField("centerings");
if (centerings_temp_MxArray) {
   const int *centerings_temp = (int *) mxGetData(centerings_temp_MxArray);
   for (unsigned int i = 0; i < numDims; ++i) {
       std::stringstream attName;
       attName << "NRRD_centerings[" << i << "]";
       std::string val = airEnumStr(nrrdCenter, centerings_temp[i]);
       itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
   }
}
// kinds
const mxArray *const kinds_temp_MxArray = msm.GetField("kinds");
if (kinds_temp_MxArray) {
   const int *const kinds_temp = (int *) mxGetData(kinds_temp_MxArray);
   for (unsigned int axIdx = 0; axIdx < numDims; ++axIdx) {
       std::stringstream attName;
       attName << "NRRD_kinds[" << axIdx << "]";
       std::string val = airEnumStr(nrrdKind, kinds_temp[axIdx]);
       itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
   }
}

// space units
itk::EncapsulateMetaData<std::string>(thisDic, std::string("space units"),
                                     std::string("\"mm\" \"mm\" \"mm\""));

if (msm.isDWIdata()) {
   // One of the directions is the gradient list instead of a spaceDim
   { // Get the DWI specific fields.
       const mxArray *const gradientdirections = msm.GetField("gradientdirections");
       if (gradientdirections != 0) {
           std::stringstream ss;
           ss.precision(10);
           const unsigned numGradients = mxGetNumberOfElements(gradientdirections) / 3;
           if (numGradients > 0) {
               const double *const gradients = static_cast<double *>(mxGetData(gradientdirections));
               for (unsigned i = 0; i < numGradients; ++i) {
                   ss.str(std::string());
                   ss << gradients[i] << " "
                      << gradients[i + numGradients] << " "
                      << gradients[i + (2 * numGradients)];
                   std::stringstream gradName;
                   gradName << GRADIENT_PREFIX << std::setw(4) << std::setfill('0') << i;
                   itk::EncapsulateMetaData<std::string>(thisDic,
                                                         gradName.str(), ss.str());
               }
               itk::EncapsulateMetaData<std::string>(thisDic, std::string("modality"),
                                                     std::string("DWMRI"));
               const double *const bvalue = static_cast<double *>(mxGetData(msm.GetField("bvalue")));
               ss << *bvalue;
               itk::EncapsulateMetaData<std::string>(thisDic, std::string("DWMRI_b-value"),
                                                     ss.str());
           }
       }
       // OK
   }
}

typename WriterType::Pointer writer = WriterType::New();
writer->SetFileName(filename);
writer->SetInput(im);
try {
   writer->Write();
}
catch (itk::ExceptionObject &excp) {
   std::string msg("itkSaveWithMetaData: can't read\n");
   msg += filename;
   msg += " ";
   msg += excp.what();
   char errBuff[NRRD_MAX_ERROR_MSG_SIZE] = {'\0'};
   snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: at line %d", msg.c_str(), __LINE__);
   mexErrMsgTxt(errBuff);
}
}// end WriteDWIFile*/


// Only used in 4D image.
template < typename TScalar >
void
WriteDWINrrd( const MatlabStructManager & msm, const char * filename, const char * voxelTypeName )
{
  using ImageType = typename itk::Image< TScalar, 4 >;
  std::string _fname( filename );
  std::string ext = itksys::SystemTools::GetFilenameLastExtension( _fname );
  std::string dataFilename;
  bool        nrrdFormat( true );
  if ( ext == ".nhdr" )
  {
    nrrdFormat = false;
    dataFilename = itksys::SystemTools::GetFilenameWithoutLastExtension( _fname );
    dataFilename += ".raw"; // dump out uncompressed.
  }
  else if ( ext != ".nrrd" )
  {
    std::string msg = "Only NRRD file types recognized, extension = ";
    msg += ext;
    mexErrMsgTxt( msg.c_str() );
    return;
  }

  const unsigned int   numDims = msm.GetNumberOfDimensions( "data" );
  const mwSize * const mSize = msm.GetDimensions( "data" );

  size_t numPixels = 1;
  {
    unsigned int axIdx;
    for ( axIdx = 0; axIdx < numDims; axIdx++ )
    {
      numPixels *= mSize[axIdx];
    }
  }

  /** spaceorigin **/
  typename ImageType::PointType itkOrigin;
  itkOrigin.Fill( 0.0 );
  const double * spaceorigin_temp = (double *)mxGetData( msm.GetField( "spaceorigin" ) );
  // WriteDWINrrd only work for 4D image.
  constexpr unsigned int spatialDims = 3; // spatialDims is different with numDims
  for ( unsigned int sdIdx = 0; sdIdx < spatialDims; sdIdx++ )
  {
    if ( sdIdx < spatialDims )
    {
      itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
    }
  }

  typename ImageType::SpacingType   itkSpacing;
  typename ImageType::DirectionType itkDirection;

  // fill out extra dimension for (usually) 4 D -- add zero to 4th
  // entry in direction vector, except for the last vector, which
  // needs to be 0 0 0 1
  itkDirection.SetIdentity();
  /** spacedirections **/
  const double * spacedirections_temp = (double *)mxGetData( msm.GetField( "spacedirections" ) );
  for ( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx )
  {
    vnl_vector_fixed< double, spatialDims > vec;
    // std::vector<double> vec;
    // vec.clear();
    for ( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
    {
      const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
      vec[sdIdx] = spacedirections_temp[sdir_offset];
      // vec.push_back(spacedirections_temp[sdir_offset]);
    }
    itkSpacing[axIdx] = vec.magnitude();
    // itkSpacing[axIdx] = getVecMagnitude(vec);
    // vec.normalize();
    // normalizeVec(vec);
    for ( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
    {
      itkDirection[sdIdx][axIdx] = vec[sdIdx];
    }
  }

  double      flipFactors[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int * nrrdSpaceDefinition = (int *)mxGetData( msm.GetField( "space" ) );
  if ( nrrdSpaceDefinition )
  {
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch ( *nrrdSpaceDefinition )
    {
      case nrrdSpaceLeftPosteriorSuperior:
        // already in ITK space format
        break;
      case nrrdSpaceRightAnteriorSuperior:
        for ( unsigned int i = 0; i < spatialDims; ++i )
        {
          itkDirection[i][0] *= -1;
          itkDirection[i][1] *= -1;
        }
        flipFactors[0] = -1.0;
        flipFactors[1] = -1.0;
        break;
      case nrrdSpaceLeftAnteriorSuperior:
        for ( unsigned int i = 0; i < spatialDims; ++i )
        {
          itkDirection[i][0] *= -1;
        }
        flipFactors[0] = -1.0;
        break;
      default:
        myMexPrintf( "nrrdSpaceDefinition failed into any preset case in WriteDWINrrd function." );
        break;
    }
  }
  std::ofstream header;
  // std::string headerFileName = outputDir + "/" + outputFileName;

  header.open( filename, std::ios::out | std::ios::binary );
  header << "NRRD0005" << std::endl << std::setprecision( 17 ) << std::scientific;

  // stamp with DWIConvert branding
  header << "# This file was created by itkSaveMetaData version 1.0" << std::endl;
  header << "#" << std::endl << "#" << std::endl;

  if ( !nrrdFormat )
  {
    header << "content: exists(" << itksys::SystemTools::GetFilenameName( dataFilename ) << ",0)" << std::endl;
  }
  header << "type: " << voxelTypeName << std::endl;
  header << "dimension: 4" << std::endl;
  header << "space: ";
  if ( !nrrdSpaceDefinition || *nrrdSpaceDefinition == nrrdSpaceLeftPosteriorSuperior )
  {
    header << "left-posterior-superior" << std::endl;
  }
  else if ( *nrrdSpaceDefinition == nrrdSpaceRightAnteriorSuperior )
  {
    header << "right-anterior-superior" << std::endl;
  }
  else if ( *nrrdSpaceDefinition == nrrdSpaceLeftAnteriorSuperior )
  {
    header << "left-anterior-superior" << std::endl;
  }

  itk::NumberToString< double > DoubleConvert;

  header << "sizes: " << mSize[0] << " " << mSize[1] << " " << mSize[2] << " " << mSize[3] << std::endl;
  header << "thicknesses:  NaN  NaN " << DoubleConvert( itkSpacing[2] ) << " NaN" << std::endl;
  // need to check
  header << "space directions: "
         << "(" << DoubleConvert( itkDirection[0][0] ) << "," << DoubleConvert( itkDirection[1][0] ) << ","
         << DoubleConvert( itkDirection[2][0] ) << ") "
         << "(" << DoubleConvert( itkDirection[0][1] ) << "," << DoubleConvert( itkDirection[1][1] ) << ","
         << DoubleConvert( itkDirection[2][1] ) << ") "
         << "(" << DoubleConvert( itkDirection[0][2] ) << "," << DoubleConvert( itkDirection[1][2] ) << ","
         << DoubleConvert( itkDirection[2][2] ) << ") none" << std::endl;
  header << "centerings: cell cell cell ???" << std::endl;
  header << "kinds: space space space list" << std::endl;

  header << "endian: little" << std::endl;
  header << "encoding: raw" << std::endl;
  header << "space units: \"mm\" \"mm\" \"mm\"" << std::endl;

  header << "space origin: "
         << "(" << DoubleConvert( itkOrigin[0] ) << "," << DoubleConvert( itkOrigin[1] ) << ","
         << DoubleConvert( itkOrigin[2] ) << ") " << std::endl;
  if ( !nrrdFormat )
  {
    header << "data file: " << itksys::SystemTools::GetFilenameName( dataFilename ) << std::endl;
  }

  {
    mxArray * mxMeasurementFrame = msm.GetField( "measurementframe" );
    if ( mxMeasurementFrame != 0 )
    {
      const mwSize * const measurementFrameSize = msm.GetDimensions( "data" );
      const double * const mxMeasurementFrame_temp = (double *)mxGetData( mxMeasurementFrame );
      using MeasurementMatType = std::vector< std::vector< double > >;
      MeasurementMatType measurementFrame;
      for ( unsigned int i = 0, count = 0; i < 3; ++i )
      {
        // if the saved measurement frame is N-d and weire at N+1-D...
        std::vector< double > tmp;
        if ( i >= measurementFrameSize[0] )
        {
          for ( unsigned int j = 0; j < 3; ++j )
          {
            if ( j == i )
            {
              tmp.push_back( 1.0 );
            }
            else
            {
              tmp.push_back( 0.0 );
            }
          }
        }
        else
        {
          for ( unsigned int j = 0; j < 3; ++j )
          {
            if ( j < measurementFrameSize[0] )
            {
              tmp.push_back( mxMeasurementFrame_temp[count] * flipFactors[j] );
              ++count;
            }
            else if ( j == i )
            {
              tmp.push_back( 1.0 );
            }
            else
            {
              tmp.push_back( 0.0 );
            }
          }
        }
        measurementFrame.push_back( tmp );
      }
      header << "measurement frame: "
             << "(" << DoubleConvert( measurementFrame[0][0] ) << "," << DoubleConvert( measurementFrame[1][0] ) << ","
             << DoubleConvert( measurementFrame[2][0] ) << ") "
             << "(" << DoubleConvert( measurementFrame[0][1] ) << "," << DoubleConvert( measurementFrame[1][1] ) << ","
             << DoubleConvert( measurementFrame[2][1] ) << ") "
             << "(" << DoubleConvert( measurementFrame[0][2] ) << "," << DoubleConvert( measurementFrame[1][2] ) << ","
             << DoubleConvert( measurementFrame[2][2] ) << ")" << std::endl;
    }
    else
    {
      header << "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl;
    }
  }

  header << "modality:=DWMRI" << std::endl;
  // this is the norminal BValue, i.e. the largest one.
  const double * const bvalue = static_cast< double * >( mxGetData( msm.GetField( "bvalue" ) ) );
  header << "DWMRI_b-value:=" << DoubleConvert( *bvalue ) << std::endl;

  {
    const mxArray * const gradientdirections = msm.GetField( "gradientdirections" );
    if ( gradientdirections != 0 )
    {
      std::stringstream ss;
      ss.precision( 10 );
      const unsigned numGradients = mxGetNumberOfElements( gradientdirections ) / 3;
      if ( numGradients > 0 )
      {
        const double * const gradients = static_cast< double * >( mxGetData( gradientdirections ) );
        for ( unsigned i = 0; i < numGradients; ++i )
        {
          header << "DWMRI_gradient_" << std::setw( 4 ) << std::setfill( '0' ) << i
                 << ":=" << DoubleConvert( gradients[i] ) << "   " << DoubleConvert( gradients[i + numGradients] )
                 << "   " << DoubleConvert( gradients[i + ( 2 * numGradients )] ) << std::endl;
        }
      }
    }
  }
  // write data in the same file is .nrrd was chosen
  header << std::endl;
  ;
  const mxArray * const dataMx = msm.GetField( "data" );
  TScalar * const       voxels = static_cast< TScalar * >( mxGetData( dataMx ) );
  if ( nrrdFormat )
  {
    header.write( reinterpret_cast< char * >( voxels ), numPixels * sizeof( TScalar ) );
  }
  else
  {
    std::ofstream datafile;
    datafile.open( dataFilename, std::ios::out | std::ios::binary );
    datafile.write( reinterpret_cast< char * >( voxels ), numPixels * sizeof( TScalar ) );
  }
  header.close();
} // end WriteDWINrrd

// itkSaveWithMetaData is called from matlab
void
itkSaveWithMetaData( int nrhs, const mxArray * prhs[] )
{
  assert( 2 == nrhs );
  const char                me[] = "itkSaveWithMetadata";
  char                      errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
  const mxArray * const     filenameMx = prhs[0];
  const mxArray * const     structMx = prhs[1];
  const MatlabStructManager msm( structMx );

  /* Metadata stuff */
  const int filenameLen = mxGetM( filenameMx ) * mxGetN( filenameMx ) + 1;
  /* managed by Matlab */
  char * filename = static_cast< char * >( mxCalloc( filenameLen, sizeof( mxChar ) ) );
  mxGetString( filenameMx, filename, filenameLen );

  /* Error checking on the data */
  if ( mxIsComplex( msm.GetField( "data" ) ) )
  {
    snprintf( errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: sorry, array must be real", me );
    mexErrMsgTxt( errBuff );
  }

  const itk::ImageIOBase::IOComponentType ntype = typeMtoITK( mxGetClassID( msm.GetField( "data" ) ) );

  if ( itk::ImageIOBase::UNKNOWNCOMPONENTTYPE == ntype )
  {
    snprintf( errBuff,
              NRRD_MAX_ERROR_MSG_SIZE,
              "%s: sorry, can't handle type %s",
              me,
              mxGetClassName( msm.GetField( "data" ) ) );
    mexErrMsgTxt( errBuff );
  }

  const unsigned int dim = msm.GetNumberOfDimensions( "data" );
  if ( dim < 1 || dim > 4 )
  {
    snprintf( errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: number of array dimensions %u outside range [1,%d]", me, dim, 4 );
    mexErrMsgTxt( errBuff );
  }

  switch ( dim )
  {
    case 2:
      switch ( ntype )
      {
        case itk::ImageIOBase::CHAR:
          WriteITKImageFromMatlabStructure< itk::Image< char, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::UCHAR:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned char, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::SHORT:
          WriteITKImageFromMatlabStructure< itk::Image< short, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::USHORT:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned short, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::INT:
          WriteITKImageFromMatlabStructure< itk::Image< int, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::UINT:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned int, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::LONG:
          WriteITKImageFromMatlabStructure< itk::Image< long, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::ULONG:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned long, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::FLOAT:
          WriteITKImageFromMatlabStructure< itk::Image< float, 2 > >( msm, filename );
          break;
        case itk::ImageIOBase::DOUBLE:
          WriteITKImageFromMatlabStructure< itk::Image< double, 2 > >( msm, filename );
          break;
        default:
          break;
      }
      break;

    case 3:
      switch ( ntype )
      {
        case itk::ImageIOBase::CHAR:
          WriteITKImageFromMatlabStructure< itk::Image< char, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::UCHAR:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned char, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::SHORT:
          WriteITKImageFromMatlabStructure< itk::Image< short, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::USHORT:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned short, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::INT:
          WriteITKImageFromMatlabStructure< itk::Image< int, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::UINT:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned int, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::LONG:
          WriteITKImageFromMatlabStructure< itk::Image< long, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::ULONG:
          WriteITKImageFromMatlabStructure< itk::Image< unsigned long, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::FLOAT:
          WriteITKImageFromMatlabStructure< itk::Image< float, 3 > >( msm, filename );
          break;
        case itk::ImageIOBase::DOUBLE:
          WriteITKImageFromMatlabStructure< itk::Image< double, 3 > >( msm, filename );
          break;
        default:
          myMexPrintf( " INVALID 3D TYPE SPECIFIED\n" );
          break;
      }
      break;

    case 4:
    {
      const mxArray * const gradientdirections = msm.GetField( "gradientdirections" );
      if ( !gradientdirections || mxGetNumberOfElements( gradientdirections ) < 1 )
      {
        switch ( ntype )
        {
          case itk::ImageIOBase::CHAR:
            WriteITKImageFromMatlabStructure< itk::Image< char, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::UCHAR:
            WriteITKImageFromMatlabStructure< itk::Image< unsigned char, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::SHORT:
            WriteITKImageFromMatlabStructure< itk::Image< short, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::USHORT:
            WriteITKImageFromMatlabStructure< itk::Image< unsigned short, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::INT:
            WriteITKImageFromMatlabStructure< itk::Image< int, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::UINT:
            WriteITKImageFromMatlabStructure< itk::Image< unsigned int, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::LONG:
            WriteITKImageFromMatlabStructure< itk::Image< long, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::ULONG:
            WriteITKImageFromMatlabStructure< itk::Image< unsigned long, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::FLOAT:
            WriteITKImageFromMatlabStructure< itk::Image< float, 4 > >( msm, filename );
            break;
          case itk::ImageIOBase::DOUBLE:
            WriteITKImageFromMatlabStructure< itk::Image< double, 4 > >( msm, filename );
            break;
          default:
            myMexPrintf( " INVALID 4D TYPE SPECIFIED\n" );
            break;
        }
      }
      else // WRITE DWI DATA
      {
        switch ( ntype )
        {
          case itk::ImageIOBase::CHAR:
            WriteDWINrrd< char >( msm, filename, "int8" );
            break;
          case itk::ImageIOBase::UCHAR:
            WriteDWINrrd< unsigned char >( msm, filename, "uchar" );
            break;
          case itk::ImageIOBase::SHORT:
            WriteDWINrrd< short >( msm, filename, "short" );
            break;
          case itk::ImageIOBase::USHORT:
            WriteDWINrrd< unsigned short >( msm, filename, "ushort" );
            break;
          case itk::ImageIOBase::INT:
            WriteDWINrrd< int >( msm, filename, "int" );
            break;
          case itk::ImageIOBase::UINT:
            WriteDWINrrd< unsigned int >( msm, filename, "uint" );
            break;
          case itk::ImageIOBase::FLOAT:
            WriteDWINrrd< float >( msm, filename, "float" );
            break;
          case itk::ImageIOBase::DOUBLE:
            WriteDWINrrd< double >( msm, filename, "double" );
            break;
          default:
            myMexPrintf( " INVALID 4D DWI Nrrd TYPE SPECIFIED\n" );
            break;
        }
      }
    }
    break;

    default:
      break;
  }
} // end itkSaveWithMetaData

void
mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] )
{
  // check input parameters
  if ( 0 != nlhs || 2 != nrhs || NULL == mxGetData( prhs[0] ) || !mxIsChar( prhs[0] ) )
  {
    std::cerr << "Parameters Error in using itkSaveWithMetadata" << std::endl;
    std::cerr << "Usage: itkSaveWithMetadata('outputFilename', matlabDataStruct)" << std::endl;
    return;
  }

  try
  {
    itkSaveWithMetaData( nrhs, prhs );
  }
  catch ( std::exception & e )
  {
    mexErrMsgTxt( e.what() );
  }
} // end mexFunction
