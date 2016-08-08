#include "itkImage.h"
#include "itksys/SystemTools.hxx"
#include "itkMetaDataObject.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageIOFactoryRegisterManager.h"
#include "nrrdCommon.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "vnl/vnl_vector_fixed.h"
#include "Debug.h"
#include "itkNumberToString.h"

bool debug = false;

void myMexPrintf(  std::string msg){
  if(debug)
    mexPrintf( msg.c_str());
}void myMexPrintf(  std::string msg, int value){
  if(debug)
    mexPrintf( msg.c_str(), value);
}void myMexPrintf(  std::string msg, std::string value){
  if(debug)
    mexPrintf( msg.c_str(), value.c_str());
}void myMexPrintf( std::string msg, double value){
  if(debug)
    mexPrintf( msg.c_str(), value);
}void myMexPrintf( std::string msg, unsigned int value){
  if(debug)
    mexPrintf( msg.c_str(), value);
}void myMexPrintf(  std::string msg, int value, int value2, int value3){
  if(debug)
    mexPrintf( msg.c_str(), value, value2, value3);
}void myMexPrintf(  std::string msg, double value, int value2, int value3){
  if(debug)
    mexPrintf( msg.c_str(), value, value2, value3);
}

inline
itk::ImageIOBase::IOComponentType typeMtoITK(const mxClassID mtype)
{

  switch( mtype ){
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
  }//end switch
  return itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;
}//end typeMtoITK

/** need a specialization that only calls
* SetNumberOfComponentsPerPixel for VectorImages
 */
typedef itk::VectorImage<double, 3> DWIImage;

template <typename TImage>
void SetNumberOfComponentsPerPixel(typename TImage::Pointer, unsigned)
{
}

template
<>
void SetNumberOfComponentsPerPixel<DWIImage>(DWIImage::Pointer im, unsigned numComponents)
{
  im->SetNumberOfComponentsPerPixel(numComponents);
}

template <typename TImage>
void
WriteFile(const MatlabStructManager & msm, const char *filename)
{
  typedef TImage                          ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  const unsigned int   numDims = msm.GetNumberOfDimensions("data");
  const mwSize * const mSize = msm.GetDimensions("data");
  myMexPrintf( " size =\n");

  typename ImageType::SizeType itkSize;
  size_t numPixels = 1;
  {
    unsigned int axIdx;
    for( axIdx = 0; axIdx < ImageType::ImageDimension; axIdx++ ){
      itkSize[axIdx] = mSize[axIdx];
      numPixels *= itkSize[axIdx];
      myMexPrintf( "%d, %d, %d\n",(int)itkSize[axIdx], (int)axIdx, (int)ImageType::ImageDimension);
    }//end for
    if( axIdx < numDims ){
      numPixels *= mSize[axIdx];
    }//end if
  }

  /** spaceorigin **/
  typename ImageType::PointType itkOrigin;
  itkOrigin.Fill(0.0);
  const double *spaceorigin_temp = (double *)mxGetData( msm.GetField("spaceorigin") );
  // TODO:  Make work for 2D, but currently only works for 3D and 3D vectors.
  const unsigned int spatialDims = 3; // msm.GetNumberOfDimensions("data");
  if( spatialDims != 3 ){
    myMexPrintf( "ERROR: ONLY 3D images (3D+gradients OK) supported.\n");
    return;
  }//end if
  myMexPrintf( "Space Origin =\n");
  for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++ ){
    if( sdIdx < spatialDims ){
      itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
    }else{
      itkOrigin[sdIdx] = 0.0;
    }//end else
    myMexPrintf( "%d, %d, %d\n",(int)itkOrigin[sdIdx], (int)sdIdx, (int)ImageType::ImageDimension);
  }//end for

  typename ImageType::SpacingType itkSpacing;
  typename ImageType::DirectionType itkDirection;
  //
  // fill out extra dimension for (usually) 4 D -- add zero to 4th
  // entry in direction vector, except for the last vector, which
  // needs to be 0 0 0 1
  itkDirection.SetIdentity();
  /** spacedirections **/
  const double *spacedirections_temp = (double *)mxGetData( msm.GetField("spacedirections") );
  for( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx ){
    vnl_vector_fixed<double, spatialDims> vec;
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx ){
      const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
      vec[sdIdx] = spacedirections_temp[sdir_offset];
    }//end for
    itkSpacing[axIdx] = vec.magnitude();
    vec.normalize();
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx ){
      itkDirection[sdIdx][axIdx] = vec[sdIdx];
    }//end for
  }//end for
  myMexPrintf( "directions =\n");
  for( unsigned int axIdx = 0; axIdx < ImageType::ImageDimension; ++axIdx ){
    for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx ){
      myMexPrintf( "%lf, %d, %d\n",itkDirection[sdIdx][axIdx], sdIdx, ImageType::ImageDimension);
    }//end for
  }//end for

  double      flipFactors[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int * nrrdSpaceDefinition = (int *)mxGetData( msm.GetField("space") );
  if( nrrdSpaceDefinition ){
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch( *nrrdSpaceDefinition ){
      case nrrdSpaceLeftPosteriorSuperior:
        // already in ITK space format
        break;
      case nrrdSpaceRightAnteriorSuperior:
        for( unsigned int i = 0; i < ImageType::ImageDimension; ++i ){
          itkDirection[i][0] *= -1;
          itkDirection[i][1] *= -1;
        }//end for
        flipFactors[0] = -1.0;
        flipFactors[1] = -1.0;
        break;
      case nrrdSpaceLeftAnteriorSuperior:
        for( unsigned int i = 0; i < ImageType::ImageDimension; ++i ){
          itkDirection[i][0] *= -1;
        }//end for
        flipFactors[0] = -1.0;
        break;
    }//end switch
  }//end if

  typename ImageType::Pointer im = ImageType::New();
  im->SetDirection(itkDirection);
  im->SetOrigin(itkOrigin);
  im->SetSpacing(itkSpacing);
  im->SetRegions(itkSize);
  size_t ComponentsPerPixel=1;
  //
  // if we have a vector image
  if( numDims > ImageType::ImageDimension ){
    ComponentsPerPixel = mSize[ImageType::ImageDimension];
    SetNumberOfComponentsPerPixel<ImageType>(im, ComponentsPerPixel);
  }//end if
  im->Allocate();
  const mxArray * const dataMx = msm.GetField("data");
  //Note: Matlab returns the internal value types, not a vector of those types.
  const typename TImage::InternalPixelType * voxels =
    static_cast<typename TImage::InternalPixelType *>( mxGetData( dataMx ) );
  {
    std::copy(voxels, voxels + numPixels*ComponentsPerPixel, im->GetBufferPointer() );
  }

  itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();
  //
  // space direction
  std::string spaceDirKey("NRRD_");
  spaceDirKey += airEnumStr(nrrdSpace, nrrdField_space);
  itk::EncapsulateMetaData<std::string>(thisDic, spaceDirKey,
                                        std::string( airEnumStr(nrrdSpace, nrrdSpaceLeftPosteriorSuperior) ) );

  //
  // fill out metadata
  // Measurement Frame
  mxArray *mxMeasurementFrame = msm.GetField("measurementframe");
  /** measurementframe **/
  if( mxMeasurementFrame ){
    const double * const mxMeasurementFrame_temp = (double *)mxGetData( mxMeasurementFrame );
    const mwSize * const measurementFrameSize = msm.GetDimensions("data");
    if( mxMeasurementFrame != 0 ){
      typedef std::vector<std::vector<double> > MeasurementMatType;
      MeasurementMatType measurementFrame;
      std::string        measurementFrameFieldName = "NRRD_";
      measurementFrameFieldName += airEnumStr(nrrdField, nrrdField_measurement_frame);
      for( unsigned int i = 0, count = 0; i < ImageType::ImageDimension; ++i ){
        // if the saved measurement frame is N-d and weire at N+1-D...
        std::vector<double> tmp;
        if( i >= measurementFrameSize[0] ){
          for( unsigned int j = 0; j < ImageType::ImageDimension; ++j ){
            if( j == i ){
              tmp.push_back(1.0);
            }else{
              tmp.push_back(0.0);
            }//end else
          }//end for
        }else{
          for( unsigned int j = 0; j < ImageType::ImageDimension; ++j ){
            if( j < measurementFrameSize[0] ){
              tmp.push_back(mxMeasurementFrame_temp[count] * flipFactors[j]);
              ++count;
            }else if( j == i ){
              tmp.push_back(1.0);
            }else{
              tmp.push_back(0.0);
            }//end else
          }//end for
        }//end else
        measurementFrame.push_back(tmp);
      }//end for
      itk::EncapsulateMetaData<MeasurementMatType>(thisDic, measurementFrameFieldName, measurementFrame);
    }//end if
  }//end if
  //
  // centerings.
  const mxArray * const centerings_temp_MxArray = msm.GetField("centerings");
  if( centerings_temp_MxArray ){
    const int *centerings_temp = (int *)mxGetData( centerings_temp_MxArray );
    for( unsigned int i = 0; i < numDims; ++i ){
      std::stringstream attName;
      attName << "NRRD_centerings[" << i << "]";
      std::string val = airEnumStr(nrrdCenter, centerings_temp[i]);
      itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
    }//end for
  }//end if
  // kinds
  const mxArray * const kinds_temp_MxArray = msm.GetField("kinds");
  if( kinds_temp_MxArray )
    {
    const int * const kinds_temp = (int *)mxGetData( kinds_temp_MxArray );
    for( unsigned int axIdx = 0; axIdx < numDims; ++axIdx )
      {
      std::stringstream attName;
      attName << "NRRD_kinds[" << axIdx << "]";
      std::string val = airEnumStr(nrrdKind, kinds_temp[axIdx]);
      itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
      }
    }

  // space units
  itk::EncapsulateMetaData<std::string>(thisDic, std::string("space units"),
                                        std::string("\"mm\" \"mm\" \"mm\"") );

  if( msm.isDWIdata() )
    {
    // One of the directions is the gradient list instead of a spaceDim
      { // Get the DWI specific fields.
      const mxArray * const gradientdirections = msm.GetField("gradientdirections");
      if( gradientdirections != 0 )
        {
        std::stringstream ss;
        ss.precision(10);
        const unsigned numGradients = mxGetNumberOfElements(gradientdirections) / 3;
        if( numGradients > 0 )
          {
          const double * const gradients = static_cast<double *>(mxGetData(gradientdirections) );
          for( unsigned i = 0; i < numGradients; ++i )
            {
            ss.str(std::string() );
            ss << gradients[i] << " "
               << gradients[i + numGradients] << " "
               << gradients[i + (2 * numGradients)];
            std::stringstream gradName;
            gradName << GRADIENT_PREFIX << std::setw(4) << std::setfill('0') << i;
            itk::EncapsulateMetaData<std::string>(thisDic,
                                                  gradName.str(), ss.str() );
            }
          itk::EncapsulateMetaData<std::string>(thisDic, std::string("modality"),
                                                std::string("DWMRI") );
          const double * const bvalue = static_cast<double *>(mxGetData(msm.GetField("bvalue") ) );
          ss << *bvalue;
          itk::EncapsulateMetaData<std::string>(thisDic, std::string("DWMRI_b-value"),
                                                ss.str() );
          }
        }
      // OK
      }
    }//end if msm.isDWIdata

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(im);
  try
    {
    writer->Write();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::string msg("itkSaveWithMetaData: can't read\n");
    msg += filename;
    msg += " ";
    msg += excp.what();
    char errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: at line %d", msg.c_str(), __LINE__);
    mexErrMsgTxt(errBuff);
    }
}// end WriteFile

template <typename TImage>
void
WriteDWIFile(const MatlabStructManager & msm, const char *filename)
{
  typedef TImage                          ImageType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  const unsigned int   numDims = msm.GetNumberOfDimensions("data");
  const mwSize * const mSize = msm.GetDimensions("data");
  myMexPrintf( " size =\n");

  typename ImageType::SizeType itkSize;
  size_t numPixels = 1;
    {
    unsigned int axIdx;
    for( axIdx = 0; axIdx < ImageType::ImageDimension; axIdx++ )
      {
      itkSize[axIdx] = mSize[axIdx];
      numPixels *= itkSize[axIdx];
      myMexPrintf( "%lf, %d, %d\n",itkSize[axIdx], axIdx, ImageType::ImageDimension);
      }
    if( axIdx < numDims )
      {
      numPixels *= mSize[axIdx];
      }
    }

  /** spaceorigin **/
  typename ImageType::PointType itkOrigin;
  itkOrigin.Fill(0.0);
  const double *spaceorigin_temp = (double *)mxGetData( msm.GetField("spaceorigin") );
  // TODO:  Make work for 2D, but currently only works for 3D and 3D vectors.
  const unsigned int spatialDims = 3; // msm.GetNumberOfDimensions("data");
  if( spatialDims != 3 )
    {
    myMexPrintf( "ERROR: ONLY 3D images (3D+gradients OK) supported.\n");
    return;
    }
  myMexPrintf( "Space Origin =\n");
  for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++ )
    {
    if( sdIdx < spatialDims )
      {
      itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
      }
    else
      {
      itkOrigin[sdIdx] = 0.0;
      }
    myMexPrintf( "%lf, %d, %d\n",itkOrigin[sdIdx], sdIdx, ImageType::ImageDimension);
    }

  typename ImageType::SpacingType itkSpacing;
  typename ImageType::DirectionType itkDirection;
  //
  // fill out extra dimension for (usually) 4 D -- add zero to 4th
  // entry in direction vector, except for the last vector, which
  // needs to be 0 0 0 1
  itkDirection.SetIdentity();
  /** spacedirections **/
  const double *spacedirections_temp = (double *)mxGetData( msm.GetField("spacedirections") );
  for( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx )
    {
    vnl_vector_fixed<double, spatialDims> vec;
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
      {
      const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
      vec[sdIdx] = spacedirections_temp[sdir_offset];
      }
    itkSpacing[axIdx] = vec.magnitude();
    vec.normalize();
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
      {
      itkDirection[sdIdx][axIdx] = vec[sdIdx];
      }
    }
  myMexPrintf( "directions =\n");
  for( unsigned int axIdx = 0; axIdx < ImageType::ImageDimension; ++axIdx )
    {
    for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx )
      {
      myMexPrintf( "%lf, %d, %d\n",itkDirection[sdIdx][axIdx], sdIdx, ImageType::ImageDimension);
      }
    }

  double      flipFactors[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int * nrrdSpaceDefinition = (int *)mxGetData( msm.GetField("space") );
  if( nrrdSpaceDefinition )
    {
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch( *nrrdSpaceDefinition )
      {
      case nrrdSpaceLeftPosteriorSuperior:
        // already in ITK space format
        break;
      case nrrdSpaceRightAnteriorSuperior:
        for( unsigned int i = 0; i < ImageType::ImageDimension; ++i )
          {
          itkDirection[i][0] *= -1;
          itkDirection[i][1] *= -1;
          }
        flipFactors[0] = -1.0;
        flipFactors[1] = -1.0;
        break;
      case nrrdSpaceLeftAnteriorSuperior:
        for( unsigned int i = 0; i < ImageType::ImageDimension; ++i )
          {
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
  size_t ComponentsPerPixel=1;
  //
  // if we have a vector image
  if( numDims > ImageType::ImageDimension )
    {
    ComponentsPerPixel = mSize[ImageType::ImageDimension];
    SetNumberOfComponentsPerPixel<ImageType>(im, ComponentsPerPixel);
    }
  im->Allocate();
  const mxArray * const dataMx = msm.GetField("data");
  //Note: Matlab returns the internal value types, not a vector of those types.
  const typename TImage::InternalPixelType * voxels =
    static_cast<typename TImage::InternalPixelType *>( mxGetData( dataMx ) );
  if( numDims > ImageType::ImageDimension )
    {
    itk::ImageRegionIterator<TImage> it(im,im->GetLargestPossibleRegion());
    while(!it.IsAtEnd())
      {
      itk::VariableLengthVector<typename TImage::InternalPixelType> vlv(voxels,ComponentsPerPixel,false);
      voxels+=ComponentsPerPixel;
      it.Set(vlv);
      }
    }

  itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();
  //
  // space direction
  std::string spaceDirKey("NRRD_space");
  spaceDirKey += airEnumStr(nrrdSpace, nrrdField_space);

  myMexPrintf( "%s\n",spaceDirKey.c_str());
  itk::EncapsulateMetaData<std::string>(thisDic, spaceDirKey,
                                        std::string( airEnumStr(nrrdSpace, nrrdSpaceLeftPosteriorSuperior) ) );

  //
  // fill out metadata
  // Measurement Frame
  mxArray *mxMeasurementFrame = msm.GetField("measurementframe");
  /** measurementframe **/
  if( mxMeasurementFrame )
    {
    const double * const mxMeasurementFrame_temp = (double *)mxGetData( mxMeasurementFrame );
    const mwSize * const measurementFrameSize = msm.GetDimensions("data");
    if( mxMeasurementFrame != 0 )
      {
      typedef std::vector<std::vector<double> > MeasurementMatType;
      MeasurementMatType measurementFrame;
      std::string        measurementFrameFieldName = "NRRD_";
      measurementFrameFieldName += airEnumStr(nrrdField, nrrdField_measurement_frame);
      for( unsigned int i = 0, count = 0; i < ImageType::ImageDimension; ++i )
        {
        // if the saved measurement frame is N-d and weire at N+1-D...
        std::vector<double> tmp;
        if( i >= measurementFrameSize[0] )
          {
          for( unsigned int j = 0; j < ImageType::ImageDimension; ++j )
            {
            if( j == i )
              {
              tmp.push_back(1.0);
              }
            else
              {
              tmp.push_back(0.0);
              }
            }
          }
        else
          {
          for( unsigned int j = 0; j < ImageType::ImageDimension; ++j )
            {
            if( j < measurementFrameSize[0] )
              {
              tmp.push_back(mxMeasurementFrame_temp[count] * flipFactors[j]);
              ++count;
              }
            else if( j == i )
              {
              tmp.push_back(1.0);
              }
            else
              {
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
  const mxArray * const centerings_temp_MxArray = msm.GetField("centerings");
  if( centerings_temp_MxArray )
    {
    const int *centerings_temp = (int *)mxGetData( centerings_temp_MxArray );
    for( unsigned int i = 0; i < numDims; ++i )
      {
      std::stringstream attName;
      attName << "NRRD_centerings[" << i << "]";
      std::string val = airEnumStr(nrrdCenter, centerings_temp[i]);
      itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
      }
    }
  // kinds
  const mxArray * const kinds_temp_MxArray = msm.GetField("kinds");
  if( kinds_temp_MxArray )
    {
    const int * const kinds_temp = (int *)mxGetData( kinds_temp_MxArray );
    for( unsigned int axIdx = 0; axIdx < numDims; ++axIdx )
      {
      std::stringstream attName;
      attName << "NRRD_kinds[" << axIdx << "]";
      std::string val = airEnumStr(nrrdKind, kinds_temp[axIdx]);
      itk::EncapsulateMetaData<std::string>(thisDic, attName.str(), val);
      }
    }

  // space units
  itk::EncapsulateMetaData<std::string>(thisDic, std::string("space units"),
                                        std::string("\"mm\" \"mm\" \"mm\"") );

  if( msm.isDWIdata() )
    {
    // One of the directions is the gradient list instead of a spaceDim
      { // Get the DWI specific fields.
      const mxArray * const gradientdirections = msm.GetField("gradientdirections");
      if( gradientdirections != 0 )
        {
        std::stringstream ss;
        ss.precision(10);
        const unsigned numGradients = mxGetNumberOfElements(gradientdirections) / 3;
        if( numGradients > 0 )
          {
          const double * const gradients = static_cast<double *>(mxGetData(gradientdirections) );
          for( unsigned i = 0; i < numGradients; ++i )
            {
            ss.str(std::string() );
            ss << gradients[i] << " "
               << gradients[i + numGradients] << " "
               << gradients[i + (2 * numGradients)];
            std::stringstream gradName;
            gradName << GRADIENT_PREFIX << std::setw(4) << std::setfill('0') << i;
            itk::EncapsulateMetaData<std::string>(thisDic,
                                                  gradName.str(), ss.str() );
            }
          itk::EncapsulateMetaData<std::string>(thisDic, std::string("modality"),
                                                std::string("DWMRI") );
          const double * const bvalue = static_cast<double *>(mxGetData(msm.GetField("bvalue") ) );
          ss << *bvalue;
          itk::EncapsulateMetaData<std::string>(thisDic, std::string("DWMRI_b-value"),
                                                ss.str() );
          }
        }
      // OK
      }
    }

  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(im);
  try
    {
    writer->Write();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::string msg("itkSaveWithMetaData: can't read\n");
    msg += filename;
    msg += " ";
    msg += excp.what();
    char errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: at line %d", msg.c_str(), __LINE__);
    mexErrMsgTxt(errBuff);
    }
}// end WriteDWIFile


template <typename TScalar>
void
WriteDWINrrd(const MatlabStructManager & msm, const char *filename, const char *voxelTypeName)
{
  typedef typename itk::Image<TScalar,4> ImageType;
  std::string _fname(filename);
  std::string ext = itksys::SystemTools::GetFilenameLastExtension(_fname);
  std::string dataFilename;
  bool nrrdFormat(true);
  if(ext == ".nhdr")
    {
    nrrdFormat = false;
    dataFilename = itksys::SystemTools::GetFilenameWithoutLastExtension(_fname);
    dataFilename += ".raw"; // dump out uncompressed.
    }
  else if(ext != ".nrrd")
    {
    std::string msg = "Only NRRD file types recognized, extension = ";
    msg += ext;
    mexErrMsgTxt(msg.c_str());
    return;
    }

  const unsigned int   numDims = msm.GetNumberOfDimensions("data");
  const mwSize * const mSize = msm.GetDimensions("data");
  myMexPrintf( " size =\n");

  size_t numPixels = 1;
    {
    unsigned int axIdx;
    for( axIdx = 0; axIdx < numDims; axIdx++ )
      {
      numPixels *= mSize[axIdx];
      myMexPrintf( "%d, %d, %d\n", (int)mSize[axIdx], (int)axIdx, (int)numDims);
      }
    if( axIdx < numDims )
      {
      numPixels *= mSize[axIdx];
      }
    }

  /** spaceorigin **/
  typename ImageType::PointType itkOrigin;
  itkOrigin.Fill(0.0);
  const double *spaceorigin_temp = (double *)mxGetData( msm.GetField("spaceorigin") );
  // TODO:  Make work for 2D, but currently only works for 3D and 3D vectors.
  const unsigned int spatialDims = 3; // msm.GetNumberOfDimensions("data");
  if( spatialDims != 3 )
    {
    myMexPrintf( "ERROR: ONLY 3D images (3D+gradients OK) supported.\n");
    return;
    }
  myMexPrintf( "Space Origin =\n");
  for( unsigned int sdIdx = 0; sdIdx < spatialDims; sdIdx++ )
    {
    if( sdIdx < spatialDims )
      {
      itkOrigin[sdIdx] = spaceorigin_temp[sdIdx];
      }
    myMexPrintf( "%d, %d, %d\n",(int)itkOrigin[sdIdx], (int)sdIdx, (int)ImageType::ImageDimension);
    }

  typename ImageType::SpacingType itkSpacing;
  typename ImageType::DirectionType itkDirection;
  //
  // fill out extra dimension for (usually) 4 D -- add zero to 4th
  // entry in direction vector, except for the last vector, which
  // needs to be 0 0 0 1
  itkDirection.SetIdentity();
  /** spacedirections **/
  const double *spacedirections_temp = (double *)mxGetData( msm.GetField("spacedirections") );
  for( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx )
    {
    vnl_vector_fixed<double, spatialDims> vec;
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
      {
      const unsigned int sdir_offset = axIdx * spatialDims + sdIdx;
      vec[sdIdx] = spacedirections_temp[sdir_offset];
      }
    itkSpacing[axIdx] = vec.magnitude();
    vec.normalize();
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
      {
      itkDirection[sdIdx][axIdx] = vec[sdIdx];
      }
    }
  myMexPrintf( "directions =\n");
  for( unsigned int axIdx = 0; axIdx < spatialDims; ++axIdx )
    {
    for( unsigned int sdIdx = 0; sdIdx < spatialDims; ++sdIdx )
      {
      myMexPrintf( "%lf, %d, %d\n",itkDirection[sdIdx][axIdx], sdIdx, spatialDims);
      }
    }

  double      flipFactors[4] = { 1.0, 1.0, 1.0, 1.0 };
  const int * nrrdSpaceDefinition = (int *)mxGetData( msm.GetField("space") );
  if( nrrdSpaceDefinition )
    {
    // only need to adjust if an image was read in using the nrrd
    // reader, and is not LPS.
    switch( *nrrdSpaceDefinition )
      {
      case nrrdSpaceLeftPosteriorSuperior:
        // already in ITK space format
        break;
      case nrrdSpaceRightAnteriorSuperior:
        for( unsigned int i = 0; i < spatialDims; ++i )
          {
          itkDirection[i][0] *= -1;
          itkDirection[i][1] *= -1;
          }
        flipFactors[0] = -1.0;
        flipFactors[1] = -1.0;
        break;
      case nrrdSpaceLeftAnteriorSuperior:
        for( unsigned int i = 0; i < spatialDims; ++i )
          {
          itkDirection[i][0] *= -1;
          }
        flipFactors[0] = -1.0;
        break;
      }
    }
  std::ofstream header;
  // std::string headerFileName = outputDir + "/" + outputFileName;

  header.open(filename, std::ios::out | std::ios::binary);
  header << "NRRD0005" << std::endl
         << std::setprecision(17) << std::scientific;

  // stamp with DWIConvert branding
  header << "# This file was created by itkSaveMetaData version 1.0" << std::endl;
  header << "#" << std::endl << "#" << std::endl;

  if( !nrrdFormat )
    {
    header << "content: exists(" << itksys::SystemTools::GetFilenameName(dataFilename) << ",0)"
           << std::endl;
    }
  header << "type: " << voxelTypeName << std::endl;
  header << "dimension: 4" << std::endl;
  header << "space: ";
  if(!nrrdSpaceDefinition || *nrrdSpaceDefinition == nrrdSpaceLeftPosteriorSuperior)
    {
    header << "left-posterior-superior" << std::endl;
    }
  else if(*nrrdSpaceDefinition == nrrdSpaceRightAnteriorSuperior)
    {
    header << "right-anterior-superior" << std::endl;
    }
  else if(*nrrdSpaceDefinition == nrrdSpaceLeftAnteriorSuperior)
    {
    header << "left-anterior-superior" << std::endl;
    }

  itk::NumberToString<double> DoubleConvert;

  header << "sizes: " << mSize[0]
         << " " << mSize[1]
         << " " << mSize[2]
         << " " << mSize[3] << std::endl;
  header << "thicknesses:  NaN  NaN " << DoubleConvert(itkSpacing[2]) << " NaN" << std::endl;
  // need to check
  header << "space directions: "
         << "("
         << DoubleConvert(itkDirection[0][0]) << ","
         << DoubleConvert(itkDirection[1][0]) << ","
         << DoubleConvert(itkDirection[2][0])
         << ") "
         << "("
         << DoubleConvert(itkDirection[0][1]) << ","
         << DoubleConvert(itkDirection[1][1]) << ","
         << DoubleConvert(itkDirection[2][1]) << ") "
         << "("
         << DoubleConvert(itkDirection[0][2]) << ","
         << DoubleConvert(itkDirection[1][2]) << ","
         << DoubleConvert(itkDirection[2][2])
         << ") none" << std::endl;
  header << "centerings: cell cell cell ???" << std::endl;
  header << "kinds: space space space list" << std::endl;

  header << "endian: little" << std::endl;
  header << "encoding: raw" << std::endl;
  header << "space units: \"mm\" \"mm\" \"mm\"" << std::endl;

  header << "space origin: "
         << "(" << DoubleConvert(itkOrigin[0])
         << "," << DoubleConvert(itkOrigin[1])
         << "," << DoubleConvert(itkOrigin[2]) << ") " << std::endl;
  if( !nrrdFormat )
    {
    header << "data file: " << itksys::SystemTools::GetFilenameName(dataFilename) << std::endl;
    }

    {
    mxArray *mxMeasurementFrame = msm.GetField("measurementframe");
    if(mxMeasurementFrame != 0)
      {
      const mwSize * const measurementFrameSize = msm.GetDimensions("data");
      const double * const mxMeasurementFrame_temp = (double *)mxGetData( mxMeasurementFrame );
      typedef std::vector<std::vector<double> > MeasurementMatType;
      MeasurementMatType measurementFrame;
      for( unsigned int i = 0, count = 0; i < 3; ++i )
        {
        // if the saved measurement frame is N-d and weire at N+1-D...
        std::vector<double> tmp;
        if( i >= measurementFrameSize[0] )
          {
          for( unsigned int j = 0; j < 3; ++j )
            {
            if( j == i )
              {
              tmp.push_back(1.0);
              }
            else
              {
              tmp.push_back(0.0);
              }
            }
          }
        else
          {
          for( unsigned int j = 0; j < 3; ++j )
            {
            if( j < measurementFrameSize[0] )
              {
              tmp.push_back(mxMeasurementFrame_temp[count] * flipFactors[j]);
              ++count;
              }
            else if( j == i )
              {
              tmp.push_back(1.0);
              }
            else
              {
              tmp.push_back(0.0);
              }
            }
          }
        measurementFrame.push_back(tmp);
        }
      header << "measurement frame: "
             << "(" << DoubleConvert(measurementFrame[0][0]) << ","
             << DoubleConvert(measurementFrame[1][0]) << ","
             << DoubleConvert(measurementFrame[2][0]) << ") "
             << "(" << DoubleConvert(measurementFrame[0][1]) << ","
             << DoubleConvert(measurementFrame[1][1]) << ","
             << DoubleConvert(measurementFrame[2][1]) << ") "
             << "(" << DoubleConvert(measurementFrame[0][2]) << ","
             << DoubleConvert(measurementFrame[1][2]) << ","
             << DoubleConvert(measurementFrame[2][2]) << ")"
             << std::endl;
      }
    else
      {
      header << "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl;
      }
    }

  header << "modality:=DWMRI" << std::endl;
  // this is the norminal BValue, i.e. the largest one.
  const double * const bvalue = static_cast<double *>(mxGetData(msm.GetField("bvalue") ) );
  header << "DWMRI_b-value:=" << DoubleConvert(*bvalue) << std::endl;

  {
  const mxArray * const gradientdirections = msm.GetField("gradientdirections");
  if( gradientdirections != 0 )
    {
    std::stringstream ss;
    ss.precision(10);
    const unsigned numGradients = mxGetNumberOfElements(gradientdirections) / 3;
    if( numGradients > 0 )
      {
      const double * const gradients = static_cast<double *>(mxGetData(gradientdirections) );
      for( unsigned i = 0; i < numGradients; ++i )
        {
        header << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << i << ":="
               << DoubleConvert(gradients[i]) << "   "
               << DoubleConvert(gradients[i + numGradients]) << "   "
               << DoubleConvert(gradients[i + (2 * numGradients)])
               << std::endl;
        }
      }
    }
  }
  // write data in the same file is .nrrd was chosen
  header << std::endl;;
  const mxArray * const dataMx = msm.GetField("data");
  TScalar * const voxels =
    static_cast<TScalar *>( mxGetData( dataMx ) );
  if( nrrdFormat )
    {
    header.write( reinterpret_cast<char *>(voxels),
                  numPixels * sizeof(TScalar) );
    }
  else
    {
    std::ofstream datafile;
    datafile.open(dataFilename, std::ios::out | std::ios::binary);
    datafile.write(reinterpret_cast<char *>(voxels),
                  numPixels * sizeof(TScalar) );
    }
  header.close();
}//end WriteDWINrrd

// itkSaveWithMetaData is called from matlab
void itkSaveWithMetaData(int, mxArray *[],
                         int nrhs, const mxArray *prhs[])
{
  const char            me[] = "itkSaveWithMetadata";
  char                  errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
  const mxArray * const filenameMx = prhs[0];

  if( !(2 == nrhs && mxIsChar(filenameMx) ) )
    {
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: requires two args: one string, one struct", me);
    mexErrMsgTxt(errBuff);
    return;
    }
  const mxArray * const     structMx = prhs[1];
  const MatlabStructManager msm(structMx);

  /* Metadata stuff */
  const int filenameLen = mxGetM(filenameMx) * mxGetN(filenameMx) + 1;
  /* managed by Matlab */
  char *filename = static_cast<char *>(mxCalloc(filenameLen, sizeof(mxChar) ) );
  mxGetString(filenameMx, filename, filenameLen);

  /* Error checking on the data */
  if( mxIsComplex(msm.GetField("data") ) )
    {
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: sorry, array must be real", me);
    mexErrMsgTxt(errBuff);
    }

  const itk::ImageIOBase::IOComponentType ntype =
    typeMtoITK(mxGetClassID(msm.GetField("data") ) );

  //myMexPrintf( "type index =\n");  myMexPrintf( ntype.toString());
  if( itk::ImageIOBase::UNKNOWNCOMPONENTTYPE == ntype )
    {
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: sorry, can't handle type %s",
             me, mxGetClassName(msm.GetField("data") ) );
    mexErrMsgTxt(errBuff);
    }

  const unsigned int dim = msm.GetNumberOfDimensions("data");
  if( dim < 1 || dim > 4 )
    {
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: number of array dimensions %u outside range [1,%d]",
             me, dim, 4);
    mexErrMsgTxt(errBuff);
    }
  myMexPrintf( "dim = %d\n",dim);
  myMexPrintf( "===============================================================\n");
  //myMexPrintf( ntype.toString());
  switch( dim )
    {
#if 0 // NOT SUPPORTING WRITING OF 2D yet, not tested
    case 2:

      switch( ntype )
        {
        case itk::ImageIOBase::CHAR:
          WriteFile<itk::Image<char, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::UCHAR:
          WriteFile<itk::Image<unsigned char, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::SHORT:
          WriteFile<itk::Image<short, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::USHORT:
          WriteFile<itk::Image<unsigned short, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::INT:
          WriteFile<itk::Image<int, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::UINT:
          WriteFile<itk::Image<unsigned int, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::LONG:
          WriteFile<itk::Image<long, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::ULONG:
          WriteFile<itk::Image<unsigned long, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::FLOAT:
          WriteFile<itk::Image<float, 2> >(msm, filename);
          break;
        case itk::ImageIOBase::DOUBLE:
          WriteFile<itk::Image<double, 2> >(msm, filename);
          break;
        default:
          break;
        }
      break;
#endif
    case 3:

      switch( ntype )
        {
        case itk::ImageIOBase::CHAR:
          myMexPrintf( " Writing 3D char\n");
          WriteFile<itk::Image<char, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::UCHAR:
          myMexPrintf( " Writing 3D uchar\n");
          WriteFile<itk::Image<unsigned char, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::SHORT:
          myMexPrintf( " Writing 3D short\n");
          WriteFile<itk::Image<short, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::USHORT:
          myMexPrintf( " Writing 3D ushort\n");
          WriteFile<itk::Image<unsigned short, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::INT:
          myMexPrintf( " Writing 3D int\n");
          WriteFile<itk::Image<int, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::UINT:
          myMexPrintf( " Writing 3D uint\n");
          WriteFile<itk::Image<unsigned int, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::LONG:
          myMexPrintf( " Writing 3D long\n");
          WriteFile<itk::Image<long, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::ULONG:
          myMexPrintf( " Writing 3D ulong\n");
          WriteFile<itk::Image<unsigned long, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::FLOAT:
          myMexPrintf( " Writing 3D float\n");
          WriteFile<itk::Image<float, 3> >(msm, filename);
          break;
        case itk::ImageIOBase::DOUBLE:
          myMexPrintf( " Writing 3D double\n");
          WriteFile<itk::Image<double, 3> >(msm, filename);
          break;
        default:
            myMexPrintf( " INVALID 3D TYPE SPECIFIED\n");
          break;
        }
      break;
    case 4:
      {
      const mxArray * const gradientdirections = msm.GetField("gradientdirections");
      if( !gradientdirections  || mxGetNumberOfElements(gradientdirections) < 1 )
        {
        switch( ntype )
          {
          case itk::ImageIOBase::CHAR:
            myMexPrintf( " Writing 4D char\n");
            WriteFile<itk::Image<char, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::UCHAR:
            myMexPrintf( " Writing 4D uchar\n");
            WriteFile<itk::Image<unsigned char, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::SHORT:
            myMexPrintf( " Writing 4D short\n");
            WriteFile<itk::Image<short, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::USHORT:
            myMexPrintf( " Writing 4D ushort\n");
            WriteFile<itk::Image<unsigned short, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::INT:
            myMexPrintf( " Writing 4D int\n");
            WriteFile<itk::Image<int, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::UINT:
            myMexPrintf( " Writing 4D uint\n");
            WriteFile<itk::Image<unsigned int, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::LONG:
            myMexPrintf( " Writing 4D long\n");
            WriteFile<itk::Image<long, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::ULONG:
            myMexPrintf( " Writing 4D ulong\n");
            WriteFile<itk::Image<unsigned long, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::FLOAT:
            myMexPrintf( " Writing 4D float\n");
            WriteFile<itk::Image<float, 4> >(msm, filename);
            break;
          case itk::ImageIOBase::DOUBLE:
            myMexPrintf( " Writing 4D double\n");
            WriteFile<itk::Image<double, 4> >(msm, filename);
            break;
          default:
            myMexPrintf( " INVALID 4D TYPE SPECIFIED\n");
            break;
          }
        }
      else
        {
        myMexPrintf( "WRITE DWI DATA");
        // WriteFile<itk::VectorImage<double,3> >(msm,filename);
        switch( ntype )
          {
          case itk::ImageIOBase::CHAR:
            myMexPrintf( " Writing 4D char\n");
            WriteDWINrrd<char>(msm, filename, "int8");
            break;
          case itk::ImageIOBase::UCHAR:
            myMexPrintf( " Writing 4D uchar\n");
            WriteDWINrrd<unsigned char >(msm, filename,"uchar");
            break;
          case itk::ImageIOBase::SHORT:
            myMexPrintf( " Writing 4D short\n");
            WriteDWINrrd<short >(msm, filename, "short");
            break;
          case itk::ImageIOBase::USHORT:
            myMexPrintf( " Writing 4D ushort\n");
            WriteDWINrrd<unsigned short>(msm, filename,"ushort");
            break;
          case itk::ImageIOBase::INT:
            myMexPrintf( " Writing 4D int\n");
            WriteDWINrrd<int>(msm, filename,"int");
            break;
          case itk::ImageIOBase::UINT:
            myMexPrintf( " Writing 4D uint\n");
            WriteDWINrrd<unsigned int>(msm, filename,"uint");
            break;
          case itk::ImageIOBase::FLOAT:
            myMexPrintf( " Writing 4D float\n");
            WriteDWINrrd<float>(msm, filename,"float");
            break;
          case itk::ImageIOBase::DOUBLE:
            myMexPrintf( " Writing 4D double\n");
            WriteDWINrrd<double>(msm, filename,"double");
            break;
          default:
            myMexPrintf( " INVALID 4D TYPE SPECIFIED\n");
            break;
          }
        }
      }
      break;
    default:
      break;
    }
  myMexPrintf( "===============================================================");
  myMexPrintf( "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
}// end itkSaveWithMetaData

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  try
    {
    itkSaveWithMetaData(nlhs, plhs, nrhs, prhs);
    }
  catch( itk::ExceptionObject & excp )
    {
    std::string msg = excp.what();
    char errBuff[NRRD_MAX_ERROR_MSG_SIZE] = { '\0' };
    snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "%s: at line %d", msg.c_str(), __LINE__);
    mexErrMsgTxt(errBuff);
    }
  catch( ... )
    {
    printf("Exception in itkSaveWithMetaData\n");
    mexErrMsgTxt("Exception in itkSaveWithMetaData");
    }
} // end mexFunction
