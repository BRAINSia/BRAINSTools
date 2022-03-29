#include "mex.h"

#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkNrrdImageIOFactory.h"
#include <iomanip>
#include "nrrdCommon.h"
#include "myMexPrintf.h" //add by HuiXie

template <typename PixelType>
mxClassID
ITKToMType()
{
  itkGenericExceptionMacro(<< "Unhandled image type");
}

template <>
mxClassID
ITKToMType<double>()
{
  myMexPrintf("UNITS = mxDOUBLE_CLASS\n");
  return mxDOUBLE_CLASS;
}

template <>
mxClassID
ITKToMType<float>()
{
  return mxSINGLE_CLASS;
}

template <>
mxClassID
ITKToMType<int>()
{
  return mxINT32_CLASS;
}

template <>
mxClassID
ITKToMType<unsigned int>()
{
  return mxUINT32_CLASS;
}

template <>
mxClassID
ITKToMType<short>()
{
  myMexPrintf("UNITS = mxINT16_CLASS\n");
  return mxINT16_CLASS;
}

template <>
mxClassID
ITKToMType<unsigned short>()
{
  myMexPrintf("UNITS = mxUINT16_CLASS\n");
  return mxUINT16_CLASS;
}

template <>
mxClassID
ITKToMType<itk::VariableLengthVector<double>>()
{
  myMexPrintf("UNITS = mxDOUBLE_CLASS\n");
  return mxDOUBLE_CLASS;
}

template <>
mxClassID
ITKToMType<itk::VariableLengthVector<float>>()
{
  return mxSINGLE_CLASS;
}

/** In order to use Reorganize below in a template function, there
 * needs to be a valid template function for the compiler to
 * instantiate even if it isn't going to be used
 */
template <typename ImageType1, typename ImageType2>
void
Reorganize(typename ImageType1::Pointer &, typename ImageType2::Pointer &)
{
  itkGenericExceptionMacro(<< "Should never Instantiate this template function");
}

/** Go from a Vector Image to a 4D image */
template <typename TPixel>
void Reorganize(typename itk::VectorImage<TPixel, 3>::Pointer & vecImage,
                typename itk::Image<TPixel, 4>::Pointer &       scalarImage)
{
  using TVecImage = itk::VectorImage<TPixel, 3>;
  using TScalarImage = itk::Image<TPixel, 4>;
  //
  // this does a full, correct job of creating the corresponding 4D
  // image from the 3D Vector Image, even though all we care about in
  // this context is the reorganization of the voxels.
  const typename TVecImage::SizeType vecImageSize = vecImage->GetLargestPossibleRegion().GetSize();
  typename TScalarImage::SizeType    size;

  const typename TVecImage::SpacingType vecSpacing = vecImage->GetSpacing();
  typename TScalarImage::SpacingType    spacing;

  const typename TVecImage::PointType vecOrigin = vecImage->GetOrigin();
  typename TScalarImage::PointType    origin;

  const typename TVecImage::DirectionType vecDirections = vecImage->GetDirection();
  typename TScalarImage::DirectionType    directions;

  for (unsigned int i = 0; i < TScalarImage::ImageDimension; ++i)
  {
    if (i < TVecImage::ImageDimension)
    {
      size[i] = vecImageSize[i];
      spacing[i] = vecSpacing[i];
      origin[i] = vecOrigin[i];
    }
    else
    {
      size[i] = vecImage->GetNumberOfComponentsPerPixel();
      spacing[i] = 1.0;
      origin[i] = 0.0;
    }
    for (unsigned int j = 0; j < TScalarImage::ImageDimension; ++j)
    {
      if (i < TVecImage::ImageDimension && j < TVecImage::ImageDimension)
      {
        directions(i, j) = vecDirections(i, j);
      }
      else
      {
        directions(i, j) = i == j ? 1.0 : 0.0;
      }
    }
  }
  scalarImage = TScalarImage::New();
  scalarImage->SetRegions(size);
  scalarImage->SetOrigin(origin);
  scalarImage->SetSpacing(spacing);
  scalarImage->SetDirection(directions);
  scalarImage->Allocate();

  typename TVecImage::IndexType    vecIndex;
  typename TScalarImage::IndexType volIndex;

  using IndexValueType = typename TScalarImage::IndexType::IndexValueType;

  // convert from vector image to 4D volume image
  for (volIndex[3] = 0; volIndex[3] < static_cast<IndexValueType>(size[3]); ++volIndex[3])
  {
    for (volIndex[2] = 0; volIndex[2] < static_cast<IndexValueType>(size[2]); ++volIndex[2])
    {
      vecIndex[2] = volIndex[2];
      for (volIndex[1] = 0; volIndex[1] < static_cast<IndexValueType>(size[1]); ++volIndex[1])
      {
        vecIndex[1] = volIndex[1];
        for (volIndex[0] = 0; volIndex[0] < static_cast<IndexValueType>(size[0]); ++volIndex[0])
        {
          vecIndex[0] = volIndex[0];
          scalarImage->SetPixel(volIndex, vecImage->GetPixel(vecIndex)[volIndex[3]]);
        }
      }
    }
  }
} // end reorganize

/** Copy image data. In order for C++ function signature matching to
 * work properly, the function is templated over type, and it needs to
 * be repeated for every possible vector image type, since template
 * functions need to be fully specified.
 */
template <typename TImage>
void
CopyImageData(typename TImage::Pointer & im, void * target, unsigned long numPixels)
{
  using PixelType = typename TImage::PixelType;
  const PixelType * data = static_cast<PixelType *>(im->GetBufferPointer());
  std::copy(data, data + numPixels, static_cast<PixelType *>(target));
}

template <>
void
CopyImageData<itk::VectorImage<double, 3>>(itk::VectorImage<double, 3>::Pointer & im,
                                           void *                                 target,
                                           unsigned long                          numPixels)
{
  using ScalarImageType = itk::Image<double, 4>;
  ScalarImageType::Pointer newImage;
  Reorganize<double>(im, newImage);
  CopyImageData<ScalarImageType>(newImage, target, numPixels);
}

template <>
void
CopyImageData<itk::VectorImage<float, 3>>(itk::VectorImage<float, 3>::Pointer & im,
                                          void *                                target,
                                          unsigned long                         numPixels)
{
  using ScalarImageType = itk::Image<float, 4>;
  ScalarImageType::Pointer newImage;
  Reorganize<float>(im, newImage);
  CopyImageData<ScalarImageType>(newImage, target, numPixels);
}

template <>
void
CopyImageData<itk::VectorImage<short, 3>>(itk::VectorImage<short, 3>::Pointer & im,
                                          void *                                target,
                                          unsigned long                         numPixels)
{
  using ScalarImageType = itk::Image<short, 4>;
  ScalarImageType::Pointer newImage;
  Reorganize<short>(im, newImage);
  CopyImageData<ScalarImageType>(newImage, target, numPixels);
}

// build the matlab DWI image structure. This should work as well for non-DWI images
template <typename TImage>
void
BuildMatlabStruct(mxArray *& structMx, typename TImage::Pointer im)
{

  using ImageType = TImage;
  using PixelType = typename TImage::PixelType;
  const itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();

  // determine if image contains bvalue & gradients
  double      bValue(0.0);
  std::string bValueString;
  using GradientListType = std::vector<std::vector<double>>;
  GradientListType gradients;

  if (itk::ExposeMetaData<std::string>(thisDic, "DWMRI_b-value", bValueString))
  {
    std::stringstream ss(bValueString);
    ss >> bValue;

    std::stringstream gradTag;
    unsigned int      gradCount(0);
    gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << gradCount;
    gradCount++;
    std::string gradString;

    while (itk::ExposeMetaData<std::string>(thisDic, gradTag.str(), gradString))
    {
      std::stringstream   grad(gradString);
      std::vector<double> curGrad(3);
      double              val;
      for (unsigned i = 0; i < 3; ++i)
      {
        grad >> val;
        curGrad[i] = val;
      }
      gradients.push_back(curGrad);
      gradTag.str("");
      gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << gradCount;
      ++gradCount;
    }
  }


  // Setup the Fieldnames for output matlab martix
  const char ** fieldnames;
  if (gradients.size() > 0)
  {
    fieldnames = static_cast<const char **>(mxCalloc(FIELDNAME_INDEX_MAXFEILDS, sizeof(*fieldnames)));
  }
  else // no gradients
  {
    fieldnames = static_cast<const char **>(mxCalloc(FIELDNAME_INDEX_MANDATORYFEILDS, sizeof(*fieldnames)));
  }

  fieldnames[FIELDNAME_INDEX_data] = "data";
  fieldnames[FIELDNAME_INDEX_space] = "space";
  fieldnames[FIELDNAME_INDEX_spacedirections] = "spacedirections";
  fieldnames[FIELDNAME_INDEX_centerings] = "centerings";
  fieldnames[FIELDNAME_INDEX_kinds] = "kinds";
  fieldnames[FIELDNAME_INDEX_spaceunits] = "spaceunits";
  fieldnames[FIELDNAME_INDEX_spacedefinition] = "spacedefinition";
  fieldnames[FIELDNAME_INDEX_spaceorigin] = "spaceorigin";
  fieldnames[FIELDNAME_INDEX_measurementframe] = "measurementframe";

  if (gradients.size() > 0) // if there are gradients
  {
    fieldnames[FIELDNAME_INDEX_modality] = "modality";
    fieldnames[FIELDNAME_INDEX_bvalue] = "bvalue";
    fieldnames[FIELDNAME_INDEX_gradientdirections] = "gradientdirections";
    structMx = mxCreateStructMatrix(1, 1, FIELDNAME_INDEX_MAXFEILDS, fieldnames);
  }
  else
  { // no gradients
    structMx = mxCreateStructMatrix(1, 1, FIELDNAME_INDEX_MANDATORYFEILDS, fieldnames);
  }
  mxFree((void *)fieldnames);


  const typename ImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();

  const mwSize mxNrrdDim = ImageType::ImageDimension;
  // matrix size has to reflect gradient vector size if present.
  mwSize numMxDimensions = ImageType::ImageDimension;
  if (im->GetNumberOfComponentsPerPixel() > 1)
  {
    ++numMxDimensions;
  } // end if
  mwSize * sizeI = static_cast<mwSize *>(mxCalloc(numMxDimensions, sizeof(mwSize)));

  size_t numPixels(1);

  /** data **/
  {
    unsigned int axIdx;
    for (axIdx = 0; axIdx < mxNrrdDim; axIdx++)
    {
      sizeI[axIdx] = size[axIdx];
      numPixels *= size[axIdx];
    } // end for
    if (axIdx < numMxDimensions)
    {
      sizeI[axIdx] = im->GetNumberOfComponentsPerPixel();
      numPixels *= sizeI[axIdx];
    } // end if
  }   // end data section

  myMexPrintf("mxNrrdDimension = %lf\n", static_cast<double>(mxNrrdDim));
  for (unsigned i = 0; i < numMxDimensions; ++i)
  {
    myMexPrintf("sizeI = %d, i = %d, numMxDimensions = %d\n", static_cast<int>(sizeI[i]), i, numMxDimensions);
  } // end for

  // create voxel data array for matlab
  const mxClassID mtype = ITKToMType<PixelType>(); // work needed here

  mxArray * data = mxCreateNumericArray(numMxDimensions, sizeI, mtype, mxREAL);
  // copy voxels from ITK image to matlab matrix
  CopyImageData<ImageType>(im, mxGetData(data), numPixels);

  // add voxel data to matlab struct
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_data, data);

  myMexPrintf("DEFAULT after setting data\n");

  // space is handled in ITK by always converting to LPS
  mwSize    space_size = 1;
  mxArray * space = mxCreateNumericArray(1, &space_size, mxINT32_CLASS, mxREAL);
  int *     space_temp = (int *)mxGetData(space);
  *space_temp = nrrdSpaceLeftPosteriorSuperior;
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_space, space);
  myMexPrintf("DEFAULT after space\n");

  /** centerings **/
  // not tracked by ITK, so defaults
  // TODO if it's a nrrd file loaded this is in the metadata
  mwSize    centerings_size = numMxDimensions;
  mxArray * centerings = mxCreateNumericArray(1, &centerings_size, mxINT32_CLASS, mxREAL);
  int *     centerings_temp = (int *)mxGetData(centerings);
  for (unsigned int axIdx = 0; axIdx < numMxDimensions; axIdx++)
  {
    std::string       val;
    std::stringstream ss;
    ss << "NRRD_centerings[" << axIdx << "]";
    if (itk::ExposeMetaData<std::string>(thisDic, ss.str(), val))
    { // NRRD case
      int centering = airEnumVal(nrrdCenter, val.c_str());
      centerings_temp[axIdx] = centering;
      myMexPrintf("NRRD for centerings %d \n", centerings_temp[axIdx]);
    }
    else
    { // default case
      centerings_temp[axIdx] =
        axIdx < 3 ? nrrdCenterCell : nrrdCenterUnknown; // what is this?? nrrdCenterCell is this ITK default???
      myMexPrintf("DEFAULT for centerings %d \n", centerings_temp[axIdx]);
    } // end default case
  }   // end for
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_centerings, centerings);

  /** kinds **/
  // defaults for ITK
  // TODO if it's a nrrd file loaded this is in the metadata
  mxArray * kinds = mxCreateNumericArray(1, &numMxDimensions, mxINT32_CLASS, mxREAL);
  int *     kinds_temp = (int *)mxGetData(kinds);
  for (unsigned int axIdx = 0; axIdx < numMxDimensions; axIdx++)
  {
    std::string       val;
    std::stringstream ss;

    ss << "NRRD_kinds[" << axIdx << "]";
    if (itk::ExposeMetaData<std::string>(thisDic, ss.str(), val))
    { // nrrd case
      int kind = airEnumVal(nrrdKind, val.c_str());
      kinds_temp[axIdx] = kind;
      myMexPrintf("NRRD kind %d\n", kinds_temp[axIdx]);
    }
    else if (gradients.size() > 0) // i.e. We are a DWI image
    {
      kinds_temp[axIdx] = axIdx < numMxDimensions - 1 ? nrrdKindSpace : nrrdKindList;
      myMexPrintf("DEFAULT DWI kind %d\n", kinds_temp[axIdx]);
    }
    else
    { // default case//debug
      kinds_temp[axIdx] = nrrdKindSpace;
      myMexPrintf("DEFAULT kind %d\n", kinds_temp[axIdx]);
    } // end default case
  }   // end for
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_kinds, kinds);


  /** spacedirections **/
  // Not sure what is actually done with this in matlab. ITK has one
  // concept of orientation, NRRD had several, and just blindly
  // assigning the directions from a NRRD file would mean that any
  // matlab program that depended on orieintation would have to take
  // the nrrd space variable into account.
  // For now, just have to assume that the ITK directions are OK as is.
  // Working for NRRD so nothing to do here.
  mxArray * spacedirections = mxCreateNumericMatrix(mxNrrdDim, mxNrrdDim, mxDOUBLE_CLASS, mxREAL);
  double *  spacedirections_temp = (double *)mxGetData(spacedirections);
  const typename ImageType::DirectionType directions = im->GetDirection();
  typename ImageType::SpacingType         itkSpacing = im->GetSpacing();

  for (unsigned int axIdx = 0, count = 0; axIdx < ImageType::ImageDimension; ++axIdx)
  {
    for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx)
    {
      spacedirections_temp[count] = directions[sdIdx][axIdx] * itkSpacing[axIdx];
      count++;
    } // end for
  }   // end for
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_spacedirections, spacedirections);
  myMexPrintf("DEFAULT space directions\n");


  /** spaceunits **/
  // Never enters NRRD case, nothing to do for NRRD
  std::string theSpaceUnitsString;
  if (itk::ExposeMetaData<std::string>(thisDic, "NRRD_spaceunits", theSpaceUnitsString))
  { // NRRD case
    mxArray * spaceunits = mxCreateString(theSpaceUnitsString.c_str());
    mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_spaceunits, spaceunits);
    myMexPrintf("NRRD spaceunits \"%s\"\n", theSpaceUnitsString.c_str());
  }
  else
  { // DEFAULT case
    mxArray * spaceunits = mxCreateCellArray(1, &mxNrrdDim);
    for (unsigned sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++)
    {
      mxArray * spaceunits_temp = mxCreateString("mm");
      mxSetCell(spaceunits, sdIdx, spaceunits_temp);
    } // end for
    mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_spaceunits, spaceunits);
    myMexPrintf("DEFAULT spaceunits \"mm\"\n");
  } // end default case

  /** spacedefinition **/
  std::string theSpaceString;
  if (itk::ExposeMetaData<std::string>(thisDic, "NRRD_space", theSpaceString))
  { // NRRD case
    myMexPrintf("NRRD space definintion \"%s\"\n", theSpaceString.c_str());
  }
  else
  { // Default case
    theSpaceString = airEnumStr(nrrdSpace, nrrdSpaceLeftPosteriorSuperior);
    myMexPrintf("DEFAULT space definintion \"%s\"\n", theSpaceString.c_str());
  } // end default case
  mxArray * spacedefinition = mxCreateString(theSpaceString.c_str());
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_spacedefinition, spacedefinition);

  /** spaceorigin **/
  mxArray *                     spaceorigin = mxCreateNumericArray(1, &mxNrrdDim, mxDOUBLE_CLASS, mxREAL);
  double *                      spaceorigin_temp = (double *)mxGetData(spaceorigin);
  typename ImageType::PointType origin = im->GetOrigin(); // from ITK::Image
  for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++)
  {
    spaceorigin_temp[sdIdx] = origin[sdIdx];
  } // end for
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_spaceorigin, spaceorigin);
  myMexPrintf("DEFAULT space origin\n");

  /** measurementframe -- only significant in NRRD files **/
  // http://teem.sourceforge.net/nrrd/format.html#measurementframe
  std::vector<std::vector<double>> measurementFrame;
  std::string                      measurementFrameFieldName = "NRRD_";
  measurementFrameFieldName += airEnumStr(nrrdField, nrrdField_measurement_frame);
  mxArray * measurementframe = mxCreateNumericMatrix(mxNrrdDim, mxNrrdDim, mxDOUBLE_CLASS, mxREAL);
  double *  measurementframe_temp = static_cast<double *>(mxGetData(measurementframe));
  // measurementFrameFieldName = "NRRD_measurement frame"
  if (itk::ExposeMetaData<std::vector<std::vector<double>>>(thisDic, measurementFrameFieldName, measurementFrame))
  { // NRRD special case
    for (unsigned int axIdx = 0, count = 0; axIdx < ImageType::ImageDimension; ++axIdx)
    {
      for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx, ++count)
      {
        measurementframe_temp[count] = measurementFrame[axIdx][sdIdx];
        myMexPrintf("NRRD measurement frame: %lf\n", measurementframe_temp[count]);
      } // end for
    }   // end for
  }
  else
  { // Default case
    for (unsigned int axIdx = 0, count = 0; axIdx < ImageType::ImageDimension; ++axIdx)
    {
      for (unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx, ++count)
      {
        measurementframe_temp[count] = 0.0;
        myMexPrintf("DEFAULT measurement frame: %lf\n", measurementframe_temp[count]);
      } // end for
    }   // end for
  }     // end else
  mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_measurementframe, measurementframe);


  /** modality, bvalue, & gradientdirections **/
  if (gradients.size() > 0)
  { // only for images with gradient
    /** modality **/
    mxArray * modality = mxCreateString("DWMRI");
    mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_modality, modality);
    myMexPrintf("DEFAULT set modality\n");

    /* use ten to parse the key/value pairs */
    /** bvalue **/
    constexpr mwSize bvalue_size = 1;
    mxArray * const  bvalue = mxCreateNumericArray(1, &bvalue_size, mxDOUBLE_CLASS, mxREAL);
    double * const   bvalue_temp = (double *)mxGetData(bvalue);
    *bvalue_temp = bValue;
    mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_bvalue, bvalue);
    myMexPrintf("DEFAULT BValue\n");

    /** gradientdirections **/
    /*First find the index that contains the graident directions! */
    mxArray * gradientdirections = mxCreateNumericMatrix(gradients.size(), mxNrrdDim, mxDOUBLE_CLASS, mxREAL);
    double *  gradientdirections_temp = (double *)mxGetData(gradientdirections);
    for (unsigned int dwiIdx = 0; dwiIdx < gradients.size(); ++dwiIdx)
    {
      std::stringstream           ss;
      const std::vector<double> & curGrad = gradients[dwiIdx];
      ss << dwiIdx << " " << curGrad[0] << " " << curGrad[1] << " " << curGrad[2] << std::endl;
      // myMexPrintf( ss.str().c_str() ); //for printing grad to matlab for debugging
      gradientdirections_temp[dwiIdx] = curGrad[0];
      gradientdirections_temp[dwiIdx + gradients.size()] = curGrad[1];
      gradientdirections_temp[dwiIdx + gradients.size() * 2] = curGrad[2];
    } // end for
    mxSetFieldByNumber(structMx, 0, FIELDNAME_INDEX_gradientdirections, gradientdirections);
    myMexPrintf("DEFAULT gradients\n");
  } // end if
} // end BuildMatlabStruct

/** do the actual load */
template <typename TImage>
void
LoadDWIImage(const std::string & filename, mxArray *& structMx)
{
  using ImageType = TImage;
  using ReaderType = itk::ImageFileReader<ImageType>;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::string msg = "Error: itk::ImageFileReader can't read " + filename + " : " + e.what();
    mexErrMsgTxt(msg.c_str());
    return; // add by Hui Xie on Nov 3th, 2016
  }
  typename ImageType::Pointer im = reader->GetOutput();
  BuildMatlabStruct<ImageType>(structMx, im);
}

void
itkLoadWithMetadata(mxArray * plhs[], const mxArray * prhs[])
{
  // Get filename
  std::string filename;
  {
    const int filenameLen = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
    char *    _filename = new char[filenameLen];
    mxGetString(prhs[0], _filename, filenameLen);
    filename = _filename;
    delete[] _filename;
  }

  // try to read image file and get ImageIO pointer
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
  if (imageIO.IsNotNull())
  {
    if (imageIO->CanReadFile(filename.c_str()))
    {
      imageIO->SetFileName(filename);
      imageIO->ReadImageInformation();
    }
    else
    {
      std::string msg = "Error: itk::ImageIOBase::Pointer can't read " + filename;
      mexErrMsgTxt(msg.c_str());
    }
  }
  else
  {
    std::cerr << "Error Possible: the input filename may be not exist or file format is incorrect." << std::endl;
    std::string msg = "Error: itk::ImageIOFactory::CreateImageIO failed for " + filename;
    mexErrMsgTxt(msg.c_str());
  }

  //=========For now support scalar images of 2 or 3 dimensions========================
  // according to different pixel type and dimension, instantiate different LoadDWIImage
  plhs[0] = NULL;
  itk::ImageIOBase::IOPixelType pixtype = imageIO->GetPixelType();
  unsigned                      imageDimension = imageIO->GetNumberOfDimensions();
  switch (pixtype)
  {
    case itk::ImageIOBase::SCALAR:
      switch (imageDimension)
      {
        case 2:
          LoadDWIImage<itk::Image<double, 2>>(filename, plhs[0]);
          break;
        case 3:
          if (imageIO->GetComponentType() == itk::ImageIOBase::SHORT ||
              imageIO->GetComponentType() == itk::ImageIOBase::USHORT)
          {
            LoadDWIImage<itk::Image<unsigned short, 3>>(filename, plhs[0]);
          }
          else if (imageIO->GetComponentType() == itk::ImageIOBase::DOUBLE) // default
          {
            LoadDWIImage<itk::Image<double, 3>>(filename, plhs[0]);
          }
          else // ITK UNIT DEFAULT DOUBLE 3D
          {
            LoadDWIImage<itk::Image<double, 3>>(filename, plhs[0]);
          }
          break;
        default:
        {
          std::stringstream ss;
          ss << filename << " has " << imageDimension << " dimensions. Only 2 or 3 supported";
          mexErrMsgTxt(ss.str().c_str());
        }
      }
      break;
    case itk::ImageIOBase::DIFFUSIONTENSOR3D:
    case itk::ImageIOBase::VECTOR:
      if (imageDimension == 3)
      {
        switch (imageIO->GetComponentType())
        {
          case itk::ImageIOBase::DOUBLE:
            LoadDWIImage<itk::VectorImage<double, 3>>(filename, plhs[0]);
            break;
          case itk::ImageIOBase::FLOAT:
            LoadDWIImage<itk::VectorImage<float, 3>>(filename, plhs[0]);
            break;
          case itk::ImageIOBase::SHORT:
            LoadDWIImage<itk::VectorImage<short, 3>>(filename, plhs[0]);
            break;
          default:
          {
            std::string errmsg =
              "Unsupported Pixel Type in file " + filename + " : " + imageIO->GetPixelTypeAsString(pixtype);
            mexErrMsgTxt(errmsg.c_str());
          }
        }
      }
      else
      {
        std::stringstream ss;
        ss << filename << " has " << imageDimension << " dimensions. Only 3 supported";
        mexErrMsgTxt(ss.str().c_str());
      }
      break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    case itk::ImageIOBase::RGB:
    case itk::ImageIOBase::RGBA:
    case itk::ImageIOBase::OFFSET:
    case itk::ImageIOBase::POINT:
    case itk::ImageIOBase::COVARIANTVECTOR:
    case itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR:
    case itk::ImageIOBase::COMPLEX:
    case itk::ImageIOBase::FIXEDARRAY:
    case itk::ImageIOBase::MATRIX:
    {
      std::string errmsg =
        "Unsupported Pixel Type in file " + filename + " : " + imageIO->GetPixelTypeAsString(pixtype);
      mexErrMsgTxt(errmsg.c_str());
    }
  }
}

void
mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  // check input parameters
  if (1 != nlhs || 1 != nrhs || NULL == mxGetData(prhs[0]) || !mxIsChar(prhs[0]))
  {
    std::cerr << "Parameters Error in using itkLoadWithMetadata" << std::endl;
    std::cerr << "Usage: outputStruct = itkLoadWithMetadata('nrrdImageFilename')" << std::endl;
    return;
  }

  try
  {
    itkLoadWithMetadata(plhs, prhs);
  }
  catch (std::exception & e)
  {
    mexErrMsgTxt(e.what());
  }
}
