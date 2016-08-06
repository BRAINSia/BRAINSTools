#include "mex.h"

#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkNrrdImageIOFactory.h"
#include "itkImageIOFactoryRegisterManager.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"
#include <iomanip>
#include <sstream>
#include <list>
#include <vector>
#include <algorithm>
#include "nrrdCommon.h"
#include <typeinfo>
#include "Debug.h"

/** ITKToMType -- use template function to map ITK/native scalar types
 * to Matlab scalar types.
 */
template <typename PixelType>
mxClassID ITKToMType()
{
    itkGenericExceptionMacro(<< "Unhandled image type");
}
template <> mxClassID ITKToMType<double>() { return mxDOUBLE_CLASS; }
template <> mxClassID ITKToMType<float>() { return mxSINGLE_CLASS; }
template <> mxClassID ITKToMType<int>() { return mxINT32_CLASS; }
template <> mxClassID ITKToMType<unsigned int>() { return mxUINT32_CLASS; }
template <> mxClassID ITKToMType<short>() { return mxINT16_CLASS; }
template <> mxClassID ITKToMType<unsigned short>() { return mxUINT16_CLASS; }
template<> mxClassID ITKToMType< itk::VariableLengthVector<double> >() { return mxDOUBLE_CLASS; }
template <>mxClassID ITKToMType< itk::VariableLengthVector<float> >() { return mxSINGLE_CLASS; }

/** In order to use Reorganize below in a template function, there
 * needs to be a valid template function for the compiler to
 * instantiate even if it isn't going to be used
 */
template <typename ImageType1,typename ImageType2>
void
Reorganize(typename ImageType1::Pointer &,
           typename ImageType2::Pointer &)
{
    itkGenericExceptionMacro(<< "Should never Instantiate this template function");
}

/** Go from a Vector Image to a 4D image */
template <typename TPixel>
void
Reorganize(typename itk::VectorImage<TPixel,3>::Pointer &vecImage,
           typename itk::Image<TPixel,4>::Pointer &scalarImage)
{
    typedef itk::VectorImage<TPixel,3> TVecImage;
    typedef itk::Image<TPixel,4> TScalarImage;
    //
    // this does a full, correct job of creating the corresponding 4D
    // image from the 3D Vector Image, even though all we care about in
    // this context is the reorganization of the voxels.
    const typename TVecImage::SizeType vecImageSize =
    vecImage->GetLargestPossibleRegion().GetSize();
    typename TScalarImage::SizeType size;

    const typename TVecImage::SpacingType  vecSpacing =
    vecImage->GetSpacing();
    typename TScalarImage::SpacingType spacing;

    const typename TVecImage::PointType vecOrigin =
    vecImage->GetOrigin();
    typename TScalarImage::PointType    origin;

    const typename TVecImage::DirectionType vecDirections =
    vecImage->GetDirection();
    typename  TScalarImage::DirectionType    directions;

    for(unsigned int i = 0; i < TScalarImage::ImageDimension; ++i)
    {
        if(i < TVecImage::ImageDimension)
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
        for(unsigned int j = 0; j < TScalarImage::ImageDimension; ++j)
        {
            if(i < TVecImage::ImageDimension && j < TVecImage::ImageDimension)
            {
                directions(i,j) = vecDirections(i,j);
            }
            else
            {
                directions(i,j) = i == j ? 1.0 : 0.0;
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

    typedef typename TScalarImage::IndexType::IndexValueType IndexValueType;

    // convert from vector image to 4D volume image
    for( volIndex[3] = 0; volIndex[3] < static_cast<IndexValueType>(size[3]); ++volIndex[3] )
    {
        for( volIndex[2] = 0; volIndex[2] < static_cast<IndexValueType>(size[2]); ++volIndex[2] )
        {
            vecIndex[2] = volIndex[2];
            for( volIndex[1] = 0; volIndex[1] < static_cast<IndexValueType>(size[1]); ++volIndex[1] )
            {
                vecIndex[1] = volIndex[1];
                for( volIndex[0] = 0; volIndex[0] < static_cast<IndexValueType>(size[0]); ++volIndex[0] )
                {
                    vecIndex[0] = volIndex[0];
                    scalarImage->SetPixel(volIndex, vecImage->GetPixel(vecIndex)[volIndex[3]]);
                }
            }
        }
    }
}/** END reorganize Vector Image to a 4D image */

/** Copy image data. In order for C++ function signature matching to
 * work properly, the function is templated over type, and it needs to
 * be repeated for every possible vector image type, since template
 * functions need to be fully specified.
 */
template<typename TImage>
void CopyImageData(typename TImage::Pointer &im, void *target, unsigned long numPixels)
{
    typedef typename TImage::PixelType PixelType;
    const PixelType *data = static_cast<PixelType *>(im->GetBufferPointer());
    std::copy(data,data + numPixels,static_cast<PixelType *>(target));
}

template<>
void CopyImageData<itk::VectorImage<double,3> >( itk::VectorImage<double,3>::Pointer &im,
                                                void *target, unsigned long numPixels)
{
    typedef itk::Image<double,4> ScalarImageType;
    ScalarImageType::Pointer newImage;
    Reorganize<double>(im,newImage);
    CopyImageData<ScalarImageType>(newImage,target,numPixels);
}

template<>
void CopyImageData<itk::VectorImage<float,3> >( itk::VectorImage<float,3>::Pointer &im,
                                               void *target, unsigned long numPixels)
{
    typedef itk::Image<float,4> ScalarImageType;
    ScalarImageType::Pointer newImage;
    Reorganize<float>(im,newImage);
    CopyImageData<ScalarImageType>(newImage,target,numPixels);
}

template<>
void CopyImageData<itk::VectorImage<short,3> >( itk::VectorImage<short,3>::Pointer &im,
                                               void *target, unsigned long numPixels)
{
    typedef itk::Image<short,4> ScalarImageType;
    ScalarImageType::Pointer newImage;
    Reorganize<short>(im,newImage);
    CopyImageData<ScalarImageType>(newImage,target,numPixels);
}/** END Copy image data.*/

/** build the matlab DWI image structure.
 *  This should work as well for non-DWI images.
 */
template <typename TImage>
void
BuildMatlabStruct(mxArray *& structMx, typename TImage::Pointer im, const std::string & filename)
{
    typedef TImage ImageType;
    typedef typename TImage::PixelType PixelType;

    DEBUG_PRINT(std::string("filename = ") + filename, nl);

    DEBUG_PRINTDICT(im.GetPointer() );
    const itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();

    DEBUG_PRINT("after image read", nl);


    // determine if image contains bvalue & gradients.
    double      bValue(0.0);
    std::string bValueString;
    typedef std::vector<std::vector<double> > GradientListType;
    GradientListType gradients;

    if( itk::ExposeMetaData<std::string>(thisDic, "DWMRI_b-value", bValueString) == true )
    {
        std::stringstream ss(bValueString);
        ss >> bValue;

        std::stringstream gradTag;
        unsigned int      gradCount(0);;
        gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << gradCount;
        gradCount++;
        std::string gradString;
        for( ; itk::ExposeMetaData<std::string>(thisDic, gradTag.str(), gradString) != false; ++gradCount )
        {
            DEBUG_PRINT(gradTag.str().c_str() ); DEBUG_PRINT("=");
            DEBUG_PRINT(gradString, nl);
            std::stringstream   grad(gradString);
            std::vector<double> curGrad(3);
            double              val;
            for( unsigned i = 0; i < 3; ++i )
            {
                grad >> val;
                curGrad[i] = val;
            }
            gradients.push_back(curGrad);
            gradTag.str("");
            gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << gradCount;
        }
    }//end check for DWMRI_b-value

    /** Setup the Fieldnames **/
    const char * *fieldnames;  /* matlab struct field names*/
    if( gradients.size() > 0 )
    {
        fieldnames = static_cast<const char * *>(mxCalloc(FIELDNAME_INDEX_MAXFEILDS, sizeof(*fieldnames) ) );
    }
    else
    {
        fieldnames = static_cast<const char * *>(mxCalloc(FIELDNAME_INDEX_MANDATORYFEILDS, sizeof(*fieldnames) ) );
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
    if( gradients.size() > 0 )
    {
        fieldnames[FIELDNAME_INDEX_modality] = "modality";
        fieldnames[FIELDNAME_INDEX_bvalue] = "bvalue";
        fieldnames[FIELDNAME_INDEX_gradientdirections] = "gradientdirections";
        structMx = mxCreateStructMatrix( 1, 1, FIELDNAME_INDEX_MAXFEILDS, fieldnames );
    }
    else
    {
        structMx = mxCreateStructMatrix( 1, 1, FIELDNAME_INDEX_MANDATORYFEILDS, fieldnames );
    }

    mxFree( (void *)fieldnames);

    DEBUG_PRINT("Created Field Names", nl);

    const typename ImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();

    const mwSize mxNrrdDim = ImageType::ImageDimension;
    // matrix size has to reflect gradient vector size if present.
    mwSize numMxDimensions = ImageType::ImageDimension;
    if( im->GetNumberOfComponentsPerPixel() > 1 )
    {
        ++numMxDimensions;
    }
    mwSize *sizeI = static_cast<mwSize *>(mxCalloc(numMxDimensions, sizeof(mwSize) ) );

    size_t numPixels(1);
    /**DATA**/
    unsigned int axIdx;
    for( axIdx = 0; axIdx < mxNrrdDim; axIdx++ )
    {
        sizeI[axIdx] = size[axIdx];
        numPixels *= size[axIdx];
    }
    if( axIdx < numMxDimensions )
    {
        sizeI[axIdx] = im->GetNumberOfComponentsPerPixel();
        numPixels *= sizeI[axIdx];
    }

    DEBUG_PRINT("dim = "); DEBUG_PRINT(mxNrrdDim, nl);
    for( unsigned i = 0; i < numMxDimensions; ++i )
    {
        DEBUG_PRINT(sizeI[i], i, numMxDimensions);
    }

    //
    // create voxel data array for matlab
    const mxClassID mtype = ITKToMType<PixelType>();

    mxArray *data = mxCreateNumericArray( numMxDimensions, sizeI, mtype, mxREAL );
    // copy voxels from ITK image to matlab matrix
    CopyImageData<ImageType>(im,mxGetData(data),numPixels);

    // add voxel data to matlab struct
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_data, data );

    DEBUG_PRINT("after setting data", nl);

    // space is handled in ITK by always converting to LPS
    mwSize   space_size = 1;
    mxArray *space = mxCreateNumericArray( 1, &space_size, mxINT32_CLASS, mxREAL );
    int *    space_temp = (int *)mxGetData(space);
    *space_temp = nrrdSpaceLeftPosteriorSuperior;
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_space, space );
    DEBUG_PRINT("after space", nl);

    /** centerings **/
    // not tracked by ITK, so defaults
    // TODO if it's a nrrd file loaded this is in the metadata
    //TRYING TO GET centerings FROM METADATA in nrrd file

    // TODO if it's a nrrd file loaded this is in the metadata
    mwSize   centerings_size = numMxDimensions;
    mxArray *centerings = mxCreateNumericArray( 1, &centerings_size, mxINT32_CLASS, mxREAL );
    int *    centerings_temp = (int *)mxGetData(centerings);
    if(DEFUALT)//NOT NRRD
    {
        for( unsigned int axIdx = 0; axIdx < centerings_size; ++axIdx )
        {
            if( axIdx < 3 )
            {
                centerings_temp[axIdx] =  nrrdCenterCell;
            }
            else  // The 4th dimension is the vector dimension for 4D DWI data.
            {
                centerings_temp[axIdx] =  nrrdCenterUnknown;
            }
        }
    }// end DEFAULT
    else{// IS A NRRD
        for( unsigned int axIdx = 0; axIdx < centerings_size; ++axIdx )
        {
            //centerings_temp[axIdx] = im->axis[axIdx].centering;// THIS MAY NOT WORK!!!
            centerings_temp[axIdx] = mxGetData(axis)[axIdx].centering;
        }
    }// end NRRD case
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_centerings, centerings );
    DEBUG_PRINT("after centering", nl);

    /** kinds **/
    // defaults for ITK
    // TODO if it's a nrrd file loaded this is in the metadata
    mxArray *kinds = mxCreateNumericArray( 1, &numMxDimensions, mxINT32_CLASS, mxREAL );
    int *    kinds_temp = (int *)mxGetData(kinds);
    if(DEFAULT)//NOT A NRRD
    {
        for( unsigned int axIdx = 0; axIdx < numMxDimensions; axIdx++ )
        {
            std::string       val;
            std::stringstream ss;
            ss << "NRRD_kinds[" << axIdx << "]";
            if( itk::ExposeMetaData<std::string>(thisDic, ss.str(), val) != false )
            {
                int kind = airEnumVal(nrrdKind, val.c_str() );
                kinds_temp[axIdx] = kind;
            }
            else
            {
                kinds_temp[axIdx] = axIdx < numMxDimensions - 1 ?
                nrrdKindSpace : nrrdKindList;
            }
        }
    }else{//IF IT IS A NRRD!!!
        for( unsigned int axIdx = 0; axIdx < nrrd->dim; axIdx++ )
        {
            //kinds_temp[axIdx] = im->axis[axIdx].kind;// THIS MAY NOT WORK
            kinds_temp[axIdx] = mxGetData(axis)[axIdx].kind;
        }
    }//end NRRD case
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_kinds, kinds );
    DEBUG_PRINT("after kind", nl);

    /** spacedirections **/
    // Not sure what is actually done with this in matlab. ITK has one
    // concept of orientation, NRRD had several, and just blindly
    // assigning the directions from a NRRD file would mean that any
    // matlab program that depended on orieintation would have to take
    // the nrrd space variable into account.
    // For now, just have to assume that the ITK directions are OK as is.
    // MUST BE OKAY SINCE IT IS NOT RAISING ANY FLAGS WHEN COMPARING WITH nrrdLoadWithMetaData
    mxArray *spacedirections = mxCreateNumericMatrix(mxNrrdDim,
                                                     mxNrrdDim, mxDOUBLE_CLASS, mxREAL );
    double *spacedirections_temp = (double *)mxGetData(spacedirections);

    const typename ImageType::DirectionType directions = im->GetDirection();

    typename ImageType::SpacingType itkSpacing = im->GetSpacing();
    for( unsigned int axIdx = 0, count = 0; axIdx < ImageType::ImageDimension; ++axIdx )
    {
        for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx )
        {
            spacedirections_temp[count] = directions[sdIdx][axIdx] * itkSpacing[axIdx];
            count++;
        }
    }
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_spacedirections, spacedirections );
    DEBUG_PRINT("after space directions", nl);

    /** spaceunits **/
    // TODO -- could read from metaDataDictionary, just default for now
    // GETTING spaceunits from mxGetData
    mxArray *spaceunits = mxCreateCellArray( 1, &mxNrrdDim );
    for( unsigned sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++ )
    {
        mxArray *spaceunits_temp = mxCreateString( mxGetData(spaceunits)[sdIdx] ); //originally simply "mm"
        //mxArray *spaceunits_temp = mxCreateString( im->spaceunits[sdIdx] );
        mxSetCell(spaceunits, sdIdx, spaceunits_temp);
    }
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_spaceunits, spaceunits );
    DEBUG_PRINT("after space units", nl);

    // TODO What the hell is this? does it matter?
    /** spacedefinition **/
    //NOT SURE WHAT THIS IS EITHER
    std::string theSpaceString;
    if( itk::ExposeMetaData<std::string>(thisDic, "NRRD_space", theSpaceString) == false )
    {
        theSpaceString = airEnumStr(nrrdSpace, nrrdSpaceLeftPosteriorSuperior);
    }
    mxArray *spacedefinition = mxCreateString( theSpaceString.c_str() );
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_spacedefinition, spacedefinition );
    DEBUG_PRINT("after space definition", nl);

    /** spaceorigin **/
    mxArray *spaceorigin = mxCreateNumericArray( 1, &mxNrrdDim, mxDOUBLE_CLASS, mxREAL );
    double * spaceorigin_temp = (double *)mxGetData(spaceorigin);
    typename ImageType::PointType origin = im->GetOrigin();
    for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; sdIdx++ )
    {
        spaceorigin_temp[sdIdx] = origin[sdIdx];
    }
    mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_spaceorigin, spaceorigin );
    DEBUG_PRINT("after space origin", nl);

    /** measurementframe -- only significant in NRRD files **/
    std::vector<std::vector<double> > measurementFrame;
    std::string                       measurementFrameFieldName = "NRRD_";
    measurementFrameFieldName += airEnumStr(nrrdField, nrrdField_measurement_frame);
    if( itk::ExposeMetaData<std::vector<std::vector<double> > >(thisDic,
                                                                measurementFrameFieldName,
                                                                measurementFrame) != false )
    {
        mxArray *measurementframe = mxCreateNumericMatrix(mxNrrdDim, mxNrrdDim, mxDOUBLE_CLASS, mxREAL );
        double * measurementframe_temp = static_cast<double *>(mxGetData(measurementframe) );
        for( unsigned int axIdx = 0, count = 0; axIdx < ImageType::ImageDimension; ++axIdx )
        {
            for( unsigned int sdIdx = 0; sdIdx < ImageType::ImageDimension; ++sdIdx, ++count )
            {
                measurementframe_temp[count] = measurementFrame[axIdx][sdIdx];
            }
        }
        mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_measurementframe, measurementframe );
    }
    DEBUG_PRINT("after measurement frame", nl);

    /** modality, bvalue, & gradientdirections **/
    if( gradients.size() > 0 )
    {
        /** modality **/
        mxArray * modality = mxCreateString( "DWMRI" );
        mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_modality, modality );

        DEBUG_PRINT("after set modality", nl);

        /* use ten to parse the key/value pairs */
        /** bvalue **/
        const mwSize    bvalue_size = 1;
        mxArray * const bvalue = mxCreateNumericArray( 1, &bvalue_size, mxDOUBLE_CLASS, mxREAL );
        double * const  bvalue_temp = (double *)mxGetData(bvalue);
        *bvalue_temp = bValue;
        mxSetFieldByNumber( structMx, 0, FIELDNAME_INDEX_bvalue, bvalue );
        DEBUG_PRINT("after BValue", nl);

        /** gradientdirections **/
        /*First find the index that contains the graident directions! */
        mxArray *gradientdirections = mxCreateNumericMatrix( gradients.size(),
                                                            mxNrrdDim,
                                                            mxDOUBLE_CLASS,
                                                            mxREAL );
        double *gradientdirections_temp = (double *)mxGetData(gradientdirections);
        for( unsigned int dwiIdx = 0; dwiIdx < gradients.size(); ++dwiIdx )
        {
            std::stringstream           ss;
            const std::vector<double> & curGrad = gradients[dwiIdx];
            ss << dwiIdx << " " << curGrad[0] << " " << curGrad[1] << " "
            << curGrad[2] << std::endl;
            DEBUG_PRINT(ss.str() );
            gradientdirections_temp[dwiIdx] = curGrad[0];
            gradientdirections_temp[dwiIdx + gradients.size()] = curGrad[1];
            gradientdirections_temp[dwiIdx + gradients.size() * 2] = curGrad[2];
        }
        mxSetFieldByNumber( structMx, 0,
                           FIELDNAME_INDEX_gradientdirections, gradientdirections );
        DEBUG_PRINT("after gradients", nl);
    }
}//DONE
//Time for actual loading
/** do the actual load */
template <typename TImage>
void
LoadDWIImage(const std::string &filename,
             mxArray *& structMx)
{
    typedef TImage ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::string msg("itkLoadWithMetaData: can't read ");
        msg += filename;
        msg += " ";
        msg += excp.what();
        mexErrMsgTxt(msg.c_str() );
    }
    typename ImageType::Pointer im = reader->GetOutput();
    BuildMatlabStruct<ImageType>(structMx, im, filename);
}
void itkLoadWithMetadata(int, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // NOTE: Useful for "printf" debugging
    // SHOULD USE TO FIGURE OUT WHERE WE ARE IN MATLAB
    // snprintf(errBuff, NRRD_MAX_ERROR_MSG_SIZE, "HERE %d: %s ", __LINE__,
    // __FILE__);
    // mexWarnMsgTxt(errBuff);

    /* Meta data stuff */
    const mxArray * const filenameMx = prhs[0];

    if( !(1 == nrhs && mxIsChar(filenameMx) ) )
    {
        mexErrMsgTxt("itkLoadWithMetadata: requires one string argument (the name of the file)");
    }

    std::string filename;
    {
        const int filenameLen = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
        char *    _filename = new char[filenameLen];
        mxGetString(prhs[0], _filename, filenameLen);
        filename = _filename;
        delete [] _filename;
    }
    //
    // have to find out what type of file this is.
    itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
    if( imageIO.IsNotNull() )
    {
        try
        {
            imageIO->SetFileName(filename);
            imageIO->ReadImageInformation();
        }
        catch( itk::ExceptionObject & excp )
        {
            std::string msg("itkLoadWithMetaData: can't read ");
            msg += filename;
            msg += " ";
            msg += excp.what();
            mexErrMsgTxt(msg.c_str() );
        }
    }
    else  // can't find proper reader for file
    {
        std::string msg = "Can't file ITK reader for ";
        msg += filename;
        mexErrMsgTxt(msg.c_str() );
    }
    mxArray * structMx = NULL;
    //
    // For now support scalar images of 2 or 3 dimensions.  Adding more
    // isn't a problem, but it complicates how we build the matlab structure.
    itk::ImageIOBase::IOPixelType pixtype = imageIO->GetPixelType();
    unsigned                      imageDimension = imageIO->GetNumberOfDimensions();

    switch( pixtype )
    {
        case itk::ImageIOBase::SCALAR:
            switch( imageDimension )
        {
            case 2:
                LoadDWIImage<itk::Image<double,2> >(filename,structMx);
                break;
            case 3:
                LoadDWIImage<itk::Image<double,3> >(filename,structMx);
                break;
            default:
            {
                std::stringstream ss;
                ss << filename << " has " << imageDimension
                << " dimensions. Only 2 or 3 supported";
                mexErrMsgTxt(ss.str().c_str() );
            }
        }
            break;
        case itk::ImageIOBase::DIFFUSIONTENSOR3D:
        case itk::ImageIOBase::VECTOR:
            if( imageDimension == 3 )
            {
                switch(imageIO->GetComponentType())
                {
                    case itk::ImageIOBase::DOUBLE:
                        LoadDWIImage<itk::VectorImage<double, 3> >(filename,structMx);
                        break;
                    case itk::ImageIOBase::FLOAT:
                        LoadDWIImage<itk::VectorImage<float, 3> >(filename,structMx);
                        break;
                    case itk::ImageIOBase::SHORT:
                        LoadDWIImage<itk::VectorImage<short, 3> >(filename,structMx);
                        break;
                    default:
                    {
                        std::string errmsg = "Unsupported Pixel Type in file ";
                        errmsg += filename;
                        errmsg += " : ";
                        errmsg += imageIO->GetPixelTypeAsString(pixtype);
                        mexErrMsgTxt(errmsg.c_str() );
                    }
                }
            }
            else
            {
                std::stringstream ss;
                ss << filename << " has " << imageDimension
                << " dimensions. Only 3 supported";
                mexErrMsgTxt(ss.str().c_str() );
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
            std::string errmsg = "Unsupported Pixel Type in file ";
            errmsg += filename;
            errmsg += " : ";
            errmsg += imageIO->GetPixelTypeAsString(pixtype);
            mexErrMsgTxt(errmsg.c_str() );
        }
    }
    plhs[0] = structMx;
    return;
}

void mexFunction(int nhls, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try
    {
        itkLoadWithMetadata(nhls, plhs, nrhs, prhs);
    }
    catch( ... )
    {
        printf("Exception in itkLoadWithMetaData\n");
        mexErrMsgTxt("Exception in itkLoadWithMetaData");
    }
}
