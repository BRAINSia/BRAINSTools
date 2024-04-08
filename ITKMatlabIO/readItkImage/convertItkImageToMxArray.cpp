/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
*/

#include <itkImage.h>
#include <itkImageFileReader.h>
#include "convertItkImageToMxArray.h"

void itkImage2MxArray(const std::string filename, mxArray* plhs[])
{
    //read image
    typedef itk::RGBPixel<double>  PixelType;
    typedef itk::Image< PixelType, 3>  ImageType;
    typedef itk::ImageFileReader<ImageType>  ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);

    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
        std::cerr << "Exception thrown while reading the image file: " << filename <<std::endl;
        std::cerr << excep << std::endl;
    }

    ImageType::Pointer image = reader->GetOutput();

    ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    if (3 != size.Dimension){
        printf("image dimension = %d, incorrect.\n", size.Dimension);
        return;
    }

    //create Mex array
    mwSize ndim =  size.Dimension+1;
    mwSize *dims = new mwSize[ndim];
    dims[0] = size[1]; //row
    dims[1] = size[0]; //column
    dims[2] = 3 ;      //RGB channel;
    dims[3] = size[2]; //slice
    plhs[0] =  mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL);
    delete[] dims;
    double* mxPointer = mxGetPr(plhs[0]);

    //convert itk 3D image into mex Array
    ImageType::IndexType itkIndex;
    for (unsigned long slice=0; slice<size[2]; ++slice)   //slice
       for (unsigned long row=0; row<size[1]; ++row)      //row
          for (unsigned long col=0; col<size[0]; ++col)  //column
    {
        itkIndex[0] = col;  itkIndex[1] = row;  itkIndex[2] = slice;
        unsigned long mxRIndex = slice*size[0]*size[1]*3 +        0          + col*size[1] + row;
        unsigned long mxGIndex = slice*size[0]*size[1]*3 + size[0]*size[1]*1 + col*size[1] + row;
        unsigned long mxBIndex = slice*size[0]*size[1]*3 + size[0]*size[1]*2 + col*size[1] + row;
        mxPointer[mxRIndex] =  image->GetPixel(itkIndex).GetRed();
        mxPointer[mxGIndex] =  image->GetPixel(itkIndex).GetGreen();
        mxPointer[mxBIndex] =  image->GetPixel(itkIndex).GetBlue();
    }

    return;
}
