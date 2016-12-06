/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
*/

#include <itkImage.h>
#include <itkImageFileReader.h>
#include "itk3DImageToMxArray.h"

void itk3DImage2MxArray(const std::string filename, mxArray* plhs[])
{
    //read image
    typedef double  PixelType;
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
        std::cerr << "Usage: outputImageArray = itkRead3DImage('Scalar3DImageFilename') " <<std::endl;
        std::cerr << "This program only supports ITK scalar 3D image." << std::endl;
        std::cerr << "Exception thrown while reading the image file: " << filename <<std::endl;
        std::cerr << excep.GetDescription() << std::endl;
        return;
    }

    ImageType::Pointer image = reader->GetOutput();

    ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    if (3 != size.Dimension){
        std::cerr << "image dimension = %d, incorrect." <<size.Dimension << std::endl;
        return;
    }

    //create Mex array
    mwSize ndim =  size.Dimension;
    mwSize *dims = new mwSize[ndim];
    dims[0] = size[1]; //row
    dims[1] = size[0]; //column
    dims[2] = size[2]; //slice
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
        unsigned long mxIndex = slice*size[0]*size[1] + col*size[1] + row;
        mxPointer[mxIndex] =  image->GetPixel(itkIndex);
     }

    return;
}