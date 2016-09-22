/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
 */


#include "mex.h"
#include "itk3DImageToMxArray.h"
#include <iostream>

// Read ITK Scalar 3D Image

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    if (1 != nlhs || 1 != nrhs || !mxIsChar(prhs[0])) {
        std::cerr << "Using Parameters error." << std::endl;
        std::cerr << "Usage: outputImageArray = itkRead3DImage('Scalar3DImageFilename') " <<std::endl;
        std::cerr << "This program only supports ITK scalar 3D image." << std::endl;
        mexErrMsgIdAndTxt("MATLAB::itkRead3DImage:ErrorParameter", "Usage: imageArray = itkReadImage('Scalar3DImageFilename') \n");
        // it returns control to the MATLAB prompt.
    }
    else {

        int length=(int)mxGetN(prhs[0])+1;
        char  *filenameBuf = new char[length];

        int status=mxGetString(prhs[0],filenameBuf,(mwSize)length);
        std::string filename(filenameBuf);
        delete[] filenameBuf;

        if (0 != status || 0 == length) {
            std::cerr << "Error in input image file name parameter." << std::endl;
            mexErrMsgIdAndTxt("MATLAB::itkRead3DImage:ErrorParameter", "illegal input image Filename.\n");
        }
        else{
            itk3DImage2MxArray(filename, plhs);
        }

    }
}