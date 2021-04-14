/*
 * Author: Hui Xie    SheenXH@gmail.com
 * date: Sep 20th, 2016
 */


#include "mex.h"
#include "convertItkImageToMxArray.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    if (1 != nlhs || 1 != nrhs || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB::itkReadImage:ErrorParameter", "Usage: imageArray = itkReadImage(imageFilename) \n");
        // it returns control to the MATLAB prompt.
    }
    else {

        int length=(int)mxGetN(prhs[0])+1;
        char  *filenameBuf = new char[length];

        int status=mxGetString(prhs[0],filenameBuf,(mwSize)length);
        std::string filename(filenameBuf);
        delete[] filenameBuf;

        if (0 != status || 0 == length) {
            mexErrMsgIdAndTxt("MATLAB::itkReadImage:ErrorParameter", "illegal input imageFilename.\n");
        }
        else{
            itkImage2MxArray(filename, plhs);
        }

    }
}