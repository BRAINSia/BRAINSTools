//
// Created by Hui Xie on 9/8/16.
//

#include "ConvertImage.h"
#include <iostream>

int main(int argc, char *argv[])
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " inputImageFile" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string inputFileName(argv[1]);
    const std::string outputFileName("/tmp/test.nrrd");

    //std::cout << myvar << std::endl;

    return ConvertImage(inputFileName,outputFileName);
}