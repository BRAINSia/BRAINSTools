//
// Created by Johnson, Hans J on 11/20/16.
//

#include <iostream>
#include <itkImageFileReader.h>

using namespace std;

template <unsigned long int ImageDimension>
void DumpImageInfo(const std::string filename)
{
    typedef double PixelType;
    typedef itk::Image<PixelType,ImageDimension> ImageType;
    typedef typename itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);
    reader->Update();

    typename ImageType::Pointer finalImage = reader->GetOutput();
    cout << finalImage << endl;
}

int main(int argc, char * argv[])
{
    if ( argc != 3 )
    {
    cout << "USAGE: " << argv[0] << " <Dimension> <Filename>" << endl;
    return EXIT_FAILURE;
    }

    const int requestedDim = atoi(argv[1]);
    switch(requestedDim)
    {
    case 2:
        DumpImageInfo<2>(argv[2]);
        break;
    case 3:
        DumpImageInfo<3>(argv[2]);
        break;
    case 4:
        DumpImageInfo<4>(argv[2]);
        break;
    default:
        cout << "ERROR: Invalid dimension choosen" << endl;
        break;
    }


    return EXIT_SUCCESS;
}
