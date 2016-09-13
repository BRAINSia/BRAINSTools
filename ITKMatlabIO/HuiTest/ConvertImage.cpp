//
// Created by Hui Xie on 9/8/16.
//

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

//#include <QuickView.h>
#include "ConvertImage.h"

static const size_t DIMENSION = 3;

//int myvar=77;

int ConvertImage(const std::string inputFileName,const std::string outputFileName)
{

    typedef itk::Image< double, DIMENSION >  ImageType;
    typedef itk::ImageFileReader<ImageType>  ReaderType;
    typedef itk::ImageFileWriter<ImageType>  WriterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inputFileName);

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFileName);
    writer->SetInput(reader->GetOutput());

    writer->Update();

    return EXIT_SUCCESS;
}