#ifndef BRAINSTOOLS_BRAINSRefacerUtilityFunctions_H
#define BRAINSTOOLS_BRAINSRefacerUtilityFunctions_H

#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkRoundImageFilter.h>
#include <itkCastImageFilter.h>

//Convenience function to write images
template< typename TImageType >
void WriteImage(std::string filename, TImageType *image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}

//Convenience function to write images
template< typename TImageType >
void WriteSmartImage(std::string filename, typename TImageType::Pointer image)
{
  std::cout << "Writing Image: " << filename << std::endl;
  typedef itk::ImageFileWriter<TImageType> FileWriterType;
  typename FileWriterType::Pointer fileWriter = FileWriterType::New();

  fileWriter->SetInput(image);
  fileWriter->SetFileName(filename);
  fileWriter->Update();
  std::cout << "\tdone writing Image: " << filename << std::endl;
}


//Convienience function to write transforms
template< typename TTransformType >
void WriteTransform(std::string transformFileName, TTransformType transform )
{
  std::cout << "Writing Transform: " << transformFileName << std::endl;
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  transformWriter->SetInput(transform);
  transformWriter->SetFileName(transformFileName);
  transformWriter->Update();
  std::cout << "\t done writing Transform: " << transformFileName << std::endl;
}


template<typename TInputImageType, typename TOutputImageType>
void RoundAndWriteImage( std::string fileName, typename TInputImageType::Pointer image)
{
  auto rounder = itk::RoundImageFilter<TInputImageType, TOutputImageType>::New();
  rounder->SetInput(image);
  WriteSmartImage<TOutputImageType>(fileName, rounder->GetOutput());
}



template<typename TInputImageType, typename TOutputImageType>
void CastAndWriteImage( std::string fileName, typename TInputImageType::Pointer image)
{
  auto caster = itk::CastImageFilter<TInputImageType, TOutputImageType>::New();
  caster->SetInput(image);
  WriteSmartImage<TOutputImageType>(fileName, caster->GetOutput());
}


template< typename TInputImageType, int TDimension>
int ConvertAndSave(std::string fileName, typename TInputImageType::Pointer image, itk::ImageIOBase::IOComponentType OutputComponentType_ENUM)
{
  switch (OutputComponentType_ENUM)
    {

    //INTEGER like types  --- round and write
    case itk::ImageIOBase::SHORT:
      {
      typedef itk::Image<short, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UCHAR:
      {
      typedef itk::Image<unsigned char, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::CHAR:
      {
      typedef itk::Image<char, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::USHORT:
      {
      typedef itk::Image<unsigned short, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UINT:
      {
      typedef itk::Image<unsigned int, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::INT:
      {
      typedef itk::Image<int, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::ULONG:
      {
      typedef itk::Image<unsigned long, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::LONG:
      {
      typedef itk::Image<long, TDimension> OutputType;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }


      //FLOAT LIKE TYPES --- cast and write
    case itk::ImageIOBase::FLOAT:
      {
      typedef itk::Image<float, TDimension> OutputType;
      CastAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::DOUBLE:
      {
      typedef itk::Image<double, TDimension> OutputType;
      CastAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }


      //ERROR producing cases
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      {
      std::cerr << "Unknown component type: ..." << std::endl;
      return EXIT_FAILURE;
      }

    default:
      {
      std::cerr << "Unprocessed Case: " << itk::ImageIOBase::GetComponentTypeAsString(OutputComponentType_ENUM) <<
      std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

#endif