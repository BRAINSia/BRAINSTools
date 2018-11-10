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
  using FileWriterType = itk::ImageFileWriter<TImageType>;
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
  using FileWriterType = itk::ImageFileWriter<TImageType>;
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
  using TransformWriterType = itk::TransformFileWriter;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  transformWriter->SetInput(transform);
  transformWriter->SetFileName(transformFileName);
#if ITK_VERSION_MAJOR >= 5
  transformWriter->SetUseCompression(true);
#endif
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
      using OutputType = itk::Image<short, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UCHAR:
      {
      using OutputType = itk::Image<unsigned char, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::CHAR:
      {
      using OutputType = itk::Image<char, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::USHORT:
      {
      using OutputType = itk::Image<unsigned short, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::UINT:
      {
      using OutputType = itk::Image<unsigned int, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::INT:
      {
      using OutputType = itk::Image<int, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::ULONG:
      {
      using OutputType = itk::Image<unsigned long, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::LONG:
      {
      using OutputType = itk::Image<long, TDimension>;
      RoundAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }


      //FLOAT LIKE TYPES --- cast and write
    case itk::ImageIOBase::FLOAT:
      {
      using OutputType = itk::Image<float, TDimension>;
      CastAndWriteImage<TInputImageType, OutputType>(fileName, image);
      break;
      }

    case itk::ImageIOBase::DOUBLE:
      {
      using OutputType = itk::Image<double, TDimension>;
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
