/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __castconverthelpers_h__
#define __castconverthelpers_h__

#include "castconvertConfigure.h"

/** Basic Image support. */
#include "itkImage.h"
#include "itkImageIORegion.h"

/** For the support of RGB voxels. */
// #include "itkRGBPixel.h"

/** Reading and writing images. */
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"

/** DICOM headers. */
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

/** One of these is used to cast the image. */
#include "itkShiftScaleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#ifdef VTK_FOUND
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkXMLImageDataReader.h"
#include "vtkMetaImageWriter.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#endif

/** Print image information from the reader and the writer. */
template <class ReaderType, class WriterType>
void PrintInfo( ReaderType reader, WriterType writer )
{
  /** Typedef's. */
  typedef itk::ImageIOBase                     ImageIOBaseType;
  typedef itk::ImageIORegion                   ImageIORegionType;
  typedef typename ImageIORegionType::SizeType SizeType;

  /** Get IOBase of the reader and extract information. */
  ImageIOBaseType::Pointer imageIOBaseIn = reader->GetImageIO();
  ImageIORegionType        iORegionIn = imageIOBaseIn->GetIORegion();

  const char * fileNameIn = imageIOBaseIn->GetFileName();
  std::string  pixelTypeIn = imageIOBaseIn->GetPixelTypeAsString( imageIOBaseIn->GetPixelType() );
  unsigned int nocIn = imageIOBaseIn->GetNumberOfComponents();
  std::string  componentTypeIn = imageIOBaseIn->GetComponentTypeAsString( imageIOBaseIn->GetComponentType() );
  unsigned int dimensionIn = imageIOBaseIn->GetNumberOfDimensions();
  SizeType     sizeIn = iORegionIn.GetSize();

  /**  Get  IOBase of  the  writer and extract information.  */
  ImageIOBaseType::Pointer imageIOBaseOut = writer->GetImageIO();
  ImageIORegionType        iORegionOut = imageIOBaseOut->GetIORegion();

  const char * fileNameOut = imageIOBaseOut->GetFileName();
  std::string  pixelTypeOut = imageIOBaseOut->GetPixelTypeAsString( imageIOBaseOut->GetPixelType() );
  unsigned int nocOut = imageIOBaseOut->GetNumberOfComponents();
  std::string  componentTypeOut = imageIOBaseOut->GetComponentTypeAsString( imageIOBaseOut->GetComponentType() );
  unsigned int dimensionOut = imageIOBaseOut->GetNumberOfDimensions();
  SizeType     sizeOut = iORegionOut.GetSize();

  /** Print information. */
  std::cout << "Information about the input image \"" << fileNameIn << "\":" << std::endl;
  std::cout << "\tdimension:\t\t" << dimensionIn << std::endl;
  std::cout << "\tpixel type:\t\t" << pixelTypeIn << std::endl;
  std::cout << "\tnumber of components:\t" << nocIn << std::endl;
  std::cout << "\tcomponent type:\t\t" << componentTypeIn << std::endl;
  std::cout << "\tsize:\t\t\t";
  for( unsigned int i = 0; i < dimensionIn; i++ )
    {
    std::cout << sizeIn[i] << " ";
    }
  std::cout << std::endl;

  /** Print information. */
  std::cout << std::endl;
  std::cout << "Information about the output image \"" << fileNameOut << "\":" << std::endl;
  std::cout << "\tdimension:\t\t" << dimensionOut << std::endl;
  std::cout << "\tpixel type:\t\t" << pixelTypeOut << std::endl;
  std::cout << "\tnumber of components:\t" << nocOut << std::endl;
  std::cout << "\tcomponent type:\t\t" << componentTypeOut << std::endl;
  std::cout << "\tsize:\t\t\t";
  for( unsigned int i = 0; i < dimensionOut; i++ )
    {
    std::cout << sizeOut[i] << " ";
    }
  std::cout << std::endl;
}  // end PrintInfo

/** The function that reads the input dicom image and writes the output image.
 * This function is templated over the image types. In the main function
 * we have to make sure to call the right instantiation.
 */
template <class InputImageType, class OutputImageType>
void ReadDicomSeriesCastWriteImage( std::string inputDirectoryName, std::string  outputFileName )
{
  /** Typedef the correct reader, caster and writer. */
  typedef typename itk::ImageSeriesReader<InputImageType> SeriesReaderType;
  typedef typename itk::ShiftScaleImageFilter<InputImageType, OutputImageType> ShiftScaleFilterType;
  typedef typename itk::ImageFileWriter<OutputImageType> ImageWriterType;

  /** Typedef dicom stuff. */
  typedef itk::GDCMImageIO         GDCMImageIOType;
  typedef itk::GDCMSeriesFileNames GDCMNamesGeneratorType;
  typedef std::vector<std::string> FileNamesContainerType;

  /** Create the dicom ImageIO. */
  typename GDCMImageIOType::Pointer dicomIO = GDCMImageIOType::New();

  /** Get a list of the filenames of the 2D input dicom images. */
  GDCMNamesGeneratorType::Pointer nameGenerator = GDCMNamesGeneratorType::New();
  nameGenerator->SetInputDirectory( inputDirectoryName.c_str() );
  FileNamesContainerType fileNames = nameGenerator->GetInputFileNames();

  /** Create and setup the seriesReader. */
  typename SeriesReaderType::Pointer seriesReader = SeriesReaderType::New();
  seriesReader->SetFileNames( fileNames );
  seriesReader->SetImageIO( dicomIO );

  /** Create and setup caster and writer. */
  // typename RescaleFilterType::Pointer caster = RescaleFilterType::New();
  typename ShiftScaleFilterType::Pointer caster = ShiftScaleFilterType::New();
  typename ImageWriterType::Pointer  writer = ImageWriterType::New();
  caster->SetShift( 0.0 );
  caster->SetScale( 1.0 );
  writer->SetFileName( outputFileName.c_str()  );

  /** Connect the pipeline. */
  caster->SetInput(  seriesReader->GetOutput()  );
  writer->SetInput(  caster->GetOutput()  );
  writer->UseCompressionOn();

  // Handle .vti files as well.
#ifdef VTK_FOUND
  if( outputFileName.rfind(".vti") == (outputFileName.size() - 4) )
    {
    typedef itk::ImageToVTKImageFilter<OutputImageType> ITKToVTKFilterType;
    typename ITKToVTKFilterType::Pointer itktovtk = ITKToVTKFilterType::New();
    caster->Update();
    itktovtk->SetInput( caster->GetOutput() );
    itktovtk->Update();

    vtkSmartPointer<vtkXMLImageDataWriter> writer_vti
      = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer_vti->SetFileName(outputFileName.c_str() );
    writer_vti->SetInput( itktovtk->GetOutput() );

    // necessary to give the data array a name others reading it fails !!
    itktovtk->GetOutput()->GetPointData()->GetScalars()->SetName("Scalars_");

    writer_vti->Write();
    std::cout << "Wrote: " << outputFileName << std::endl;
    return;
    }
#endif

  /**  Do the actual  conversion.  */
  writer->Update();

  /**  Print  information. */
  PrintInfo( seriesReader, writer );
}  // end ReadDicomSeriesCastWriteImage

/** The function that reads the input image and writes the output image.
 * This function is templated over the image types. In the main function
 * we have to make sure to call the right instantiation.
 */
template <class InputImageType, class OutputImageType>
void ReadCastWriteImage( std::string inputFileName, std::string outputFileName )
{
  /**  Typedef the correct reader, caster and writer. */
  typedef typename itk::ImageFileReader<InputImageType> ImageReaderType;
  typedef typename itk::ShiftScaleImageFilter<InputImageType, OutputImageType> ShiftScaleFilterType;
  typedef typename itk::ImageFileWriter<OutputImageType> ImageWriterType;
#ifdef VTK_FOUND
  typedef itk::VTKImageToImageFilter<InputImageType> VTKToITKFilterType;
  typename VTKToITKFilterType::Pointer vtktoitk = NULL;
#endif
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  typename ShiftScaleFilterType::Pointer caster = ShiftScaleFilterType::New();
  typename ImageWriterType::Pointer  writer = ImageWriterType::New();
  caster->SetShift( 0.0 );
  caster->SetScale( 1.0 );

  /** Create and setup the reader. */
  reader->SetFileName( inputFileName.c_str() );

#ifdef VTK_FOUND
  if( inputFileName.rfind(".vti") == (inputFileName.size() - 4) )
    {
    // Handle .vti files as well.
    vtkSmartPointer<vtkXMLImageDataReader> reader_vti
      = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader_vti->SetFileName(inputFileName.c_str() );
    reader_vti->Update();
    std::cout << "Read: " << inputFileName << std::endl;

    vtktoitk = VTKToITKFilterType::New();
    vtktoitk->SetInput( reader_vti->GetOutput() );
    vtktoitk->Update();
    caster->SetInput(vtktoitk->GetOutput() );
    }
  else
    {
#endif

  // typename RescaleFilterType::Pointer caster = RescaleFilterType::New();
  caster->SetInput( reader->GetOutput() );

#ifdef VTK_FOUND
}

if( outputFileName.rfind(".vti") == (outputFileName.size() - 4) )
  {
  // Handle .vti files as well.
  typedef itk::ImageToVTKImageFilter<OutputImageType> ITKToVTKFilterType;
  typename ITKToVTKFilterType::Pointer itktovtk = ITKToVTKFilterType::New();
  caster->Update();
  itktovtk->SetInput( caster->GetOutput() );
  itktovtk->Update();

  vtkSmartPointer<vtkXMLImageDataWriter> writer_vti
    = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer_vti->SetFileName(outputFileName.c_str() );
  writer_vti->SetInput( itktovtk->GetOutput() );

  // necessary to give the data array a name others reading it fails !!
  itktovtk->GetOutput()->GetPointData()->GetScalars()->SetName("Scalars_");

  writer_vti->Write();
  std::cout << "Wrote: " << outputFileName << std::endl;
  return;
  }
#endif

  writer->UseCompressionOn();
  writer->SetFileName( outputFileName.c_str()  );
  writer->SetInput( caster->GetOutput() );
  writer->Update();

  if( inputFileName.rfind(".vti") != (inputFileName.size() - 4) )
    {
    /** Print information. */
    PrintInfo( reader, writer );
    }
}  // end ReadWriteImage

// Read VTI files and write out ITK/VTI images. We have some templated macros
// to deal with since VTK does not have templated images while ITK does.
#ifdef VTK_FOUND

// ITK's support for 64 bit types, long long etc is poor, not as exhaustive
// as VTK. Define an alternate macro here that ignores those types. Nobody
// will use them anyway.
#define vtkitkTemplateMacro(call)                                           \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call);                           \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);                             \
  vtkTemplateMacroCase(VTK_LONG, long, call);                               \
  vtkTemplateMacroCase(VTK_UNSIGNED_LONG, unsigned long, call);             \
  vtkTemplateMacroCase(VTK_INT, int, call);                                 \
  vtkTemplateMacroCase(VTK_UNSIGNED_INT, unsigned int, call);               \
  vtkTemplateMacroCase(VTK_SHORT, short, call);                             \
  vtkTemplateMacroCase(VTK_UNSIGNED_SHORT, unsigned short, call);           \
  vtkTemplateMacroCase(VTK_CHAR, char, call);                               \
  vtkTemplateMacroCase(VTK_SIGNED_CHAR, signed char, call);                 \
  vtkTemplateMacroCase(VTK_UNSIGNED_CHAR, unsigned char, call)

template <class TInputPixelType, class TOutputPixelType>
int
ReadVTICastWriteImage( std::string inputFileName,
                       std::string outputFileName,
                       TInputPixelType, TOutputPixelType )
{
  typedef itk::Image<TInputPixelType,  3> InputImageType;
  typedef itk::Image<TOutputPixelType, 3> OutputImageType;
  ReadCastWriteImage<InputImageType, OutputImageType>(
    inputFileName, outputFileName );
  return 1;
}

template <class TOutputPixelType>
int
ReadVTICastWriteImage( std::string inputFileName,
                       std::string outputFileName,
                       int dimension )
{
  int retval = 0;

  if( inputFileName.rfind(".vti") == (inputFileName.size() - 4) && dimension == 3 )
    {
    vtkSmartPointer<vtkXMLImageDataReader> reader =
      vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName( inputFileName.c_str() );
    reader->Update();
    vtkImageData *   image = reader->GetOutput();
    TOutputPixelType op = static_cast<TOutputPixelType>(0);

    switch( image->GetScalarType() )
      {
      vtkitkTemplateMacro( retval =
                             ReadVTICastWriteImage( inputFileName, outputFileName, static_cast<VTK_TT>(0), op ) );

      default:
        {
        std::cout << "VTI conversion.. Unknown scalar type" << std::endl;
        break;
        }
      }
    }

  return retval;
}

#else
template <class TOuptutPixelType>
int
ReadVTICastWriteImage( std::string,
                       std::string,
                       int )
{
  return 0;
}

#endif

/** Macros are used in order to make the code in main() look cleaner. */

/** callCorrectReadDicomWriterMacro:
 * A macro to call the dicom-conversion function.
 */

#define callCorrectReadDicomWriterMacro(typeIn, typeOut) \
  if( inputPixelComponentType == #typeIn && outputPixelComponentType == #typeOut ) \
    { \
    typedef  itk::Image<typeIn, 3>  InputImageType; \
    typedef  itk::Image<typeOut, 3> OutputImageType; \
    ReadDicomSeriesCastWriteImage<InputImageType, OutputImageType>( inputDirectoryName, outputFileName ); \
    }

/** callCorrectReadWriterMacro:
 * A macro to call the conversion function.
 */

#define callCorrectReadWriterMacro(typeIn, typeOut, dim) \
  if( inputPixelComponentType == #typeIn && outputPixelComponentType == #typeOut && inputDimension == dim ) \
    { \
    if( !ReadVTICastWriteImage<typeOut>( inputFileName, outputFileName, dim ) ) \
      { \
      typedef  itk::Image<typeIn, dim>  InputImageType; \
      typedef  itk::Image<typeOut, dim> OutputImageType; \
      ReadCastWriteImage<InputImageType, OutputImageType>( inputFileName, outputFileName ); \
      } \
    }

#endif // __castconverthelpers_h__
