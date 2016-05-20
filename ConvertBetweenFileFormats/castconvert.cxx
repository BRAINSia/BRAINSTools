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
/**
 * castconvert
 *
 * This program converts and possibly casts images.
 *
 * This is done by reading in an image, possibly casting of the image,
 * and subsequently writing the image to some format.
 * With converting we mean changing the extension of the image,
 * such as bmp, mhd, etc. With casting we mean changing the component
 * type of a voxel, such as short, unsigned long, float.
 *
 * Casting is currently done using the ShiftScaleImageFilter,
 * where values are mapped to itself, leaving the intensity range
 * the same. NOTE that when casting to a component type with a
 * smaller dynamic range, information might get lost. In this case
 * we might use the RescaleIntensityImageFilter to linearly
 * rescale the image values.
 *
 * Currently only supported are the SCALAR pixel types.
 * Input images can be in all file formats ITK supports and for which
 * the ImageFileReader works, and additionally single 3D dicom series
 * using the ImageSeriesReader. The pixel component type should
 * of course be a component type supported by the file format.
 * Output images can be in all file formats ITK supports and for which
 * the ImageFileReader works, so no dicom output is currently supported.
 *
 * authors:       Marius Staring and Stefan Klein
 *
 * Thanks to Hans J. Johnson for a modification to this program. This
 * modification breaks down the program into smaller compilation units,
 * so that the compiler does not overflow.
 *
 */

#include <iostream>
#include "itkImageFileReader.h"

/** DICOM headers. */
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"

/** In order to determine if argv[1] is a directory or a file,
 * so that we can distinguish between dicom and other files.
 */
#include <itksys/SystemTools.hxx>

extern int      FileConverterScalar(const std::string & inputPixelComponentType,
                                    const std::string & outputPixelComponentType, const std::string & inputFileName,
                                    const std::string & outputFileName, int inputDimension);

extern int DicomFileConverterScalar(const std::string & inputPixelComponentType,
                                    const std::string & outputPixelComponentType, const std::string & inputFileName,
                                    const std::string & outputFileName, int inputDimension);

extern int DicomFileConverterScalarA(const std::string & inputPixelComponentType,
                                     const std::string & outputPixelComponentType, const std::string & inputFileName,
                                     const std::string & outputFileName, int inputDimension);

// -------------------------------------------------------------------------------------
#include "itkGE4ImageIOFactory.h"
#include "itkGE5ImageIOFactory.h"
#include "itkGEAdwImageIOFactory.h"


int  main(  int  argc,  char *argv[] )
{
  try { //Whole program try-catch block for getting failures
    itk::GE4ImageIOFactory::RegisterOneFactory();
    itk::GE5ImageIOFactory::RegisterOneFactory();
    itk::GEAdwImageIOFactory::RegisterOneFactory();

    /** TASK 1:
     * Check arguments.
     * *******************************************************************
     */
    if( argc  < 3 || (strcmp(argv[1], "--help") == 0) )
      {
      std::cout << "Usage:"  << std::endl;
      std::cout << "\tcastconvert inputfilename outputfilename [outputPixelComponentType]" << std::endl;
      std::cout << "\tcastconvert dicomDirectory outputfilename [outputPixelComponentType]" << std::endl;
      std::cout << "\twhere outputPixelComponentType is one of:" << std::endl;
      std::cout << "\t\t- unsigned_char" << std::endl;
      std::cout << "\t\t- char" << std::endl;
      std::cout << "\t\t- unsigned_short" << std::endl;
      std::cout << "\t\t- short" << std::endl;
      std::cout << "\t\t- unsigned_int" << std::endl;
      std::cout << "\t\t- int" << std::endl;
      std::cout << "\t\t- unsigned_long" << std::endl;
      std::cout << "\t\t- long" << std::endl;
      std::cout << "\t\t- float" << std::endl;
      std::cout << "\t\t- double" << std::endl;
      std::cout << "\tprovided that the outputPixelComponentType is supported by the output file format." << std::endl;
      std::cout << "\tBy default the outputPixelComponentType is set to the inputPixelComponentType." << std::endl;
      return EXIT_FAILURE;
      }

    /**  Get  the  inputs. */
    std::string input = argv[1];
    std::string outputFileName = argv[2];
    std::string outputPixelComponentType = "";
    if( argc >= 4 )
      {
      outputPixelComponentType = argv[3];
      }

    /** Make sure last character of input != "/".
     * Otherwise FileIsDirectory() won't work.
     */
    if( input.rfind( "/" ) == input.size() - 1 )
      {
      input.erase( input.size() - 1, 1 );
      }

    /** Check if input is a file or a directory. */
    bool        exists = itksys::SystemTools::FileExists( input.c_str() );
    bool        isDir = itksys::SystemTools::FileIsDirectory( input.c_str() );
    bool        isVTI = (input.rfind(".vti") == (input.size() - 4) );
    bool        isDICOM = false;
    std::string inputFileName, inputDirectoryName;

    if( exists && !isDir )
      {
      /** Input is a file, and we use the ImageFileReader. */
      inputFileName = input;
      }
    else if( exists && isDir )
      {
      /** Input is a directory, and we use the ImageSeriesReader. */
      inputDirectoryName = input;
      isDICOM = true;
      }
    else
      {
      /** Something is wrong. */
      std::cerr << "ERROR: first input argument does not exist!" << std::endl;
      return EXIT_FAILURE;
      }

    /** Check outputPixelType. */
    if( outputPixelComponentType != ""
      && outputPixelComponentType != "unsigned_char"
      && outputPixelComponentType != "char"
      && outputPixelComponentType != "unsigned_short"
      && outputPixelComponentType != "short"
      && outputPixelComponentType != "unsigned_int"
      && outputPixelComponentType != "int"
      && outputPixelComponentType != "unsigned_long"
      && outputPixelComponentType != "long"
      && outputPixelComponentType != "float"
      && outputPixelComponentType != "double" )
      {
      /** In this case an illegal outputPixelComponentType is given. */
      std::cerr << "The given outputPixelComponentType is \"" << outputPixelComponentType
        << "\", which is not supported." << std::endl;
      return EXIT_FAILURE;
      }

    /** TASK 2:
     * Typedefs and test reading to determine correct image types.
     * *******************************************************************
     */

    /** Initial image type. */
    const unsigned int Dimension  =  3;
    typedef short PixelType;

    /** Some typedef's. */
    typedef itk::Image<PixelType, Dimension> ImageType;
    typedef itk::ImageFileReader<ImageType>  ReaderType;
    typedef itk::ImageIOBase                 ImageIOBaseType;
    typedef itk::GDCMImageIO                 GDCMImageIOType;
    typedef itk::GDCMSeriesFileNames         GDCMNamesGeneratorType;
    typedef std::vector<std::string>         FileNamesContainerType;

    /** Create a testReader. */
    ReaderType::Pointer testReader = ReaderType::New();

    /** Setup the testReader. */
    if( !isDICOM && !isVTI )
      {
      /** Set the inputFileName in the testReader. */
      testReader->SetFileName( inputFileName.c_str() );
      }
    else if( !isVTI )
      {
      /** Get a name of a 2D image. */
      GDCMNamesGeneratorType::Pointer nameGenerator = GDCMNamesGeneratorType::New();
      nameGenerator->SetInputDirectory( inputDirectoryName.c_str() );
      FileNamesContainerType fileNames = nameGenerator->GetInputFileNames();
      std::string            fileName = fileNames[0];

      /** Create a dicom ImageIO and set it in the testReader. */
      GDCMImageIOType::Pointer dicomIO = GDCMImageIOType::New();
      testReader->SetImageIO( dicomIO );

      /** Set the name of the 2D dicom image in the testReader. */
      testReader->SetFileName( fileName.c_str() );
      } // end isDICOM

    // The defaults are arbitrary. They are computed from GenerateOutputInformation
    // below. For VTI files, these are computed from the actual input image itself.
    unsigned int inputDimension = 3;
    unsigned int numberOfComponents = 1;
    std::string  inputPixelComponentType = "short";
    std::string  pixelType = "scalar";

    /** Get the component type, number of components, dimension and pixel type. */
    if( !isVTI )
      {
      /** Generate all information. */
      testReader->Update();

      /** Extract the ImageIO from the testReader. */
      ImageIOBaseType::Pointer testImageIOBase = testReader->GetImageIO();

      numberOfComponents = testImageIOBase->GetNumberOfComponents();
      inputDimension = testImageIOBase->GetNumberOfDimensions();
      inputPixelComponentType = testImageIOBase->GetComponentTypeAsString(
        testImageIOBase->GetComponentType() );
      pixelType = testImageIOBase->GetPixelTypeAsString( testImageIOBase->GetPixelType() );
      }

    /** TASK 3:
     * Do some preparations.
     * *******************************************************************
     */

    /** Check inputPixelType. */
    if( inputPixelComponentType != "unsigned_char"
      && inputPixelComponentType != "char"
      && inputPixelComponentType != "unsigned_short"
      && inputPixelComponentType != "short"
      && inputPixelComponentType != "unsigned_int"
      && inputPixelComponentType != "int"
      && inputPixelComponentType != "unsigned_long"
      && inputPixelComponentType != "long"
      && inputPixelComponentType != "float"
      && inputPixelComponentType != "double" )
      {
      /** In this case an illegal inputPixelComponentType is found. */
      std::cerr << "The found inputPixelComponentType is \"" << inputPixelComponentType
        << "\", which is not supported." << std::endl;
      return EXIT_FAILURE;
      }

    /** Check outputPixelType. */
    if( outputPixelComponentType == "" )
      {
      /** In this case this option is not given, and by default
       * we set it to the inputPixelComponentType.
       */
      outputPixelComponentType = inputPixelComponentType;
      }

    /** Get rid of the "_" in inputPixelComponentType and outputPixelComponentType. */
    std::basic_string<char>::size_type              pos = inputPixelComponentType.find( "_" );
    static const std::basic_string<char>::size_type npos = std::basic_string<char>::npos;
    if( pos != npos )
      {
      inputPixelComponentType.replace( pos, 1, " " );
      }
    pos = outputPixelComponentType.find( "_" );
    if( pos != npos )
      {
      outputPixelComponentType.replace( pos, 1, " " );
      }

    /** TASK 4:
     * Now we are ready to check on image type and subsequently call the
     * correct ReadCastWrite-function.
     * *******************************************************************
     */

    try
      {
      if( !isDICOM )
        {
        /**
         * ****************** Support for SCALAR pixel types. **********************************
         */
        if( pixelType == "scalar" && numberOfComponents == 1 )
          {
          const int ret_value = FileConverterScalar(
            inputPixelComponentType, outputPixelComponentType, inputFileName,
            outputFileName, inputDimension );
          if( ret_value != 0 )
            {
            return ret_value;
            }
          }
        else
          {
          std::cerr << "Pixel type is " << pixelType
            << ", component type is " << inputPixelComponentType
            << " and number of components equals " << numberOfComponents << "." << std::endl;
          std::cerr << "ERROR: This image type is not supported." << std::endl;
          return EXIT_FAILURE;
          }
        } // end NonDicom image
      else
        {
        /** In this case input is a DICOM series, from which we only support
         * SCALAR pixel types, with component type:
         * DICOMImageIO2: (unsigned) char, (unsigned) short, float
         * GDCMImageIO: (unsigned) char, (unsigned) short, (unsigned) int, double
         * It is also assumed that the dicom series consist of multiple
         * 2D images forming a 3D image.
         */

        if( pixelType == "scalar" && numberOfComponents == 1 )
          {
          const int ret_value = DicomFileConverterScalar(
            inputPixelComponentType, outputPixelComponentType,
            inputDirectoryName, outputFileName, inputDimension )
            ||                  DicomFileConverterScalarA(
              inputPixelComponentType, outputPixelComponentType,
              inputDirectoryName, outputFileName, inputDimension );
          if( ret_value != 0 )
            {
            return ret_value;
            }
          }
        else
          {
          std::cerr << "Pixel type is " << pixelType
            << ", component type is " << inputPixelComponentType
            << " and number of components equals " << numberOfComponents << "." << std::endl;
          std::cerr << "ERROR: This DICOM image type is not supported." << std::endl;
          return EXIT_FAILURE;
          }
        } // end isDICOM
      }   // end try
    catch( itk::ExceptionObject  &  err  )
      {
      /** If any errors have occurred, catch and print the exception and return false. */
      std::cerr << "ExceptionObject caught !"  << std::endl;
      std::cerr << err <<  std::endl;
      return EXIT_FAILURE;
      }
  }
  catch( itk::ExceptionObject  &  err  )
    {
    /** If any errors have occurred, catch and print the exception and return false. */
    std::cerr << "ExceptionObject caught !"  << std::endl;
    std::cerr << err <<  std::endl;
    return EXIT_FAILURE;
    }
  /** End  program. Return succes. */
  return EXIT_SUCCESS;
}  // end main
