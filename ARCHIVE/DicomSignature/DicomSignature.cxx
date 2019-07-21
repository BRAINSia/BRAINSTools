//
//  This program is designed to make a unique signature of a single dicom file.
//  The intended purpose of this applicaiton is to allow identification of
//  redundant dicom files across directories even if some anonimiation to the
//  dicom header files has occured.
//
//  This example illustrates how to read a single DICOM slice and write it back
//  as another DICOM slice. In the process an intensity rescaling is also
//  applied.
//
//  In order to read and write the slice we use here the \doxygen{GDCMImageIO}
//  class that encapsulates a connection to the underlying GDCM library. In
//  this way we gain access from ITK to the DICOM functionalities offered by
//  GDCM. The GDCMImageIO object is connected as the ImageIO object to be used
//  by the \doxygen{ImageFileWriter}.
//
//  We should first include the following header files.
//
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGDCMImageIO.h"

#include "itkMetaDataObject.h"

#include "gdcmSerieHelper.h"
#include "gdcmImageHelper.h"
#include "gdcmFileExplicitFilter.h"
#include "gdcmImageChangeTransferSyntax.h"
#include "gdcmDataSetHelper.h"
#include "gdcmStringFilter.h"
#include "gdcmImageApplyLookupTable.h"
#include "gdcmImageChangePlanarConfiguration.h"
#include "gdcmRescaler.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmUIDGenerator.h"
#include "gdcmAttribute.h"
#include "gdcmGlobal.h"
#include "gdcmMediaStorage.h"


#include <list>
#include <fstream>
#include <openssl/evp.h>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "DicomSignatureCLP.h"
#include "BRAINSThreadControl.h"

// Taken from comments posted at http://www.digitalpeer.com/id/simple
std::vector< std::string >
Tokenize( const std::string & str, const std::string & delimiters )
{
  std::vector< std::string > tokens;
  std::string::size_type     delimPos = 0, tokenPos = 0, pos = 0;

  if ( str.length() < 1 )
  {
    return tokens;
  }
  while ( 1 )
  {
    delimPos = str.find_first_of( delimiters, pos );
    tokenPos = str.find_first_not_of( delimiters, pos );

    if ( std::string::npos != delimPos )
    {
      if ( std::string::npos != tokenPos )
      {
        if ( tokenPos < delimPos )
        {
          tokens.push_back( str.substr( pos, delimPos - pos ) );
        }
        else
        {
          tokens.push_back( "" );
        }
      }
      else
      {
        tokens.push_back( "" );
      }
      pos = delimPos + 1;
    }
    else
    {
      if ( std::string::npos != tokenPos )
      {
        tokens.push_back( str.substr( pos ) );
      }
      else
      {
        tokens.push_back( "" );
      }
      break;
    }
  }

  return tokens;
}

using DictionaryType = itk::MetaDataDictionary;

static std::string
GetDicomString( const std::string Key, const DictionaryType & dictionary )
{
  using MetaDataStringType = itk::MetaDataObject< std::string >;
  const DictionaryType::ConstIterator end = dictionary.End();
  DictionaryType::ConstIterator       tagItr = dictionary.Find( Key );
  if ( tagItr != end )
  {
    MetaDataStringType::ConstPointer entryvalue =
      dynamic_cast< const MetaDataStringType * >( tagItr->second.GetPointer() );
    // If the dynamic cast succeed, then we can print out the values of the
    // label,
    // the tag and the actual value.
    if ( entryvalue )
    {
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      return tagvalue;
    }
  }
  return std::string( "" );
}

std::string
GetDigestString( EVP_MD_CTX & mdctx )
{
  unsigned char md_value[EVP_MAX_MD_SIZE];
  unsigned int  md_len;

  EVP_DigestFinal_ex( &mdctx, md_value, &md_len );
  EVP_MD_CTX_cleanup( &mdctx );

#define MYHEX( x ) std::setw( 2 ) << std::setfill( '0' ) << std::hex << static_cast< int >( x )

#define MYOCT( x ) std::setw( 2 ) << std::setfill( '0' ) << std::oct << static_cast< int >( x )

  std::ostringstream myDigest( "" );
  // myDigest << digestMode << "_";
  //  myDigest << std::setfill('0') << std::setw(2);//<< std::oct;
  // printf("DIGEST: ");
  for ( unsigned int i = 0; i < md_len; i++ )
  {
    // char tempHex[512];
    // sprintf(tempHex,"%02x", md_value[i]);
    // printf("%02x", md_value[i]);
    myDigest << MYHEX( md_value[i] );
    // std::cout << myDigest.str() << std::endl;
  }
#if 0
  printf("\nDIGEST: ");
  for( unsigned int i = 0; i < md_len; i++ )
    {
    printf("%02x", md_value[i]);
    }
  printf("\n");
#endif
  return myDigest.str();
}

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );
  // Then we declare the pixel type and image dimension, and use them for
  // instantiating the image type to be read.

  using InputPixelType = float;
  constexpr unsigned int InputDimension = 2;

  const std::string fileNamePath = ::itksys::SystemTools::GetFilenamePath( inputVolume );
  const std::string fileNameName = ::itksys::SystemTools::GetFilenameName( inputVolume );

  using InputImageType = itk::Image< InputPixelType, InputDimension >;

  // With the image type we can instantiate the type of the reader, create one,
  // and set the filename of the image to be read.

  using ReaderType = itk::ImageFileReader< InputImageType >;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume );

  // GDCMImageIO is an ImageIO class for reading and writing DICOM v3 and
  // ACR/NEMA images. The GDCMImageIO object is constructed here and connected
  // to
  // the ImageFileReader.

  using ImageIOType = itk::GDCMImageIO;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  reader->SetImageIO( gdcmImageIO );

  // At this point we can trigger the reading process by invoking the Update()
  // method.  Since this reading process may eventually throw an exception, we
  // place the invocation inside a try/catch block.
  //
  bool FileReadSuccessfullyByITK = false;
  try
  {
    reader->Update();
    // This is not reached if an exception is caught.
    FileReadSuccessfullyByITK = true;
  }
  catch ( itk::ExceptionObject & e )
  {
    // Quietly catch.
  }

  OpenSSL_add_all_digests();
  const EVP_MD * md = EVP_get_digestbyname( digestMode.c_str() );
  if ( !md )
  {
    std::cout << "Unknown message digest " << digestMode << std::endl;
    exit( 1 );
  }

  EVP_MD_CTX mdctx;
  EVP_MD_CTX_init( &mdctx );
  EVP_DigestInit_ex( &mdctx, md, NULL );

  if ( FileReadSuccessfullyByITK == true )
  {
    const DictionaryType &  dictionary = gdcmImageIO->GetMetaDataDictionary();
    InputImageType::Pointer myImage = reader->GetOutput();

    std::ostringstream mySpacingInfo;
    mySpacingInfo << myImage->GetSpacing() << std::endl;
    std::ostringstream myOriginInfo;
    myOriginInfo << myImage->GetOrigin() << std::endl;
    std::ostringstream myDirectionInfo;
    myDirectionInfo << myImage->GetDirection() << std::endl;

    // BUILD UP THE DIGEST:  Digest should be the same if all pixel values,
    // spacing, origin, and direction are the same
    // This is to determine if an all zero slice is the same for multiple
    // images.
    // Should add scan date at some point too.
    EVP_DigestUpdate( &mdctx, mySpacingInfo.str().c_str(), mySpacingInfo.str().size() );
    EVP_DigestUpdate( &mdctx, myOriginInfo.str().c_str(), myOriginInfo.str().size() );
    EVP_DigestUpdate( &mdctx, myDirectionInfo.str().c_str(), myDirectionInfo.str().size() );
    // INFO:  Add EVP_DigestUpdate for gdcm strings for
    // <element tag="0008,0021" vr="DA" vm="1" len="8"
    // name="SeriesDate">20090430</element>
    const std::string SeriesDateKey( "0008|0021" );
    const std::string SeriesDateValue = GetDicomString( SeriesDateKey, dictionary );
    EVP_DigestUpdate( &mdctx, SeriesDateValue.c_str(), SeriesDateValue.size() );

#if 0 // Time is not a good indicator for a signature, it is rounded off
      // differently, and often anonymized through truncation
      // <element tag="0008,0031" vr="TM" vm="1" len="14"
      // name="SeriesTime">092135.312000</element>
    const std::string SeriesTimeKey("0008|0031");
    const std::string SeriesTimeValue = GetDicomString(SeriesTimeKey, dictionary);
    EVP_DigestUpdate( &mdctx, SeriesTimeValue.c_str(), SeriesTimeValue.size() );
#else
    // <element tag="0020,0011" vr="IS" vm="1" len="2"
    // name="SeriesNumber">1</element>
    const std::string SeriesNumberKey( "0020|0011" );
    const std::string SeriesNumberValue = GetDicomString( SeriesNumberKey, dictionary );
    EVP_DigestUpdate( &mdctx, SeriesNumberValue.c_str(), SeriesNumberValue.size() );
#endif

    // <element tag="0020,1041" vr="DS" vm="1" len="6"
    // name="SliceLocation">-67.5</element>
    const std::string SliceLocationKey( "0020|1041" );
    const std::string SliceLocationValue = GetDicomString( SliceLocationKey, dictionary );
    EVP_DigestUpdate( &mdctx, SliceLocationValue.c_str(), SliceLocationValue.size() );
    EVP_DigestUpdate( &mdctx,
                      (unsigned char *)myImage->GetBufferPointer(),
                      myImage->GetBufferedRegion().GetNumberOfPixels() * sizeof( InputPixelType ) );

    std::string myDigestString = GetDigestString( mdctx );
    std::cout << myDigestString << "," << digestMode << "," << myImage->GetSpacing()[0] << ","
              << myImage->GetSpacing()[1] << "," << myImage->GetOrigin()[0] << "," << myImage->GetOrigin()[1] << ","
              << myImage->GetDirection()[0][0] << "," << myImage->GetDirection()[0][1] << ","
              << myImage->GetDirection()[1][0] << "," << myImage->GetDirection()[1][1] << "," << SeriesDateValue << ","
#if 0 // Time is not a good indicator for a signature, it is rounded off
      // differently, and often anonymized through truncation
              << SeriesTimeValue << ","
#else
              << SeriesNumberValue << ","
#endif
              << SliceLocationValue << "," << fileNameName << "," << fileNamePath << "," << std::endl;
  }
#define SKIP_DICOM_REPORTS 1
#if SKIP_DICOM_REPORTS

  std::cout << "-1"
            << "," << digestMode << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
            << "-1"
            << ","
#  if 0 // Time is not a good indicator for a signature, it is rounded off
        // differently, and often anonymized through truncation
              << "-1" << ","
#  else
            << "-1"
            << ","
#  endif
            << "-1"
            << "," << fileNameName << "," << fileNamePath << "," << std::endl;
  return EXIT_FAILURE;
#else
  else // Perhaps the file is a dicom report with out image data
  {
    gdcm::SerieHelper * header = new gdcm::SerieHelper;
    header->AddFileName( inputVolume );
    header->SetLoadMode( gdcm::LD_NOSEQ | gdcm::LD_NOSHADOW );
    bool headerLoaded = header->Load();
    if ( !headerLoaded )
    {
#  if 0
      std::cout << "Could not load header from file: "
                << inputVolume << std::endl
                << "Reason: "
                << itksys::SystemTools::GetLastSystemError();
#  endif
      return EXIT_FAILURE;
    }
#  if 0
    if( !header->IsReadable() )
      {
      std::cout << "Could not read header from file: "
                << inputVolume << std::endl
                << "Reason: "
                << itksys::SystemTools::GetLastSystemError();
      return EXIT_FAILURE;
      }
#  endif

    const std::string StudyDate = header->GetEntryValue( 0x0008, 0x0020 );
    if ( StudyDate == "gdcm::Unfound" )
    {
      return EXIT_FAILURE;
    }
    EVP_DigestUpdate( &mdctx, StudyDate.c_str(), StudyDate.size() );
    const std::string SeriesDate = header->GetEntryValue( 0x0008, 0x0021 );
    EVP_DigestUpdate( &mdctx, SeriesDate.c_str(), SeriesDate.size() );
    const std::string StudyTime = header->GetEntryValue( 0x0008, 0x0030 );
    EVP_DigestUpdate( &mdctx, StudyTime.c_str(), StudyTime.size() );
    const std::string SeriesTime = header->GetEntryValue( 0x0008, 0x0031 );
    EVP_DigestUpdate( &mdctx, SeriesTime.c_str(), SeriesTime.size() );
    const std::string ContentTime = header->GetEntryValue( 0x0008, 0x0033 );
    EVP_DigestUpdate( &mdctx, ContentTime.c_str(), ContentTime.size() );

    const std::string EchoTime = header->GetEntryValue( 0x0018, 0x0081 );
    EVP_DigestUpdate( &mdctx, EchoTime.c_str(), EchoTime.size() );

    const std::string EchoNumber = header->GetEntryValue( 0x0018, 0x0086 );
    EVP_DigestUpdate( &mdctx, EchoNumber.c_str(), EchoNumber.size() );

    const std::string InstanceNumber = header->GetEntryValue( 0x0020, 0x0013 );
    EVP_DigestUpdate( &mdctx, InstanceNumber.c_str(), InstanceNumber.size() );

    const std::string SOPClassName = header->GetEntryValue( 0x0008, 0x0016 );
    EVP_DigestUpdate( &mdctx, SOPClassName.c_str(), SOPClassName.size() );
    const std::string Modality = header->GetEntryValue( 0x0008, 0x0060 );
    EVP_DigestUpdate( &mdctx, Modality.c_str(), Modality.size() );
    const std::string SliceLocationKey = header->GetEntryValue( 0x0020, 0x1041 );
    EVP_DigestUpdate( &mdctx, SliceLocationKey.c_str(), SliceLocationKey.size() );
    const std::string SpacingKey = header->GetEntryValue( 0x0028, 0x0030 );
    EVP_DigestUpdate( &mdctx, SpacingKey.c_str(), SpacingKey.size() );
    const std::string OriginKey = header->GetEntryValue( 0x0020, 0x0032 );
    EVP_DigestUpdate( &mdctx, OriginKey.c_str(), OriginKey.size() );
    const std::string DirectionKey = header->GetEntryValue( 0x0020, 0x0037 );
    EVP_DigestUpdate( &mdctx, DirectionKey.c_str(), DirectionKey.size() );

    std::vector< std::string > SpacingVector = Tokenize( SpacingKey, "\\" );
    while ( SpacingVector.size() < 2 )
    {
      SpacingVector.push_back( "0.01" );
    }

    std::vector< std::string > OriginVector = Tokenize( OriginKey, "\\" );
    while ( OriginVector.size() < 2 )
    {
      OriginVector.push_back( "0.02" );
    }

    std::vector< std::string > DirectionVector = Tokenize( DirectionKey, "\\" );
    while ( DirectionVector.size() < 4 )
    {
      DirectionVector.push_back( "0.03" );
    }

    std::string myDigestString = GetDigestString( mdctx );
    std::cout << myDigestString << "," << digestMode << "," << SpacingVector[0] << "," << SpacingVector[1] << ","
              << OriginVector[0] << "," << OriginVector[1] << "," << DirectionVector[0] << "," << DirectionVector[1]
              << "," << DirectionVector[2] << "," << DirectionVector[3] << "," << SeriesDate << "," << SeriesTime << ","
              << "-1.12345"
              << "," << fileNameName << "," << fileNamePath << "," << std::endl;
  }
#endif

  return EXIT_SUCCESS;
}
