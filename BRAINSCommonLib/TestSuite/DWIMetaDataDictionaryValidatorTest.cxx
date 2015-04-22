//TODO Convert itkDwiToVectorImageFilter.hxx to use this class.

/**
 * \author Ali Ghayoor
 * This test file creates an ITK vector image,
 * and writes the image to the disk in NRRD format.
 */

#include <itkMetaDataObject.h>

#include <cmath>
#include <cstdio>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkVectorImage.h>

#include "DWIMetaDataDictionaryValidator.h"

//Force printing of strings nicely
static std::string ForceConvert( const itk::MetaDataObjectBase *  myMetaDataObjectBase)
{
  const itk::MetaDataObject<std::string> * myMetaDataObject = dynamic_cast<const itk::MetaDataObject<std::string> * >(myMetaDataObjectBase);
  if(myMetaDataObject)
    {
    //const std::string temp  =myMetaDataObject->GetMetaDataObjectTypeName();
    const std::string temp("std::string");
    return myMetaDataObject->GetMetaDataObjectValue()+"         (type): " + temp; //assuming you define print for matrix
    }
  else{
    const std::string temp  = myMetaDataObjectBase->GetMetaDataObjectTypeName();
    return "???_conversion         (type): " + temp;
  }
}


static void PrintDictionaryHelper(const itk::MetaDataDictionary & dictPrint)
{

  //  Print keys and values?
  std::vector<std::string> keys = dictPrint.GetKeys();


  std::cout << "----------------" << std::endl;
  itk::MetaDataDictionary::ConstIterator end = dictPrint.End();
  for ( itk::MetaDataDictionary::ConstIterator it = dictPrint.Begin(); it != end; ++it)
    {
    std::cout << ' ' << it->first << ":=" << ForceConvert( it->second ) << std::endl;
    }
}


#define KEY_PREFIX "NRRD_"

typedef short                           PixelValue;
typedef itk::VectorImage<PixelValue, 3> ImageType;
typedef std::vector<std::string>        StringVectorType;

static ImageType::Pointer CreateVolume(const size_t numOfComponents)
{
  const int imageSize = 10; // each image component has size of imageSize^3

  ImageType::IndexType start;
  start.Fill(0);

  ImageType::SizeType size;
  size.Fill(imageSize);

  ImageType::RegionType region(start, size);

  ImageType::Pointer nrrdVolume = ImageType::New();
  nrrdVolume->SetRegions(region);
  nrrdVolume->SetVectorLength(numOfComponents);
  nrrdVolume->Allocate();

  typedef itk::VariableLengthVector<short> VariableVectorType;
  VariableVectorType variableLengthVector;
  variableLengthVector.SetSize(numOfComponents);
  for( size_t i = 0; i < numOfComponents; ++i )
    {
    if( i==0 || i==4 )
      {
      variableLengthVector[i] = 0; // assumed as b0 images
      }
    else
      {
      variableLengthVector[i] = i*10; // assumed as 6 gradient components
      }
    }

  nrrdVolume->FillBuffer(variableLengthVector);
  // nrrdVolume = CreateMetaData(nrrdVolume);
  return nrrdVolume;
}

/*
 ImageType::Pointer CreateMetaData(ImageType::Pointer nrrdVolume)
 {
 nrrdVolume->SetMetaDataDictionary(referenceDict);
 ...
 return nrrdVolume
 }
 */

int main(int , char * [])
{
  bool allTestPass = true;

  const size_t numOfComponents=8;
  // Create a vector image
  ImageType::Pointer nrrdVolume = CreateVolume(numOfComponents);
  //itk::MetaDataDictionary & referenceDict = nrrdVolume->GetMetaDataDictionary();
  DWIMetaDataDictionaryValidator   bldValidator;

  try {
    //spaceUnits testing
    {
    std::vector<std::string> tempSpaceUnits(3,std::string("mm"));
    bldValidator.SetSpaceUnits(tempSpaceUnits);
    const std::vector<std::string> outSpaceUnits = bldValidator.GetSpaceUnits();
    if(tempSpaceUnits != outSpaceUnits)
      {
      std::cout << "ERROR: SpaceUnits not preserved" << std::endl;
      for(size_t i = 0 ; i< outSpaceUnits.size(); ++i)
        {
        std::cout << "out outSpaceUnits " << outSpaceUnits[i] << std::endl;
        }
      allTestPass=false;
      }
    }

    //Centerings testing
    {
      std::vector<std::string> tempCenterings(4,std::string("cell"));
      tempCenterings[3] = "???";
      bldValidator.SetCenterings(tempCenterings);
      const std::vector<std::string> outCenterings = bldValidator.GetCenterings();
      if(tempCenterings != outCenterings)
        {
        std::cout << "ERROR: outCenterings not preserved" << std::endl;
        for(size_t i = 0 ; i< outCenterings.size(); ++i)
          {
          std::cout << "Out outCenterings " << outCenterings[i] << std::endl;
          }
        allTestPass=false;
        }
    }
    //Centerings testing
    {
    std::vector<std::string> tempCenterings(4,std::string("cell"));
    tempCenterings[3] = "???";
    bldValidator.SetCenterings(tempCenterings);
    const std::vector<std::string> outCenterings = bldValidator.GetCenterings();
    if(tempCenterings != outCenterings)
      {
      std::cout << "ERROR: outCenterings not preserved" << std::endl;
      for(size_t i = 0 ; i< outCenterings.size(); ++i)
        {
        std::cout << "Out outCenterings " << outCenterings[i] << std::endl;
        }
      allTestPass=false;
      }
    }
    //thickness testing
    {
    std::vector<double> tempThickness(4,std::nan(""));
      tempThickness[3] = 2.123;
      bldValidator.SetThicknesses(tempThickness);
      const std::vector<double> outThicknesses = bldValidator.GetThicknesses();
      if(tempThickness != outThicknesses)
        {
        std::cout << "ERROR: outThicknesses not preserved" << std::endl;
        for(size_t i = 0 ; i< outThicknesses.size(); ++i)
          {
          std::cout << "Out outThicknesses " << outThicknesses[i] << std::endl;
          }
        allTestPass=false;
        }
    }
    // Measurement Frame
    {
      std::vector<std::vector<double> > msrFrame(3);
      for( unsigned int saxi = 0; saxi < 3; saxi++ )
        {
        msrFrame[saxi].resize(3);
        for( unsigned int saxj = 0; saxj < 3; saxj++ )
          {
          msrFrame[saxi][saxj] = 0.0;
          }
        }
      msrFrame[0][0] = 1.0; msrFrame[1][1] = 1.0; msrFrame[2][2] = 1.0;
      bldValidator.SetMeasurementFrame(msrFrame);

      const std::vector<std::vector<double> > outMsr = bldValidator.GetMeasurementFrame();
      if(msrFrame != outMsr)
        {
        std::cout << "ERROR: outMsr not preserved" << std::endl;
        for(size_t i = 0 ; i< outMsr.size(); ++i)
          {
          for(size_t j=0; j< outMsr[i].size(); ++j)
          {
          std::cout << "Out outMsr " << i << " " << j << " " << outMsr[i][j] << std::endl;
          }
          }
        allTestPass=false;
        }
    }
    // Space
    {
    std::string tempSpace("left-posterior-superior"); //The only valid DWI modality
    bldValidator.SetSpace(tempSpace);
    const std::string outSpace = bldValidator.GetSpace();
    if(tempSpace != outSpace)
      {
      std::cout << "ERROR: outSpace not preserved" << std::endl;
      std::cout << "Out outSpace " << outSpace << std::endl;
      allTestPass=false;
      }
    }
    // Modality
    {
    std::string tempModality("DWMRI"); //The only valid DWI modality
    bldValidator.SetModality(tempModality);
    const std::string outModality = bldValidator.GetModality();
    if(tempModality != outModality)
      {
      std::cout << "ERROR: outModality not preserved" << std::endl;
      std::cout << "Out outModality " << outModality << std::endl;
      allTestPass=false;
      }
    }
    //B-Value
    {
    const double bValue=10000.0/3; //The only valid DWI modality
    bldValidator.SetBValue(bValue);
    const double outBvalue = bldValidator.GetBValue();
    if(bValue != outBvalue)
      {
      std::cout << "ERROR: outBvalue not preserved" << std::endl;
      std::cout << "Out outBvalue " << outBvalue << std::endl;
      allTestPass=false;
      }
    }
    //Gradient-Directions
    {
    /* We should apply direction vcl_cosines to gradient directions if requested by
     the user */
    std::vector<std::array<double, 3> > GradientTable(6);
    GradientTable[0][0] =1;
    GradientTable[0][1] =0;
    GradientTable[0][2] =0;

    GradientTable[1][0] =0;
    GradientTable[1][1] =1;
    GradientTable[1][2] =0;

    GradientTable[2][0] =0;
    GradientTable[2][1] =0;
    GradientTable[2][2] =1;

    GradientTable[3][0] =2;
    GradientTable[3][1] =0;
    GradientTable[3][2] =0;

    GradientTable[4][0] =0;
    GradientTable[4][1] =2;
    GradientTable[4][2] =0;


    GradientTable[5][0] =0;
    GradientTable[5][1] =0;
    GradientTable[5][2] =2;

    bldValidator.SetGradientTable(GradientTable);

    const std::vector<std::array<double, 3> > outGT=bldValidator.GetGradientTable();
    if(GradientTable != outGT)
      {
      std::cout << "ERROR: outGT not preserved" << std::endl;
      for(size_t i = 0 ; i< outGT.size(); ++i)
        {
        for(size_t j=0; j< outGT[i].size(); ++j)
          {
          std::cout << "Out outMsr " << i << " " << j << " " << outGT[i][j] << std::endl;
          }
        }
      allTestPass=false;
      }
    }
    PrintDictionaryHelper(bldValidator.GetMetaDataDictionary());


    std::cout << "\n\n\nTesting IO based copy from reference data\n" << std::endl;
    //Write DWI Image To Disk
    {
    // Add meta data to nrrd volume
    nrrdVolume->SetMetaDataDictionary(bldValidator.GetMetaDataDictionary());

    // Write Nrrd volume to disk
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseCompressionOn();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput( nrrdVolume );
    nrrdWriter->SetFileName( "./testNrrdImage.nrrd" );
    nrrdWriter->Update();
  }

    //Test Replicating DWI Image Image From Disk
    //  --Read Image
    //  --Replicate Image, Store Meta in Validator
    //  --Write Image
    //  --Read WrittenImage, Compare to previous metatdata Validator
    {
    const std::string referenceTestFileName="/Users/johnsonhj/src/NEP-11/BRAINSTools-build/ExternalData/TestData/DWI_TestData_OUTPUTS/PhilipsAchieva2.nrrd";

    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer nrrdReader = ReaderType::New();
    nrrdReader->SetFileName( referenceTestFileName );
    nrrdReader->Update();

    ImageType::Pointer refVolume =nrrdReader->GetOutput();

    DWIMetaDataDictionaryValidator refVolumeValidator;
    refVolumeValidator.SetMetaDataDictionary(refVolume->GetMetaDataDictionary());

    std::cout << "**** Reference Dictionary ****" << std::endl;
    PrintDictionaryHelper(refVolumeValidator.GetMetaDataDictionary());
    std::cout << "**** Reference Dictionary ****" << std::endl;

    //Now copy over the dictionaries elementby element to test storage/retrieval
    DWIMetaDataDictionaryValidator manualVolumeValidator;


    //spaceUnits testing
    manualVolumeValidator.SetSpaceUnits(refVolumeValidator.GetSpaceUnits());

    //Centerings testing
    manualVolumeValidator.SetCenterings(refVolumeValidator.GetCenterings());
    //Thicknesses
    manualVolumeValidator.SetThicknesses(refVolumeValidator.GetThicknesses());
    //Kinds testing
    manualVolumeValidator.SetKinds(refVolumeValidator.GetKinds());
    // Measurement Frame
    manualVolumeValidator.SetMeasurementFrame(refVolumeValidator.GetMeasurementFrame());
    // Modality
    manualVolumeValidator.SetModality(refVolumeValidator.GetModality());
    // Space
    manualVolumeValidator.SetSpace(refVolumeValidator.GetSpace());
    //B-Value
    manualVolumeValidator.SetBValue(refVolumeValidator.GetBValue());
    //Gradient-Directions
      {
      DWIMetaDataDictionaryValidator::GradientTableType temp = refVolumeValidator.GetGradientTable();
      manualVolumeValidator.SetGradientTable(temp);
      }

    std::cout << "=============================================================" << std::endl;
    std::cout << "**** Manual Dictionary ****" << std::endl;
    PrintDictionaryHelper(manualVolumeValidator.GetMetaDataDictionary());
    std::cout << "**** Manual Dictionary ****" << std::endl;

    //Now reset MetaDataDictionary from validator
    refVolume->SetMetaDataDictionary((manualVolumeValidator.GetMetaDataDictionary()));

    // Write Nrrd volume to disk
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseCompressionOn();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput( refVolume );
    nrrdWriter->SetFileName( "./outSameAsRefVolume.nrrd");
    nrrdWriter->Update();

    // Add meta data to nrrd volume
    nrrdVolume->SetMetaDataDictionary(bldValidator.GetMetaDataDictionary());
    }
  }
  catch (...)
  {
  throw;
  }

  if ( allTestPass )
    {
    std::cout << "SUCCESS!" << std::endl;
    return EXIT_SUCCESS;
    }
  std::cout << "FAILRE!" << std::endl;
  return EXIT_FAILURE;
}
