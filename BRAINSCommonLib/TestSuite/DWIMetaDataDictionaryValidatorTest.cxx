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
#include <itkImageRegionIteratorWithIndex.h>

#include "DWIMetaDataDictionaryValidator.h"

//Force printing of key values nicely
static std::string ForceConvert( const itk::MetaDataObjectBase *  myMetaDataObjectBase)
{
  const itk::MetaDataObject<std::string> * myMetaDataObject = dynamic_cast<const itk::MetaDataObject<std::string> * >(myMetaDataObjectBase);
  if(myMetaDataObject)
    {
    //const std::string temp  =myMetaDataObject->GetMetaDataObjectTypeName();
    const std::string temp("std::string");
    return myMetaDataObject->GetMetaDataObjectValue()+"         (type): " + temp; //assuming you define print for matrix
    }
  else
    {
    // double type for thickness field
    const itk::MetaDataObject<double> * doubleMetaDataObject = dynamic_cast<const itk::MetaDataObject<double> * >(myMetaDataObjectBase);
    if(doubleMetaDataObject)
      {
      const std::string temp("double");
      // convert double value to string
      std::ostringstream strs;
      strs << doubleMetaDataObject->GetMetaDataObjectValue();
      std::string str = strs.str();

      return str+"         (type): " + temp;
      }
    else
      {
      const std::string temp  = myMetaDataObjectBase->GetMetaDataObjectTypeName();
      return "???_conversion         (type): " + temp;
      }
    }
}

static void PrintDictionaryHelper(const itk::MetaDataDictionary & dictPrint)
{
  std::cout << "----------------" << std::endl;
  itk::MetaDataDictionary::ConstIterator end = dictPrint.End();
  for ( itk::MetaDataDictionary::ConstIterator it = dictPrint.Begin(); it != end; ++it)
    {
    if( it->first.find("NRRD_measurement frame") != std::string::npos )
      {
      std::cout << ' ' << it->first << ":=" << std::endl;
      typedef std::vector<std::vector<double> > msrFrameType;
      const itk::MetaDataObject<msrFrameType> * msrFrameMetaDataObject =
                                                    dynamic_cast<const itk::MetaDataObject<msrFrameType> * >(it->second.GetPointer());
      const msrFrameType outMsr = msrFrameMetaDataObject->GetMetaDataObjectValue();
      for(size_t i = 0 ; i< outMsr.size(); ++i)
        {
        std::cout << "  ";
        for(size_t j=0; j< outMsr[i].size(); ++j)
          {
          std::cout << outMsr[i][j] << " ";
          }
        std::cout << std::endl;
        }
      }
    else
      {
      std::cout << ' ' << it->first << ":=" << ForceConvert( it->second ) << std::endl;
      }
    }
  std::cout << "----------------" << std::endl;
}

// Create a vector image
typedef short                           PixelType;
typedef itk::VectorImage<PixelType, 3>  VectorImageType;

static VectorImageType::Pointer CreateVolume(const size_t numOfComponents)
{
  const int imageSize = 11; // each image component has size of imageSize^3

  VectorImageType::IndexType start;
  start.Fill(0);

  VectorImageType::SizeType size;
  size.Fill(imageSize);

  VectorImageType::RegionType region(start, size);

  VectorImageType::Pointer nrrdVolume = VectorImageType::New();
  nrrdVolume->SetRegions(region);
  nrrdVolume->SetVectorLength(numOfComponents);
  nrrdVolume->Allocate();

  itk::VariableLengthVector< PixelType > ZeroPixel( numOfComponents );
  ZeroPixel.Fill( itk::NumericTraits< PixelType >::Zero );
  nrrdVolume->FillBuffer(ZeroPixel);

  itk::VariableLengthVector<PixelType> f( numOfComponents );
  for( size_t i = 0; i < numOfComponents; ++i )
    {
    if( i==0 || i==4 )
      {
      f[i] = 255; // assumed as b0 images
      }
    else
      {
      f[i] = i*10; // assumed as 6 gradient components
      }
    }
  // define a sub region
  start.Fill(3);
  size.Fill(4);
  VectorImageType::RegionType subRegion( start, size );
  typedef itk::ImageRegionIterator< VectorImageType > IteratorType;
  IteratorType it( nrrdVolume, subRegion );
  it.GoToBegin();
  while( !it.IsAtEnd() )
    {
    it.Set( f );
    ++it;
    }

  return nrrdVolume;
}

// TEST PROGRAM

int main( int argc, char *argv[] )
{
  if( argc < 4)
    {
    std::cout << argv[0] << "outputArtificialTestNrrdImage inputReferenceImage outputReplicatedReferenceImage"
              << std::endl;
    return EXIT_FAILURE;
    }

  bool allTestPass = true;

  // TEST #1
  /*
   * FIRST TEST:
   *  - Sets individual fields in the validator to create an output MetaDataDictionary
   *  - Creates a vector image with default values
   *  - Writes the created volume with output MetaDataDictionary
   */

  const size_t numOfComponents=8;

  // Create a vector image
  VectorImageType::Pointer nrrdVolume = CreateVolume(numOfComponents);

  // Instantiate a validator object
  DWIMetaDataDictionaryValidator   bldValidator;

  // Set and test validator fields individually
  try
  {
    /*
     NOTE: "centerings" and "thickness" should be created based on "volume interleaved".
           If the image is "pixel interleave" (like vectorImage in ITK), the NrrdIO
           will automatically handle the correct permutation.
    */
    //Centerings testing
    {
    std::vector<std::string> tempCenterings(4,std::string("cell"));
    tempCenterings[3] = "???";
    bldValidator.SetCenterings(tempCenterings);
    const std::vector<std::string> outCenterings = bldValidator.GetCenterings();
    if(tempCenterings != outCenterings)
      {
      std::cout << "ERROR: outCenterings not preserved" << std::endl;
      for(size_t i = 0 ; i < outCenterings.size(); ++i)
        {
        std::cout << "Out outCenterings " << outCenterings[i] << std::endl;
        }
      allTestPass=false;
      }
    }

    //thickness testing
    {
    bool thicknessPass = true;
    std::vector<double> tempThickness(4,std::nan(""));
    tempThickness[2] = 2.123;
    bldValidator.SetThicknesses(tempThickness);
    const std::vector<double> outThicknesses = bldValidator.GetThicknesses();
    if(tempThickness.size() != outThicknesses.size())
      {
      thicknessPass = false;
      }
    else
      {
      for(size_t i = 0 ; i < outThicknesses.size(); ++i)
         {
         if( std::isnan(outThicknesses[i]) )
           {
           if( !std::isnan(tempThickness[i]) )
             {
             thicknessPass = false;
             }
           }
         else
           {
           if( outThicknesses[i] != tempThickness[i] )
             {
             thicknessPass = false;
             }
           }
         }
      }
    if( thicknessPass == false )
      {
      std::cout << "ERROR: outThicknesses not preserved" << std::endl;
      for(size_t i = 0 ; i < outThicknesses.size(); ++i)
        {
        std::cout << "Output Thicknesses " << outThicknesses[i] << std::endl;
        }
      allTestPass = false;
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
    std::vector<std::array<double, 3> > GradientTable( numOfComponents );
    GradientTable[0][0] =0; //first b0 image
    GradientTable[0][1] =0;
    GradientTable[0][2] =0;

    GradientTable[1][0] =1;
    GradientTable[1][1] =0;
    GradientTable[1][2] =0;

    GradientTable[2][0] =0;
    GradientTable[2][1] =1;
    GradientTable[2][2] =0;

    GradientTable[3][0] =0;
    GradientTable[3][1] =0;
    GradientTable[3][2] =1;

    GradientTable[4][0] =0; //second b0 image
    GradientTable[4][1] =0;
    GradientTable[4][2] =0;

    GradientTable[5][0] =2;
    GradientTable[5][1] =0;
    GradientTable[5][2] =0;

    GradientTable[6][0] =0;
    GradientTable[6][1] =2;
    GradientTable[6][2] =0;

    GradientTable[7][0] =0;
    GradientTable[7][1] =0;
    GradientTable[7][2] =2;

    bldValidator.SetGradientTable(GradientTable);

    const std::vector<std::array<double, 3> > outGT = bldValidator.GetGradientTable();
    if(GradientTable != outGT)
      {
      std::cout << "ERROR: outGT not preserved! Output outGT:" << std::endl;
      for(size_t i = 0 ; i< outGT.size(); ++i)
        {
        for(size_t j=0; j< outGT[i].size(); ++j)
          {
          std::cout << outGT[i][j] << " ";
          }
        std::cout << std::endl;
        }
      allTestPass=false;
      }
    }
    std::cout << "\n****\\begin Artificial Dictionary ****" << std::endl;
    PrintDictionaryHelper(bldValidator.GetMetaDataDictionary());
    std::cout << "****\\end Artificial Dictionary ****" << std::endl;

    std::cout << "\nWrite the artificial vector image to the disk using the created artificial MetaDataDictionary...\n" << std::endl;
    //Write DWI Image To Disk
    {
    // Add meta data to nrrd volume
    nrrdVolume->SetMetaDataDictionary(bldValidator.GetMetaDataDictionary());

    // Write Nrrd volume to disk
    typedef itk::ImageFileWriter<VectorImageType> WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseCompressionOn();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput( nrrdVolume );
    nrrdWriter->SetFileName( argv[1] );
    nrrdWriter->Update();
    }

    // TEST #2
    /*
     * SECOND TEST:
     *  - Read a reference DWI image as a vector image in ITK
     *  - Set individual fields of a manual validator from the reference image metaDataDictionary
     *  - Writes the itk reference image with the manual MetaDataDictionary
     *  - Compares the input reference image with the replicated reference volume
     */

    std::cout << "\n\n\n>>>Testing IO based copy from reference data:\n" << std::endl;
    //Test Replicating DWI Image Image From Disk
    //  --Read Reference Image
    //  --Replicate its meta data dictionary using the validator
    //  --Write replicated reference image using the new meta data
    //  --Compare the replicated image with input reference volume
    {
    typedef itk::ImageFileReader<VectorImageType> ReaderType;
    ReaderType::Pointer nrrdReader = ReaderType::New();
    nrrdReader->SetFileName( argv[2] );
    nrrdReader->Update();

    VectorImageType::Pointer refVolume =nrrdReader->GetOutput();

    // Create a reference validator by setting its MetaDataDictionary from the input reference volume
    DWIMetaDataDictionaryValidator refVolumeValidator;
    refVolumeValidator.SetMetaDataDictionary(refVolume->GetMetaDataDictionary());

    std::cout << "****\\begin Reference Dictionary ****" << std::endl;
    PrintDictionaryHelper(refVolumeValidator.GetMetaDataDictionary());
    std::cout << "****\\end Reference Dictionary ****" << std::endl;

    // Now copy over the dictionaries element by element to test storage/retrieval
    DWIMetaDataDictionaryValidator manualVolumeValidator;

    /* Fields that need to be set:
       - thickness
       - centerings
       - measurement frame
       - modality
       - b-value
       - gradients
     */
    //Centerings testing
    manualVolumeValidator.SetCenterings(refVolumeValidator.GetCenterings());
    //Thicknesses
    manualVolumeValidator.SetThicknesses(refVolumeValidator.GetThicknesses());
    // Measurement Frame
    manualVolumeValidator.SetMeasurementFrame(refVolumeValidator.GetMeasurementFrame());
    // Modality
    manualVolumeValidator.SetModality(refVolumeValidator.GetModality());
    //B-Value
    manualVolumeValidator.SetBValue(refVolumeValidator.GetBValue());
    //Gradient-Directions
    {
    DWIMetaDataDictionaryValidator::GradientTableType temp = refVolumeValidator.GetGradientTable();
    manualVolumeValidator.SetGradientTable(temp);
    }

    std::cout << "\n=============================================================\n" << std::endl;
    std::cout << "****\\begin Manual Dictionary ****" << std::endl;
    PrintDictionaryHelper(manualVolumeValidator.GetMetaDataDictionary());
    std::cout << "****\\end Manual Dictionary ****" << std::endl;

    //Now reset MetaDataDictionary from validator
    refVolume->SetMetaDataDictionary(manualVolumeValidator.GetMetaDataDictionary());

    // Write Nrrd volume to disk
    typedef itk::ImageFileWriter<VectorImageType> WriterType;
    WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseCompressionOn();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput( refVolume );
    nrrdWriter->SetFileName( argv[3] );
    nrrdWriter->Update();
    }
  }
  catch(...)
  {
    throw;
  }

  if( allTestPass )
    {
    std::cout << "SUCCESS!" << std::endl;
    return EXIT_SUCCESS;
    }
  std::cout << "FAILURE!" << std::endl;
  return EXIT_FAILURE;
}
