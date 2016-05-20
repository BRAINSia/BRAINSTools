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
#include "itkOutputWindow.h"
#include "itkTextOutput.h"
#include "itkTimeProbe.h"
#include "itkImageDuplicator.h"
#include <exception>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVersion.h"
#include "itksys/SystemTools.hxx"
#include "PrettyPrintTable.h"
#include "DenoiseFiltering.h"

#include "mu.h"
#include "EMSParameters.h"
#include "AtlasDefinition.h"
#include <vector>
#include "Log.h"

#include <itksys/SystemTools.hxx>

#include <iostream>
#include <string>
#include <sstream>
#include <map>

#include <cstdlib>

#include <StandardizeMaskIntensity.h>
#include "BRAINSABCCLP.h"
#include "BRAINSThreadControl.h"

// Use manually instantiated classes for the big program chunks
#define MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.h"
#include "AtlasRegistrationMethod.h"

//
// filter types
typedef AtlasRegistrationMethod<float, float> AtlasRegType;

typedef EMSegmentationFilter<AtlasRegType::OutputImageType,
                             AtlasRegType::ProbabilityImageType> SegFilterType;
//
// image types
typedef AtlasRegType::OutputImageType FloatImageType;
typedef AtlasRegType::ByteImageType   ByteImageType;
typedef itk::Image<short, 3>          ShortImageType;

typedef FloatImageType::Pointer FloatImagePointer;
typedef ByteImageType::Pointer  ByteImagePointer;
typedef ShortImageType::Pointer ShortImagePointer;

#undef MU_MANUAL_INSTANTIATION

/**
 * \def FindPathFromAtlasXML
 * \param The encoded file name, if abosolute, then use it, else prepend atlasDefinitionPath
 * \param atlasDefinitionPath the directory path for the XML file
 * \return either the absolute path, or the prepended path
 */
  static std::string
FindPathFromAtlasXML( const std::string xmlCodedPath, std::string atlasDefinitionPath)
{
  if( xmlCodedPath[0] != '/' ) // If not an absolute path, assume relative path
    {
    return atlasDefinitionPath+"/"+xmlCodedPath;
    }
  return xmlCodedPath;
}

//
// Get a version of the filename that does not include the preceeding path, or
// the image file extensions.
static std::string
GetStrippedImageFileNameExtension(const std::string & ImageFileName)
{
  std::vector<std::string> ExtensionsToRemove(6);

  ExtensionsToRemove[0] = ".gz";
  ExtensionsToRemove[1] = ".nii";
  ExtensionsToRemove[2] = ".hdr";
  ExtensionsToRemove[3] = ".img";
  ExtensionsToRemove[4] = ".gipl";
  ExtensionsToRemove[5] = ".nrrd";

  std::string returnString = itksys::SystemTools::GetFilenameName( ImageFileName );
  for(auto & elem : ExtensionsToRemove)
    {
    size_t rfind_location = returnString.rfind(elem);
    if( ( rfind_location != std::string::npos )
        && ( rfind_location == ( returnString.size() - elem.size() ) )
      )
      {
      returnString.replace(rfind_location, elem.size(), "");
      }
    }
  return returnString;
}


static AtlasRegType::MapOfFloatImageVectors
RescaleFunctionLocal( AtlasRegType::MapOfFloatImageVectors& localList)
{
  AtlasRegType::MapOfFloatImageVectors rval;

  for(auto & elem : localList)
    {
    for(auto imIt = elem.second.begin();
        imIt != elem.second.end(); ++imIt)
      {
      typedef itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>
        RescaleType;
      RescaleType::Pointer rescaler = RescaleType::New();
      rescaler->SetOutputMinimum(1);
      rescaler->SetOutputMaximum(MAX_IMAGE_OUTPUT_VALUE);
// #define INPLACE_RESCALER 0 // HACK Test this out
// #if defined( INPLACE_RESCALER )
//       rescaler->SetInPlace(true);
// #endif
      FloatImageType::Pointer tmp = (*imIt);
      rescaler->SetInput(tmp);
      rescaler->Update();
      rval[elem.first].push_back(rescaler->GetOutput());
// #if !defined( INPLACE_RESCALER)
//       (*imIt)  = rescaler->GetOutput();
// #endif
      }
    }
  return rval;
}


//
// utility method for constructing map of vectors
AtlasRegType::MapOfStringVectors
CreateTypedMap(const AtlasRegType::StringVector &keys, const AtlasRegType::StringVector &values)
{
  AtlasRegType::MapOfStringVectors rval;
  auto keyIt(keys.begin());
  auto valueIt(values.begin());
  for( ; keyIt != keys.end() && valueIt != values.end(); ++keyIt, ++valueIt)
    {
    rval[(*keyIt)].push_back((*valueIt));
    }
  return rval;
}

//For the BRAINSABC program, we also need to minimize num threads used by TBB
#include "tbb/task_scheduler_init.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  // Construct TBB task scheduler with matching threads to ITK threads
  tbb::task_scheduler_init init( itk::MultiThreader::GetGlobalDefaultNumberOfThreads() );

  // TODO:  Need to figure out how to conserve memory better during the running
  // of this application:  itk::DataObject::GlobalReleaseDataFlagOn();
  itk::OutputWindow::SetInstance( itk::TextOutput::New() );

  typedef itk::ImageFileReader<ByteImageType>                  ReaderType;
  typedef ReaderType::Pointer                                  ReaderPointer;


  // Check the parameters for valid values
  bool AllSimpleParameterChecksValid = true;
  if( maxIterations < 1 )
    {
    muLogMacro( << "Warning: --maxIterations set to 0, so only "
                << "initialization with priors will be completed." << std::endl );
    }
  if( input_Volumes.size() == 0 )
    {
    muLogMacro( <<  "ERROR: Must specify --inputVolumes" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( output_Volumes.size() != 1 &&
      output_Volumes.size() != input_Volumes.size() )
    {
    std::cerr << input_Volumes.size() << " images in input volumeslist, but "
              << output_Volumes.size() << " names in output volumes list"
              << "OR it must be exactly 1, and be the template for writing files."
              << std::endl;
    return EXIT_FAILURE;
    }
  if( input_VolumeTypes.size() != input_Volumes.size() )
    {
    muLogMacro( <<  "ERROR: "
                << "--inputVolumeTypes and --inputVolumes must"
                << " have the same number of elements" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( atlasDefinition == "" )
    {
    muLogMacro( <<  "Error: "
                << "--atlasDefinition <xml atlas def> required"
                << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( outputDir == "" )
    {
    muLogMacro( <<  "ERROR: "
                << "outputDir must be specified" << std::endl );
    AllSimpleParameterChecksValid = false;
    }
  if( AllSimpleParameterChecksValid == false )
    {
    muLogMacro( << "ERROR:  Commanline arguments are not valid." << std::endl );
    GENERATE_ECHOARGS;
    return EXIT_FAILURE;
    }
  const std::string atlasDefinitionPath = itksys::SystemTools::GetParentDirectory( atlasDefinition.c_str() );
  AtlasDefinition atlasDefinitionParser;
  try
    {
    atlasDefinitionParser.InitFromXML(atlasDefinition);
    }
  catch( ... )
    {
    muLogMacro( <<  "Error reading Atlas Definition from "
                << atlasDefinition
                << std::endl );
    return EXIT_FAILURE;
    }
  ;
  atlasDefinitionParser.DebugPrint();

  AtlasRegType::MapOfStringVectors inputVolumeMap =
    CreateTypedMap(input_VolumeTypes,input_Volumes);
  AtlasRegType::MapOfStringVectors outputVolumeMap;
  if(output_Volumes.size() > 1)
    {
    outputVolumeMap = CreateTypedMap(input_VolumeTypes,output_Volumes);
    }
  // Create and start a new timer (for the whole process)
  //  EMSTimer* timer = new EMSTimer();
  itk::TimeProbe timer;
  timer.Start();

  // Directory separator string
  std::string separator = std::string("/");
  //  separator[0] = MU_DIR_SEPARATOR; // always a '/' -- windows can handle
  // that fine...

  // Make sure last character in output directory string is a separator
  if( outputDir[outputDir.size() - 1] != '/' /* MU_DIR_SEPARATOR */ )
    {
    outputDir += separator;
    }

  // Create the output directory, stop if it does not exist
  // if(!mu::create_dir(outputDir.c_str()))
  if( !itksys::SystemTools::MakeDirectory( outputDir.c_str() ) )
    {
    muLogMacro( << "ERROR: Could not create requested output directory " << outputDir << std::endl );
    return EXIT_FAILURE;
    }

  // Set up the logger
  const std::string logfn = outputDir + defaultSuffix + ".log";
  ( mu::Log::GetInstance() )->EchoOn();
  ( mu::Log::GetInstance() )->SetOutputFileName( logfn.c_str() );

  // Set up suffix string for images
  std::string fmt = outputFormat;
  std::string outext = ".mha";
  if( itksys::SystemTools::Strucmp(fmt.c_str(), "Nrrd") == 0 )
    {
    outext = ".nrrd";
    }
  else if( itksys::SystemTools::Strucmp(fmt.c_str(), "Meta") == 0 )
    {
    outext = ".mha";
    }
  else if( itksys::SystemTools::Strucmp(fmt.c_str(), "NIFTI") == 0 )
    {
    outext = ".nii.gz";
    }
  else
    {
    muLogMacro(<< "WARNING: output format unrecognized, using Meta format\n");
    }

  const std::string suffstr
    = std::string("_") + std::string(defaultSuffix) + outext;

  muLogMacro(<< "mu::brainseg\n");
  muLogMacro(<< "========================================\n");
  muLogMacro(<< "Program compiled on: " << __DATE__ << std::endl );
  muLogMacro(<< std::endl );

  muLogMacro(<< "Hans J. Johnson - hans-johnson@uiowa.edu, has made significant"
             << " edits to this front end of the BRAINSABC system.\n");
  muLogMacro(<< "Original application was written by Marcel Prastawa - "
             << "prastawa@sci.utah.edu, and is maintained as a separate program.\n");
  muLogMacro(<< "This software is provided for research purposes only\n");
  muLogMacro(<< std::endl );

  muLogMacro(<< "Using ITK version "
             << itk::Version::GetITKMajorVersion() << "."
             << itk::Version::GetITKMinorVersion() << "."
             << itk::Version::GetITKBuildVersion() << std::endl );
  muLogMacro(<< std::endl );

  // Write input parameters
  muLogMacro(<< "=== Parameters ===\n");
  muLogMacro(<< "Suffix: " << defaultSuffix << std::endl );
  muLogMacro(<< "Output Directory: " << outputDir << std::endl );
  muLogMacro(<< "Output Format: " << outputFormat << std::endl );
  muLogMacro(<< "Input images: \n");
  muLogMacro(<< "Non-linear filtering, method: " << filterMethod << ", "
             << filterIteration
             << " iterations, dt = " << filterTimeStep << std::endl );

  const AtlasDefinition::TissueTypeVector & PriorNames
    = atlasDefinitionParser.TissueTypes();

  PrettyPrintTable AtlasDefTable;
  AtlasDefTable.add(0, 0, "Prior Names");
  AtlasDefTable.add(0, 1, ": [");
  AtlasDefTable.add(0, PriorNames.size() + 2 + 1, "]");
  SegFilterType::VectorType priorsWeightList;
  priorsWeightList.set_size( PriorNames.size() );
  AtlasDefTable.add(1, 0, "Prior weight scales");
  AtlasDefTable.add(1, 1, ": [");
  AtlasDefTable.add(1, PriorNames.size() + 2 + 1, "]");

  unsigned int currentRow = 0;
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    AtlasDefTable.add(currentRow, 2 + pwi, PriorNames[pwi]);
    }

  currentRow++;
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorsWeightList[pwi] = atlasDefinitionParser.GetWeight(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorsWeightList[pwi]);
    }

  currentRow++;
  SegFilterType::IntVectorType priorLabelCodeVector;
  priorLabelCodeVector.set_size( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior Label Codes");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorLabelCodeVector[pwi] = atlasDefinitionParser.GetLabelCode(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorLabelCodeVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::BoolVectorType priorIsForegroundPriorVector;
  priorIsForegroundPriorVector.resize( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior IsForeground");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorIsForegroundPriorVector[pwi] = atlasDefinitionParser.GetIsForegroundPrior(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorIsForegroundPriorVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::IntVectorType priorGaussianClusterCountVector;
  priorGaussianClusterCountVector.set_size( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior Clusters");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorGaussianClusterCountVector[pwi] = atlasDefinitionParser.GetGaussianClusterCount(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorGaussianClusterCountVector[pwi], "%d");
    }

  currentRow++;
  SegFilterType::BoolVectorType priorUseForBiasVector;
  priorUseForBiasVector.resize( PriorNames.size() );
  AtlasDefTable.add(currentRow, 0, "Prior For Bias");
  AtlasDefTable.add(currentRow, 1, ": [");
  AtlasDefTable.add(currentRow, PriorNames.size() + 2 + 1, "]");
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    priorUseForBiasVector[pwi] = atlasDefinitionParser.GetUseForBias(PriorNames[pwi]);
    AtlasDefTable.add(currentRow, 2 + pwi, priorUseForBiasVector[pwi], "%d");
    }

  { // Print out the ranges.
  currentRow++;
  for( unsigned int pwi = 0; pwi < PriorNames.size(); pwi++ )
    {
    orderedmap<std::string, AtlasDefinition::BoundsType> temp_range_List;
    for( unsigned int tt = 0; tt < input_VolumeTypes.size(); tt++ )
      {
      AtlasDefTable.add(currentRow + tt * 2 + 0, 0, std::string(input_VolumeTypes[tt]) + std::string(" Lower") );
      AtlasDefTable.add(currentRow + tt * 2 + 0, 1, ": [");
      AtlasDefTable.add(currentRow + tt * 2 + 0, PriorNames.size() + 2 + 1, "]");

      AtlasDefTable.add(currentRow + tt * 2 + 1, 0, std::string(input_VolumeTypes[tt]) + std::string(" Upper") );
      AtlasDefTable.add(currentRow + tt * 2 + 1, 1, ": [");
      AtlasDefTable.add(currentRow + tt * 2 + 1, PriorNames.size() + 2 + 1, "]");

      temp_range_List[input_VolumeTypes[tt]] = atlasDefinitionParser.GetBounds(PriorNames[pwi], input_VolumeTypes[tt]);
      AtlasDefTable.add(currentRow + tt * 2 + 0, 2 + pwi, temp_range_List[input_VolumeTypes[tt]].GetLower() );
      AtlasDefTable.add(currentRow + tt * 2 + 1, 2 + pwi, temp_range_List[input_VolumeTypes[tt]].GetUpper() );
      }
    }
  }

  {
  std::ostringstream oss;
  AtlasDefTable.Print(oss);
  muLogMacro( << oss.str() );
  }
  muLogMacro(<< "Max bias polynomial degree: " << maxBiasDegree << std::endl );
  muLogMacro(<< "Atlas warping: " << !atlasWarpingOff << std::endl );
  muLogMacro(<< "Atlas warp spline grid size: " << gridSize[0] << " X "
             << gridSize[1] << " X "
             << gridSize[2] << std::endl );
  muLogMacro(<< std::endl );
  muLogMacro(<< "=== Start Registration ===\n");

  GenericTransformType::Pointer atlasToSubjectPreSegmentationTransform = ITK_NULLPTR;

  AtlasRegType::MapOfFloatImageVectors atlasOriginalImageList;
  ByteImagePointer               atlasBrainMask;
  { // Read template images needed for atlas registration
  // muLogMacro(<< "Read template mask");
  const std::string templateMask = FindPathFromAtlasXML(atlasDefinitionParser.GetTemplateBrainMask(),atlasDefinitionPath);
  if( templateMask.size() < 1 )
    {
    muLogMacro( <<  "No template mask specified" << std::endl );
    return EXIT_FAILURE;
    }
  muLogMacro( << "Reading mask : " << templateMask << "...\n");

  ReaderPointer imgreader = ReaderType::New();
  imgreader->SetFileName( templateMask.c_str() );

  try
    {
    imgreader->Update();
    }
  catch( ... )
    {
    muLogMacro( << "ERROR:  Could not read image " << templateMask << "." << std::endl );
    return EXIT_FAILURE;
    }
  atlasBrainMask = imgreader->GetOutput();
  }

  AtlasRegType::MapOfFloatImageVectors intraSubjectRegisteredImageMap;
  AtlasRegType::MapOfFloatImageVectors intraSubjectRegisteredRawImageMap;
  std::vector<std::string>       priorfnlist;

  AtlasRegType::MapOfStringVectors templateVolumes;
  const AtlasDefinition::TemplateMap &tm = atlasDefinitionParser.GetTemplateVolumes();
  for(auto & elem : inputVolumeMap)
    {
    auto ti = tm.find(elem.first);
    std::string temp;
    if( ti != tm.end() )
      {
      temp = ti->second;
      std::cerr << "STATUS:  Atlas image of type: " << elem.first
                << " added with filename: " << temp << std::endl;
      }
    else if( elem.first.compare("IDWI") == 0 ) // do not throw when IDWI is an input type
      {
      std::cerr << "WARNING: No atlas image of type IDWI." << std::endl;
      }
    else
      {
      std::cerr << "ERROR:  Atlas image of type: " << elem.first
                << " not found in xml file." << std::endl;
      throw;
      }
    for(AtlasRegType::StringVector::const_iterator imIt = elem.second.begin();
        imIt != elem.second.end(); ++imIt)
      {
      templateVolumes[elem.first].push_back(temp);
      }
    }
  std::vector<bool> candidateDuplicatesList;
  {
  AtlasRegType::Pointer atlasreg = AtlasRegType::New();

  if( debuglevel > 0 )
    {
    atlasreg->DebugOn();
    atlasreg->SetDebugLevel(debuglevel);
    }

  atlasreg->SetSuffix(defaultSuffix);
  // Compute list of file names for the atlasOriginalPriors
  for(auto & PriorName : PriorNames)
    {
    priorfnlist.push_back( atlasDefinitionParser.GetPriorFilename( PriorName ) );
    }

  {
  AtlasRegType::MapOfFloatImageVectors intraSubjectRawImageMap;
  AtlasRegType::MapOfFloatImageVectors intraSubjectNoiseRemovedImageMap;


  { // StartOriginalImagesList
  const std::string suffixstr = "";
  { // Read subject images needed for atlas registration
  // muLogMacro(<< "Read subject images");
  if( input_Volumes.size() < 1 )
    {
    muLogMacro( <<  "No data images specified" << std::endl );
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<FloatImageType> LocalReaderType;
  typedef LocalReaderType::Pointer             LocalReaderPointer;
  FloatImageType::Pointer KeyImageFirstRead=ITK_NULLPTR;

  AtlasRegType::MapOfStringVectors intraSubjectTransformFileNames;
  unsigned int i = 0;
  for(auto typeIt = inputVolumeMap.begin();
      typeIt != inputVolumeMap.end(); ++typeIt)
    {
    for( AtlasRegType::StringVector::const_iterator imIt = typeIt->second.begin();
         imIt != typeIt->second.end(); ++imIt,++i )
      {
      muLogMacro(<< "\n***Reading image " << ": " << (*imIt) << "...\n");

      LocalReaderPointer imgreader = LocalReaderType::New();
      imgreader->SetFileName( (*imIt).c_str() );

      try
        {
        imgreader->Update();
        }
      catch( ... )
        {
        muLogMacro( << "ERROR:  Could not read image " << (*imIt) << "." << std::endl );
        return EXIT_FAILURE;
        }
      // Initialize with file read in
      FloatImagePointer typewiseEqualizedToFirstImage = imgreader->GetOutput();

      // Normalize Image Intensities:
      muLogMacro( << "Standardizing Intensities: ...\n" );
      FloatImagePointer curImage =
        StandardizeMaskIntensity<FloatImageType, ByteImageType>(typewiseEqualizedToFirstImage,
                                                                ITK_NULLPTR,
                                                                0.0005, 1.0 - 0.0005,
                                                                1, 0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                                0, MAX_IMAGE_OUTPUT_VALUE);
      intraSubjectRawImageMap[typeIt->first].push_back(curImage);
      muLogMacro( << "done.\n" );
      {
      std::vector<unsigned int> unused_gridSize;
      double                    localFilterTimeStep = filterTimeStep;
      if( localFilterTimeStep <= 0 )
        {
        FloatImageType::SpacingType::ValueType minPixelSize =
          std::numeric_limits<FloatImageType::SpacingType::ValueType>::max();
        const FloatImageType::SpacingType & imageSpacing = curImage->GetSpacing();
        for( FloatImageType::ImageDimensionType is = 0; is < FloatImageType::ImageDimension; ++is )
          {
          minPixelSize = std::min( minPixelSize, imageSpacing[is]);
          }
        localFilterTimeStep =
          ( (minPixelSize - std::numeric_limits<FloatImageType::SpacingType::ValueType>::epsilon() )
            / ( std::pow(2.0, static_cast<double>(FloatImageType::ImageDimension + 1) ) )
            );
        }
      FloatImagePointer denoisedImage =
        DenoiseFiltering<FloatImageType>(curImage, filterMethod, filterIteration,
                                         localFilterTimeStep, unused_gridSize);
      intraSubjectNoiseRemovedImageMap[typeIt->first].push_back(denoisedImage);
      if ( KeyImageFirstRead.IsNull() ) //The very first image nees to be the key image.
        {
        KeyImageFirstRead=denoisedImage;
        }

      if( debuglevel > 1 )
        {
        // DEBUG:  This code is for debugging purposes only;
        typedef itk::ImageFileWriter<FloatImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->UseCompressionOn();

        std::stringstream template_index_stream("");
        template_index_stream << i;
        const std::string fn = outputDir + "/DENOISED_INDEX_" + template_index_stream.str() + ".nii.gz";
        writer->SetInput(denoisedImage);
        writer->SetFileName(fn.c_str() );
        writer->Update();
        muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl);
        }
      }
      std::string intraSubjectTransformFileName  = outputDir
        + GetStrippedImageFileNameExtension((*imIt).c_str()) + "_to_"
        + GetStrippedImageFileNameExtension(GetMapVectorFirstElement(inputVolumeMap)) + suffixstr
        + ".h5";
      intraSubjectTransformFileNames[typeIt->first].push_back(intraSubjectTransformFileName);
      }
    }
  atlasreg->SetKeySubjectImage(KeyImageFirstRead);
  atlasreg->SetIntraSubjectOriginalImageList(intraSubjectNoiseRemovedImageMap);
  atlasreg->SetIntraSubjectTransformFileNames(intraSubjectTransformFileNames);
  }

  { // Read template images needed for atlas registration
  // muLogMacro(<< "Read template images");
  if( templateVolumes.empty() )
    {
    muLogMacro( <<  "No data images specified" << std::endl );
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<FloatImageType> LocalReaderType;
  typedef LocalReaderType::Pointer             LocalReaderPointer;


  for(auto mapIt = templateVolumes.begin();
      mapIt != templateVolumes.end(); ++mapIt)
    {
    const std::string curAtlasName = FindPathFromAtlasXML(*(mapIt->second.begin()),atlasDefinitionPath);
    muLogMacro(<< "\n***Reading atlas image " << mapIt->first << ": " << curAtlasName << "...\n");
    LocalReaderPointer imgreader = LocalReaderType::New();
    imgreader->SetFileName(curAtlasName.c_str());
    try
      {
      imgreader->Update();
      }
    catch( ... )
      {
      muLogMacro( << "ERROR:  Could not read image " << curAtlasName << "." << std::endl );
      return EXIT_FAILURE;
      }
    muLogMacro( << "Standardizing Intensities: ..." );
    FloatImagePointer img_i =
      StandardizeMaskIntensity<FloatImageType, ByteImageType>(imgreader->GetOutput(),
                                                              atlasBrainMask,
                                                              0.0005, 1.0 - 0.0005,
                                                              1,
                                                              0.95 * MAX_IMAGE_OUTPUT_VALUE,
                                                              0, MAX_IMAGE_OUTPUT_VALUE);
    muLogMacro( << "done." << std::endl );
    // the atlas pointers are all the same and parallel the input
    // image map of lists structure
    for(auto nameIt = mapIt->second.begin();
        nameIt != mapIt->second.end(); ++nameIt)
      {
      atlasOriginalImageList[mapIt->first].push_back(img_i);
      }
    if( debuglevel > 7 )
      {
      typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
      FloatWriterType::Pointer writer = FloatWriterType::New();

      std::stringstream write_atlas_index_stream("");
      const std::string fn = outputDir + "RenormalizedAtlasTemplate_" + mapIt->first  + suffstr;

      writer->SetInput(img_i);
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();
      writer->Update();
      }
    }

  atlasreg->SetAtlasOriginalImageList(atlasOriginalImageList);

  std::string atlasTransformFileName = "";
  if( atlasToSubjectTransform != "" )
    {
    // If this final transform filename exists,
    // it will be just read in and will be used directly
    // without doing registration.
    const bool doesExist = itksys::SystemTools::FileExists( atlasToSubjectTransform.c_str() );
    if( doesExist )
      {
      atlasTransformFileName = atlasToSubjectTransform;
      }
    else
      {
      atlasTransformFileName = GetStrippedImageFileNameExtension(atlasToSubjectTransform)
        + std::string("PreSegmentation.h5");
      }
    }
  else
    {
    atlasTransformFileName = outputDir
      + GetStrippedImageFileNameExtension(GetMapVectorFirstElement(templateVolumes))
      + std::string("_to_")
      + GetStrippedImageFileNameExtension(input_Volumes[0]) + suffixstr
      + std::string("PreSegmentation.h5");
    }

  atlasreg->SetAtlasToSubjectTransformFileName(atlasTransformFileName);
  }

  //    atlasreg->SetOutputDebugDir(outputDir);

  if( !( ( atlasToSubjectTransformType.compare("Identity") == 0 )
         || ( atlasToSubjectTransformType.compare("Rigid") == 0 )
         || ( atlasToSubjectTransformType.compare("Affine") == 0 )
         || ( atlasToSubjectTransformType.compare("BSpline") == 0 )
         || ( atlasToSubjectTransformType.compare("SyN") == 0 ) )
    )
    {
    muLogMacro("ERROR:  Invalid atlasToSubjectTransformType specified"
               << atlasToSubjectTransformType << std::endl);
    return EXIT_FAILURE;
    }

  if( !( ( subjectIntermodeTransformType.compare("Identity") == 0 )
         || ( subjectIntermodeTransformType.compare("Rigid") == 0 ) )
    )
    {
    muLogMacro("ERROR:  Invalid subjectIntermodeTransformType specified"
               << subjectIntermodeTransformType << std::endl
               << "Only Identity or Rigid transforms are allowed." << std::endl);
    return EXIT_FAILURE;
    }

  if( atlasToSubjectInitialTransform != "" )
    {
    muLogMacro(<< "atlasToSubjectInitialTransform specified." << std::endl);
    if( atlasToSubjectTransformType.compare("Identity") == 0 )
      {
      // Error because we're applying an identity transform by an initial transform was supplied.
      muLogMacro(<< "ERROR:  atlasToSubjectTransformType is Identity but an "
                 << "initial transform supplied." << std::endl);
      return EXIT_FAILURE;
      }
    GenericTransformType::Pointer atlasToSubjectCurrentGenericTransform;
    try
      {
      atlasToSubjectCurrentGenericTransform =
        itk::ReadTransformFromDisk(atlasToSubjectInitialTransform);
      }
    catch( itk::ExceptionObject & /* excp */ )
      {
      muLogMacro("ERROR:  Invalid atlasToSubjectInitialTransform specified"
                 << atlasToSubjectInitialTransform << std::endl);
      //muLogMacro( excp << std::endl);
      return EXIT_FAILURE;
      }

    const std::string initialTransformFileType = atlasToSubjectCurrentGenericTransform->GetNameOfClass();
    if( initialTransformFileType == "VersorRigid3DTransform"
        || initialTransformFileType == "ScaleVersor3DTransform"
        || initialTransformFileType == "ScaleSkewVersor3DTransform"
        || initialTransformFileType == "CompositeTransform")
      {
      if( !( (atlasToSubjectTransformType.compare("Rigid") == 0 )
             || ( atlasToSubjectTransformType.compare("Affine") == 0 )
             || ( atlasToSubjectTransformType.compare("BSpline") == 0 )
             || ( atlasToSubjectTransformType.compare("SyN") == 0 ) )
        )
        {
        muLogMacro(
          << "Error: initialAtlasToSubjectTransform is a "
          << initialTransformFileType
          << " but atlasToSubjectTransfromType is not Rigid, Affine, or BSpline."
          << std::endl);
        return EXIT_FAILURE;
        }
      }
    else if( initialTransformFileType == "AffineTransform" )
      {
      if( !( ( atlasToSubjectTransformType.compare("Affine") == 0 )
             || ( atlasToSubjectTransformType.compare("BSpline") == 0 )
             || ( atlasToSubjectTransformType.compare("SyN") == 0 ) )
        )
        {
        muLogMacro(<< "Error: initialAtlasToSubjectTransform is "
                   << "a AffineTransform but atlasToSubjectTransfromType "
                   << "is not Affine, or BSpline."
                   << std::endl);
        return EXIT_FAILURE;
        }
      }
    else if( initialTransformFileType == "BSplineTransform" )
      {
      if( !( ( atlasToSubjectTransformType.compare("BSpline") == 0 ) )
        )
        {
        muLogMacro(<< "Error: initialAtlasToSubjectTransform "
                   << "is a BSplineTransform but "
                   << "atlasToSubjectTransfromType is not BSpline."
                   << std::endl);
        return EXIT_FAILURE;
        }
      }
    else
      {
      muLogMacro( << "ERROR:  Invalid transform initializer type found:  "
                  << initialTransformFileType );
      return EXIT_FAILURE;
      }
    atlasreg->SetAtlasToSubjectInitialTransform(atlasToSubjectCurrentGenericTransform);
    }

  if( !restoreState.empty() )
    {
    typedef AtlasRegType::CompositeTransformType  CompositeTransformType;
    CompositeTransformType::Pointer restoreStateTransform;
    try
      {
      restoreStateTransform = dynamic_cast<CompositeTransformType *>( itk::ReadTransformFromDisk(restoreState).GetPointer() );
      if( restoreStateTransform.IsNull() )
        {
        muLogMacro("ERROR:  Could not restore transform type for file: "
                   << restoreState << std::endl);
        return EXIT_FAILURE;
        }
      }
    catch( itk::ExceptionObject & /* excp */ )
      {
      muLogMacro("ERROR:  Invalid restore state specified"
                 << restoreState << std::endl);
      return EXIT_FAILURE;
      }
    const std::string restoreStateFileType = restoreStateTransform->GetNameOfClass();
    if( restoreStateFileType != "CompositeTransform" )
      {
      muLogMacro("ERROR:  Invalid restore state file type. It must be a composite transform"
                 << restoreState << std::endl);
      return EXIT_FAILURE;
      }
    atlasreg->SetRestoreState(restoreStateTransform);
    }

  if( !saveState.empty() )
    {
    atlasreg->SetSaveState(saveState);
    }

  atlasreg->SetAtlasLinearTransformChoice(atlasToSubjectTransformType);
  atlasreg->SetImageLinearTransformChoice(subjectIntermodeTransformType);

  atlasreg->SetWarpGrid(gridSize[0], gridSize[1], gridSize[2]);
  muLogMacro(<< "Registering and resampling images..." << std::endl);
  itk::TimeProbe regtimer;
  regtimer.Start();
  atlasreg->SetOutputDebugDir(outputDir);
  try
    {
    atlasreg->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  regtimer.Stop();
  itk::RealTimeClock::TimeStampType elapsedTime
    = regtimer.GetTotal();
  muLogMacro(<< "Registration took " << elapsedTime
             << " " << regtimer.GetUnit() << std::endl);

  AtlasRegType::MapOfTransformLists intraSubjectTransforms =
    atlasreg->GetIntraSubjectTransforms();

  // ::ResampleImages()
  {
  // muLogMacro(<< "ResampleImages");
  intraSubjectRegisteredImageMap =
    ResampleInPlaceImageList(resamplerInterpolatorType, intraSubjectNoiseRemovedImageMap,
                             intraSubjectTransforms);

  intraSubjectRegisteredRawImageMap =
    ResampleInPlaceImageList(resamplerInterpolatorType, intraSubjectRawImageMap,
                             intraSubjectTransforms);
  //TODO: The maps size needs to be the same, but so do the lists within the maps.
  assert( intraSubjectRegisteredImageMap.size() == intraSubjectNoiseRemovedImageMap.size() );
  assert( intraSubjectRegisteredImageMap.size() == intraSubjectRawImageMap.size() );
  intraSubjectNoiseRemovedImageMap.clear();
  intraSubjectRawImageMap.clear();
  } // End registering data

  } // EndOriginalImagesList

  atlasToSubjectPreSegmentationTransform =
    atlasreg->GetAtlasToSubjectTransform();
  /*
  if( debuglevel > 9 )
    { // NOTE THIS IS REALLY ANNOYING FOR LARGE BSPLINE REGISTRATIONS!
    muLogMacro( << __FILE__ << " " << __LINE__ << " "
                << atlasToSubjectPreSegmentationTransform->GetFixedParameters() << std::endl );
    muLogMacro( << __FILE__ << " " << __LINE__ << " "
                << atlasToSubjectPreSegmentationTransform->GetParameters() << std::endl );
    }
  */
  if( debuglevel > 4 )
    {
    // Write the registered template and images
    if( !writeLess )
      {
      muLogMacro(<< "Writing registered template...\n");

      typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleType;
      typedef ResampleType::Pointer                                    ResamplePointer;
      ResamplePointer resampler = ResampleType::New();

      resampler->SetInput(GetMapVectorFirstElement(atlasOriginalImageList));
      resampler->SetTransform(atlasToSubjectPreSegmentationTransform);

      resampler->SetOutputParametersFromImage
        (*(intraSubjectRegisteredRawImageMap[atlasOriginalImageList.begin()->first].begin()));
      resampler->SetDefaultPixelValue(0);
      resampler->Update();
      typedef itk::CastImageFilter<FloatImageType, ShortImageType>
        ShortRescaleType;
      ShortRescaleType::Pointer rescaler = ShortRescaleType::New();
      rescaler->SetInput( resampler->GetOutput() );
      rescaler->Update();

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      std::string fn
        = outputDir + std::string("AtlasToSubjectInitialization_") + suffstr;

      writer->SetInput( rescaler->GetOutput() );
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();
      writer->Update();
      }
    for(auto & elem : intraSubjectRegisteredRawImageMap)
      {

      for(unsigned i = 0; i < elem.second.size(); ++i)
        {
        typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType>
          ShortRescaleType;

        ShortRescaleType::Pointer rescaler = ShortRescaleType::New();
        rescaler->SetOutputMinimum(0);
        rescaler->SetOutputMaximum(MAX_IMAGE_OUTPUT_VALUE);
        rescaler->SetInput(elem.second[i]);
        rescaler->Update();

        std::string fn = outputDir
          + GetStrippedImageFileNameExtension(inputVolumeMap[elem.first][i])
          + std::string("_registered") + suffstr;

        typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
        ShortWriterType::Pointer writer = ShortWriterType::New();

        writer->SetInput( rescaler->GetOutput() );
        writer->SetFileName( fn.c_str() );
        writer->UseCompressionOn();
        writer->Update();
        }
      }
    }
  }
  } // end atlas reg block

  std::cerr<< "Before RescaleFunctionLocal" << std::endl;
  PrintMapOfImageVectors(intraSubjectRegisteredImageMap);
  muLogMacro(<< "Rescale intensity of filtered images...\n");
  {
  intraSubjectRegisteredImageMap = RescaleFunctionLocal(intraSubjectRegisteredImageMap);
  intraSubjectRegisteredRawImageMap = RescaleFunctionLocal(intraSubjectRegisteredRawImageMap);
  }
  /*
   ------------------------------------
   *** BRAINSABC Input Images Flow: ***
   ------------------------------------

    {imageVolumeFileNames} -> [reader] -> {typewiseEqualizedToFirstImage} -> [StandardizedMaskIntensity]
              -> {currImage} --->   [DenoiseFiltering]  --->   {denoisedImage}
                      |                                               |
                      ~                                               ~
 push_back to:{intraSubjectRawImageMap}     push_back to:{intraSubjectNoiseRemovedImageMap} => <atlasreg>:{IntraSubjOriginalImageList}
                        |                                                       |
                        ~                                                       ~
     [ResampleInPlaceImageList] <--  <atlasreg>:{intraSubjTransforms} --> [ResampleInPlaceImageList]
                        |                                                       |
                        ~                                                       ~
  {intraSubjectRegisteredRawImageMap}  <->  [RescaleFunctionLocal]  <->  {intraSubjectRegisteredImageMap}
                        |                                                       |
                        ~                                                       ~
            set to <segfilter>:{RawInputImages}                     set to <segfilter>:{InputImages}

   - Notations:
   * class: <>
   * functional unit: []
   * file: {}
   ------------------------------------
   */

  ////////////////////////////////////
  // Start the segmentation process.
  ///////////////////////////////////
  bool usePurePlugs = false;
  if( purePlugsThreshold > 0 )
    {
    usePurePlugs = true;
    }

  for( unsigned int segmentationLevel = 0; segmentationLevel < 1; segmentationLevel++ )
    {
    SegFilterType::Pointer segfilter = SegFilterType::New();
    segfilter->SetUseKNN(useKNN);

    segfilter->SetUsePurePlugs(usePurePlugs);
    segfilter->SetPurePlugsThreshold(purePlugsThreshold);
    segfilter->SetNumberOfSubSamplesInEachPlugArea(
      numberOfSubSamplesInEachPlugArea[0],
      numberOfSubSamplesInEachPlugArea[1],
      numberOfSubSamplesInEachPlugArea[2]);

    std::vector<FloatImagePointer> atlasOriginalPriors( PriorNames.size() );
    unsigned int                   AirIndex = 10000;
    for( unsigned int i = 0; i < PriorNames.size(); i++ )
      {
      typedef itk::ImageFileReader<FloatImageType> LocalReaderType;
      typedef LocalReaderType::Pointer             LocalReaderPointer;
      LocalReaderPointer priorReader = LocalReaderType::New();
      const std::string curPriorAtlasName = FindPathFromAtlasXML(
        atlasDefinitionParser.GetPriorFilename(PriorNames[i]),
        atlasDefinitionPath);
      priorReader->SetFileName( curPriorAtlasName );
      priorReader->Update();
      FloatImageType::Pointer temp = priorReader->GetOutput();
      atlasOriginalPriors[i] = temp;
      // Set the index for the background values.
      if( PriorNames[i] == std::string("AIR") )
        {
        AirIndex = i;
        }
      }
    std::cout << "MESSAGE: USING AIR INDEX of :" << AirIndex << std::endl;
    segfilter->SetAirIndex(AirIndex);
    segfilter->SetPriors(atlasOriginalPriors);

    SegFilterType::RangeDBType myRanges;
    for(auto & PriorName : PriorNames)
      {
      orderedmap<std::string, AtlasDefinition::BoundsType> temp_range_List;
      for(auto & input_VolumeType : input_VolumeTypes)
        {
        temp_range_List[input_VolumeType] = atlasDefinitionParser.GetBounds(PriorName, input_VolumeType);
        }
      myRanges[PriorName] = temp_range_List;
      }
    segfilter->SetTissueTypeThresholdMapsRange(myRanges);
    segfilter->SetPriorNames(PriorNames);

    // //Static component that does not depend on priors-consolidation.
    muLogMacro(<< "Start segmentation...\n");
    segfilter->SetOutputDebugDir(outputDir);

    if( debuglevel > 0 )
      {
      segfilter->DebugOn();
      segfilter->SetDebugLevel(debuglevel);
      }

    PrintMapOfImageVectors(intraSubjectRegisteredImageMap);
    segfilter->SetInputImages(intraSubjectRegisteredImageMap);
    segfilter->SetRawInputImages(intraSubjectRegisteredRawImageMap);

    segfilter->SetMaximumIterations( maxIterations );
    segfilter->SetOriginalAtlasImages(atlasOriginalImageList);
    segfilter->SetTemplateBrainMask(atlasBrainMask);
    segfilter->SetTemplateGenericTransform(atlasToSubjectPreSegmentationTransform);

    segfilter->SetPriorWeights(priorsWeightList);
    segfilter->SetPriorLabelCodeVector(priorLabelCodeVector);
    segfilter->SetPriorUseForBiasVector(priorUseForBiasVector);
    segfilter->SetPriorIsForegroundPriorVector(priorIsForegroundPriorVector);

    segfilter->SetMaxBiasDegree(maxBiasDegree);
    // TODO: Expose the transform type to the BRAINSABC command line
    // segfilter->SetAtlasTransformType("SyN"); // atlasTransformType);

    if( !atlasWarpingOff )
      {
      // If the final transform filename exists, it will be read in and will be used directly,
      // so no update is run during iterations.
      if( atlasToSubjectTransform != "" && itksys::SystemTools::FileExists( atlasToSubjectTransform.c_str() ) )
        {
        segfilter->UpdateTransformationOff();
        }
      else
        {
        segfilter->UpdateTransformationOn();
        }
      }
    else
      {
      segfilter->UpdateTransformationOff();
      }
    segfilter->SetWarpGrid(
      gridSize[0],
      gridSize[1],
      gridSize[2]);

    segfilter->Update();

    // Write the secondary outputs
    if( !writeLess )
      {
      muLogMacro(<< "Writing filtered and bias corrected images...\n");
      // std::vector<FloatImagePointer> imgset = segfilter->GetCorrected();
      AtlasRegType::MapOfFloatImageVectors imgset = segfilter->GetRawCorrected();
      AtlasRegType::MapOfStringVectors outFileNames;

      if( output_Volumes.size() == 1 )
        {
        for(AtlasRegType::MapOfFloatImageVectors::const_iterator mapIt = imgset.begin();
            mapIt != imgset.end(); ++mapIt)
          {
          for( unsigned i = 0; i < mapIt->second.size(); i++ )
            {
            char buf[8192];
            sprintf(buf, output_Volumes[0].c_str(), mapIt->first.c_str(), i);
            outFileNames[mapIt->first].push_back(buf);
            }
          }
        }
      else if( TotalMapSize(imgset) != TotalMapSize(outputVolumeMap))
        {
        std::cerr << TotalMapSize(imgset) << " images in filter output, but "
                  << TotalMapSize(outputVolumeMap) << " names in output volumes list"
                  << std::endl;
        return EXIT_FAILURE;
        }
      else
        {
        outFileNames = outputVolumeMap;
        }
      for(AtlasRegType::MapOfFloatImageVectors::const_iterator mapIt = imgset.begin();
          mapIt != imgset.end(); ++mapIt)
        {
        for( unsigned i = 0; i < mapIt->second.size(); i++ )
          {
          // typedef itk::RescaleIntensityImageFilter<FloatImageType,
          // ShortImageType> ShortRescaleType;
          typedef itk::CastImageFilter<FloatImageType, ShortImageType> RescaleType;
          RescaleType::Pointer caster = RescaleType::New();

          caster->SetInput(mapIt->second[i]);
          caster->Update();
          // std::string fn
          //  = outputDir + GetStrippedImageFileNameExtension(names[i]) + std::string("_corrected")
          //    + suffstr;

          typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
          ShortWriterType::Pointer writer = ShortWriterType::New();

          writer->SetInput( caster->GetOutput() );
          writer->SetFileName(outFileNames[mapIt->first][i]);
          writer->UseCompressionOn();
          writer->Update();
          }
        }
      }

    AtlasRegType::MapOfFloatImageVectors imgset = segfilter->GetRawCorrected();
    // Average together all the input images of a given type
    for(auto & elem : imgset)
      {
      std::string volumeType = elem.first;
      // Can't average images of type other since it's really a mix of types.
      std::transform(volumeType.begin(), volumeType.end(), volumeType.begin(), tolower);
      if( volumeType == "other" )
        {
        continue;
        }

      FloatImagePointer avgImage = AverageImageList<FloatImageType>(elem.second);
      // Write out average image.
      typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType> RescaleType;
      RescaleType::Pointer rescaleInstensityFilter = RescaleType::New();
      rescaleInstensityFilter->SetOutputMinimum(0);
      rescaleInstensityFilter->SetOutputMaximum(4096);
      rescaleInstensityFilter->SetInput(avgImage);
      rescaleInstensityFilter->Update();

      std::string avgFileName = outputDir + volumeType + std::string("_average") + suffstr;

      muLogMacro(<< "Writing averaged corrected input images... " << avgFileName << std::endl );

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      writer->SetInput( rescaleInstensityFilter->GetOutput() );
      writer->SetFileName(avgFileName);
      writer->UseCompressionOn();
      writer->Update();
      }

    // Write warped template and bspline trafo
    if( !atlasWarpingOff )
      {
      SegFilterType::MapOfInputImageVectors WarpedAtlasList = segfilter->GenerateWarpedAtlasImages();
      for(auto & elem : WarpedAtlasList)
        {
        for( unsigned int index = 0; index < elem.second.size(); index++ )
          {
          typedef itk::RescaleIntensityImageFilter<FloatImageType, ByteImageType>
            ByteRescaleType;

          ByteRescaleType::Pointer rescaler = ByteRescaleType::New();
          rescaler->SetOutputMinimum(0);
          rescaler->SetOutputMaximum(255);
          rescaler->SetInput(elem.second[index]);
          rescaler->Update();

          typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
          ByteWriterType::Pointer writer = ByteWriterType::New();

          const std::string fn = outputDir
            + GetStrippedImageFileNameExtension(templateVolumes[elem.first][index]) + std::string("_to_")
            + GetStrippedImageFileNameExtension( ( inputVolumeMap[elem.first][index] ) )
            + std::string("_warped") + std::string(".nii.gz");

          muLogMacro(<< "Writing warped template images... " << fn << std::endl );
          writer->SetInput( rescaler->GetOutput() );
          writer->SetFileName( fn.c_str() );
          writer->UseCompressionOn();
          writer->Update();
          }
        }

      // If this final transform filename exists, it has been used directly
      // with no pre-segmentation registration and no updating during iterations,
      // so there is no need to overwrite that.
      if( atlasToSubjectTransform != "" && !itksys::SystemTools::FileExists( atlasToSubjectTransform.c_str() ) )
        {
        const std::string postSegmentationTransformFileName = atlasToSubjectTransform;
        // NOTE:  Aliasing of smart-pointers up the polymorphic tree OK here
        // because
        // the primary
        // smart pointer is gaunteed to exists during lifespan of all aliases.
        muLogMacro(<< "Writing final atlas to subject template... " << postSegmentationTransformFileName << std::endl );
        GenericTransformType::Pointer atlasToSubjectPostSegmentationTransform =
          segfilter->GetTemplateGenericTransform();
        itk::WriteTransformToDisk<double,float>(atlasToSubjectPostSegmentationTransform, postSegmentationTransformFileName);
        }
      }

    // Write the labels
    muLogMacro(<< "Writing labels...\n");
    {
    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    writer->SetInput( segfilter->GetOutput() );
    std::string fn;
    if( outputLabels == "" )
      {
      fn = outputDir;
      fn += GetStrippedImageFileNameExtension(inputVolumeMap.begin()->second[0]);
      fn += std::string("_DirtyLabels");
      fn += suffstr;
      }
    else
      {
      fn = outputDirtyLabels;
      }
    writer->SetFileName( fn.c_str() );
    writer->UseCompressionOn();
    writer->Update();
    }
    {
    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    std::string fn;
    if( outputLabels == "" )
      {
      fn = outputDir;
      fn += GetStrippedImageFileNameExtension(inputVolumeMap.begin()->second[0]);
      fn += std::string("_labels");
      fn += suffstr;
      }
    else
      {
      fn = outputLabels;
      }

    writer->SetInput( segfilter->GetCleanedOutput() );
    writer->SetFileName( fn.c_str() );
    writer->UseCompressionOn();

    try
      {
      writer->Update();

      fn = outputDir;
      fn += "thresholded_labels.nii.gz";
      writer->SetInput( segfilter->GetThresholdedOutput() );
      writer->SetFileName( fn );
      writer->UseCompressionOn();
      writer->Modified();
      writer->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      muLogMacro( <<  e << std::endl );
      return -1;
      }
    catch( std::exception & e )
      {
      muLogMacro( <<  "Exception: " << e.what() << std::endl );
      return -1;
      }
    catch( std::string & s )
      {
      muLogMacro( <<  "Exception: " << s << std::endl );
      return -1;
      }
    catch( ... )
      {
      muLogMacro( <<  "Unknown exception" << std::endl );
      muLogMacro( << "failed to write image " << fn << std::endl );
      return -1;
      }
    }
    // Write final Posteriors
    // NOTE :  Priors and Posteriors should correspond, so use the PriorNames
    // to get the names.
    for( unsigned int probabilityIndex = 0; probabilityIndex < PriorNames.size(); probabilityIndex++ )
      {
      std::string fn;
      if( posteriorTemplate == "" )
        {
        fn = outputDir;
        fn += GetStrippedImageFileNameExtension(inputVolumeMap.begin()->second[0]);
        fn += "_POSTERIOR_";
        fn += PriorNames[probabilityIndex];
        fn += suffstr;
        }
      else
        {
        char buf[8192];
        sprintf(buf, posteriorTemplate.c_str(),
                PriorNames[probabilityIndex].c_str() );
        fn = buf;
        }
      typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
      FloatWriterType::Pointer writer = FloatWriterType::New();

      FloatImageType::Pointer currPosterior = segfilter->GetPosteriors()[probabilityIndex];
      writer->SetInput( currPosterior );
      writer->SetFileName( fn.c_str() );
      writer->UseCompressionOn();
      writer->Update();
      }
    }
  timer.Stop();
  muLogMacro(<< "All segmentation processes took " << timer.GetTotal() << " " << timer.GetUnit() << std::endl );
  return EXIT_SUCCESS;
}
