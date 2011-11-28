#include <NetConfiguration.h>
#include <itksys/SystemTools.hxx>
#include <itkMultiThreader.h>

#include <list>
#include <vector>
#include <iostream>
#include <itkIO.h>
#include "Utilities.h"
#include "NetConfigurationParser.h"

static int
WaitForProcess(itksysProcess * *processes, const unsigned processID)
{
  itksysProcess_WaitForExit(processes[processID], 0);
  int rval = itksysProcess_GetExitValue(processes[processID]);
  itksysProcess_Delete(processes[processID]);
  std::cerr << "CreateTransformFile returns " << rval << std::endl;
  processes[processID] = 0;
  return rval;
}

int GenerateRegistrations(NetConfiguration & ANNXMLObject,
                          bool reverse,
                          bool apply,
                          const unsigned int numThreads)
{
  std::cout << "GenerateRegistrations" << std::endl;
  // Get the atlas dataset
  DataSet *atlasDataSet = ANNXMLObject.GetAtlasDataSet();
  // find out the registration parameters
  RegistrationConfigurationParser *regParams =
    ANNXMLObject.Get<RegistrationConfigurationParser>("RegistrationConfiguration");
  const std::string imageTypeToUse
    ( regParams->GetAttribute<StringValue>(
      "ImageTypeToUse") );
  const std::string regID
    ( regParams->GetAttribute<StringValue>(
      "ID") );
  const double roiAUtoDilateSize
    ( regParams->GetAttribute<IntValue>("BRAINSROIAutoDilateSize") );

  // Get Atlas Image Name
  std::string AtlasImageFilename( atlasDataSet->GetImageFilenameByType(
                                    imageTypeToUse) );

  // We will use just original image without smoothing
  InternalImageType::SizeType radius;

  radius[0] = 0; radius[1] = 0; radius[2] = 0;
  // InternalImageType::Pointer AtlasImage = ReadMedianFilteredImage<InternalImageType>(AtlasImageFilename.c_str(),
  // radius);
  InternalImageType::Pointer AtlasImage =
    ReadMedianFilteredImage<InternalImageType>(AtlasImageFilename.c_str(
                                                 ), radius);

  // For each dataset, create deformation fields if they don't already exist.
  std::list<DataSet *> dataSets;
  if( apply )
    {
    dataSets = ANNXMLObject.GetApplyDataSets();
    }
  else
    {
    dataSets = ANNXMLObject.GetTrainDataSets();
    }
  // TODO Hard coded Max Processes == 9 ????
#define MAXPROCESSES 9
  itksysProcess *processes[MAXPROCESSES];
  for( unsigned i = 0; i < MAXPROCESSES; i++ )
    {
    processes[i] = 0;
    }

  unsigned process_count = 0;
  // this is a little squirrelly, so ...
  // you loop spawning processes until you fill all slots, setting
  // elements of processes non-zero. A non zero element means there's
  // a running process.
  // once all the slots are full, wait for a process to finish, and then
  // re-use the slot.
  // once all datasets have been gone through reap remaining processes.
  // three exit cases:
  // #of datasets < 4: unused slots = 0, which are ignored.
  // #of datasets % 4 == 0: 4 processes to wait for
  // #of datasets % 4 4 != 0: still 4 processes to wait for because
  // you re-use the slots one at a time, keeping them full until no
  // more datasets, then wait for all 4.
  for( std::list<DataSet *>::iterator it = dataSets.begin();
       it != dataSets.end(); ++it )
    {
    if( processes[process_count] != 0 )
      {
      int rval = WaitForProcess(processes, process_count);
      std::cerr << "CreateTransformFile returns " << rval << std::endl;
      }
    // get subject image
    std::string SubjectImage( ( *it )->GetImageFilenameByType(imageTypeToUse) );
    // Get SubjtoAtlasRegistrationFilename
    // RegistrationType *reg = (*it)->Get<RegistrationType>("Registration");
    const RegistrationType *reg = ( *it )->GetRegistrationWithID(regID);
    if( reg == NULL )
      {
      itkGenericExceptionMacro(<< "ERROR:  Invalid ID in RegistrationConfigurationParser.  ("
                               << reg
                               << ")"
                               << "with " << SubjectImage);
      }
    const std::string SubjToAtlasRegistrationFilename
      ( reg->GetAttribute<StringValue>(
        "SubjToAtlasRegistrationFilename") );
    const std::string AtlasToSubjRegistrationFilename
      ( reg->GetAttribute<StringValue>(
        "AtlasToSubjRegistrationFilename") );
    const std::string AtlasBinaryFilename
      ( reg->GetAttribute<StringValue>(
        "AtlasBinaryFilename") );
    const std::string SubjectBinaryFilename
      ( reg->GetAttribute<StringValue>(
        "SubjectBinaryFilename") );
#if 0
    const std::string landmarkType( reg->GetAttribute<StringValue>(
                                      "LandmarkType") );
    // Get Subject Landmark Name
    const std::string SubjectLandmark( ( *it )->GetAtlasFilenameByType(
                                         landmarkType) );
    const std::string AtlasLandmark( atlasDataSet->GetAtlasFilenameByType(
                                       landmarkType) );
#else
    const std::string landmarkType("");
    const std::string SubjectLandmark("");
    const std::string AtlasLandmark("");
#endif
    if( reverse == false )
      {
      if( !itksys::SystemTools::FileExists( SubjToAtlasRegistrationFilename.
                                            c_str() ) )
        {
        CreateTransformFile(SubjectImage,
                            AtlasImageFilename,
                            SubjToAtlasRegistrationFilename,
                            AtlasBinaryFilename,
                            SubjectBinaryFilename,
                            roiAUtoDilateSize,
                            true);
        ++process_count;
        }
      }
    else
      {
      if( !itksys::SystemTools::FileExists( AtlasToSubjRegistrationFilename.
                                            c_str() ) )
        {
        CreateTransformFile(AtlasImageFilename,
                            SubjectImage,
                            AtlasToSubjRegistrationFilename,
                            SubjectBinaryFilename,
                            AtlasBinaryFilename,
                            roiAUtoDilateSize,
                            true);
        ++process_count;
        }
      }
    if( process_count == MAXPROCESSES || process_count == numThreads )
      {
      process_count = 0;
      }
    }
  for( unsigned i = 0; i < MAXPROCESSES; i++ )
    {
    if( processes[i] != 0 )
      {
      int rval = WaitForProcess(processes, i);
      std::cerr << "CreateTransformFile returns " << rval << std::endl;
      }
    }
  return 0;
}
