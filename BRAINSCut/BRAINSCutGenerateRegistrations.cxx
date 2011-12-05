#include "BRAINSCutGenerateRegistrations.h"
#include "Utilities.h"

BRAINSCutGenerateRegistrations
::BRAINSCutGenerateRegistrations(  std::string netConfigurationFilename)
  : BRAINSCutPrimary( netConfigurationFilename )
{
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  SetAtlasDataSet();
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;
  SetAtlasFilename();
  std::cout << __LINE__ << "::" << __FILE__ << std::endl;

  SetRegistrationParametersFromNetConfiguration();
}

void
BRAINSCutGenerateRegistrations
::SetAtlasToSubjectRegistrationOn( bool onOff)
{
  atlasToSubjectRegistraionOn = onOff;
}

void
BRAINSCutGenerateRegistrations
::SetSubjectDataSet( bool applyDataSet )
{
  /** if applyDataSEt==false, then use training dataset */
  if( applyDataSet )
    {
    subjectDataSets = BRAINSCutNetConfiguration.GetApplyDataSets();
    }
  else
    {
    subjectDataSets = BRAINSCutNetConfiguration.GetTrainDataSets();
    }
}

void
BRAINSCutGenerateRegistrations
::GenerateRegistrations()
{
  for( std::list<DataSet *>::iterator subjectIt = subjectDataSets.begin();
       subjectIt != subjectDataSets.end();
       ++subjectIt )
    {
    const std::string subjectFilename( (*subjectIt)->GetImageFilenameByType(registrationImageTypeToUse) );

    const RegistrationType *subjectRegistration = (*subjectIt)->GetRegistrationWithID(registrationID);

    const std::string SubjectToAtlasRegistrationFilename
      ( subjectRegistration->GetAttribute<StringValue>("SubjToAtlasRegistrationFilename") );
    const std::string AtlasToSubjRegistrationFilename
      ( subjectRegistration->GetAttribute<StringValue>("AtlasToSubjRegistrationFilename") );
    const std::string AtlasBinaryFilename
      ( subjectRegistration->GetAttribute<StringValue>("AtlasBinaryFilename") );
    const std::string SubjectBinaryFilename
      ( subjectRegistration->GetAttribute<StringValue>("SubjectBinaryFilename") );

    if( atlasToSubjectRegistraionOn &&
        (!itksys::SystemTools::FileExists( SubjectToAtlasRegistrationFilename.c_str() ) ) )
      {
      CreateTransformFile(  subjectFilename,                    // moving image
                            atlasFilename,                      // fixed image
                            SubjectToAtlasRegistrationFilename,
                            AtlasBinaryFilename,                // fixed ROI
                            SubjectBinaryFilename,              // moving ROI
                            roiAutoDilateSize,
                            false );
      }
    else if( (!atlasToSubjectRegistraionOn) &&
             (!itksys::SystemTools::FileExists( AtlasToSubjRegistrationFilename.c_str() ) ) )
      {
      CreateTransformFile(  atlasFilename,                      // moving image
                            subjectFilename,                    // fixed image
                            AtlasToSubjRegistrationFilename,
                            SubjectBinaryFilename,              // fixed ROI
                            AtlasBinaryFilename,                // moving ROI
                            roiAutoDilateSize,
                            false );
      }
    }
}
