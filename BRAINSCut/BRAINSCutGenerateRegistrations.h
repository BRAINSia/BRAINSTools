#ifndef BRAINSCutGenerateRegistrations_h
#define BRAINSCutGenerateRegistrations_h

#include "BRAINSCutPrimary.h"

typedef itk::Image<unsigned char, DIMENSION> BinaryImageType;
typedef BinaryImageType::Pointer             BinaryImagePointer;

class BRAINSCutGenerateRegistrations : public BRAINSCutPrimary
{
public:
  BRAINSCutGenerateRegistrations( std::string netConfigurationFilename);

  void SetAtlasToSubjectRegistrationOn(bool atalsToSubjectRegistration );

  void SetSubjectDataSet( bool applyDataSet );

  void GenerateRegistrations();

private:
  // std::string atlasImage @ BRAINSCutPrimary.h;
  std::string atlasBinaryFilename;

  bool atlasToSubjectRegistraionOn;

  std::list<DataSet *> subjectDataSets;

  /** private functions */

  void  CreateTransformFile(const std::string & MovingImageFilename, const std::string & FixedImageFilename,
                            const std::string & MovingBinaryImageFilename, const std::string & FixedBinaryImageFilename,
                            const std::string & OutputRegName, bool verbose);
};

#endif
