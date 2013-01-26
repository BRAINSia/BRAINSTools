#ifndef BRAINSCutGenerateRegistrations_h
#define BRAINSCutGenerateRegistrations_h

#include "BRAINSCutDataHandler.h"

typedef itk::Image<unsigned char, DIMENSION> BinaryImageType;
typedef BinaryImageType::Pointer             BinaryImagePointer;

class BRAINSCutGenerateRegistrations
{
public:
  BRAINSCutGenerateRegistrations( BRAINSCutDataHandler& dataHandler );

  void SetAtlasToSubjectRegistrationOn(bool atalsToSubjectRegistration );

  void SetDataSet( bool applyDataSet );

  void GenerateRegistrations();

private:
  BRAINSCutDataHandler* myDataHandler;
  bool                  atlasToSubjectRegistraionOn;

  std::list<DataSet *> subjectDataSets;

  /** private functions */

  void  CreateTransformFile(const std::string & MovingImageFilename, const std::string & FixedImageFilename,
                            const std::string & MovingBinaryImageFilename, const std::string & FixedBinaryImageFilename,
                            const std::string & OutputRegName, bool verbose);
};

#endif
