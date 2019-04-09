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
#ifndef BRAINSCutGenerateRegistrations_h
#define BRAINSCutGenerateRegistrations_h

#include "BRAINSCutDataHandler.h"

using BinaryImageType = itk::Image<unsigned char, DIMENSION>;
using BinaryImagePointer = BinaryImageType::Pointer;

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
  std::list<DataSet *>  subjectDataSets;

  /** private functions */

  void  CreateTransformFile(const std::string & MovingImageFilename, const std::string & FixedImageFilename,
                            const std::string & MovingBinaryImageFilename, const std::string & FixedBinaryImageFilename,
                            const std::string & OutputRegName, bool verbose);
};

#endif
