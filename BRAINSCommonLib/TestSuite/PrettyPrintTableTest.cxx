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
#include <iostream>
#include "PrettyPrintTable.h"
#include "BRAINSToolsVersion.h"

int main(int, char * *)
{
  PrettyPrintTable p;

  p.setTablePad(5);
  p.add(0, 0, "String");
  p.add(0, 1, "SecondColumn");
  p.add(0, 2, "4C");
  p.add(0, 3, "5C");
  p.add(1, 0, "Integers");
  p.add(1, 1, 1, "%d"); // Similar to %d
  p.add(1, 2, 2, "%d");
  p.add(1, 3, 3, "%d");
  p.add(1, 0, "ZeroPadInt");
  p.add(1, 1, 1, "%02d"); // Similar to %02d in printf
  p.add(1, 2, 2, "%02d");
  p.add(1, 3, 3, "%02d");
  p.add(2, 0, "FloatingPoint");
  p.add(2, 1, 1.0F, "%+5.2f"); // Similar to %5.2f in printf
  p.add(2, 2, 2.0F, "%+5.2f");
  p.add(2, 3, 3.0F, "%+5.2f");
  p.Print(std::cout);
  p.rightJustify();
  p.Print(std::cout);

  std::cout << "MajorVersion:  " << BRAINSTools::Version::MajorVersion() << std::endl;
  std::cout << "MinorVersion:  " << BRAINSTools::Version::MinorVersion() << std::endl;
  std::cout << "PatchVersion:  " << BRAINSTools::Version::PatchVersion() << std::endl;
  std::cout << "TweakVersion:  " << BRAINSTools::Version::TweakVersion() << std::endl;
  std::cout << "VersionString: " << BRAINSTools::Version::VersionString() << std::endl;
  std::cout << "BuildDate:     " << BRAINSTools::Version::BuildDate() << std::endl;

  std::cout << "ExtendedVersionString: " << BRAINSTools::Version::ExtendedVersionString() << std::endl;
  std::cout << "ToString:              " << BRAINSTools::Version::ToString() << std::endl;

  return EXIT_SUCCESS;
}
