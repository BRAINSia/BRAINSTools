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
#ifndef RegistrationParams_h
#define RegistrationParams_h
#include "CompoundObjectBase.h"
#include "StringValue.h"

class RegistrationParams : public CompoundObjectBase
{
public:
  using SuperClass = CompoundObjectBase;
  virtual int
  PrintSelf( std::ostream & os, int indent ) const
  {
    indent += SuperClass::PrintSelf( os, indent );
    os << this->PrintSpaces( indent ) << "=== RegistrationParams ===" << std::endl;
    return indent + 2;
  }

  RegistrationParams()
    : CompoundObjectBase( "RegistrationParams" )
  {
    this->Add( new StringValue( "ImageTypeToUse", "" ), "ImageTypeToUse" );
    this->Add( new StringValue( "ID", "" ), "ID" );
  }
};

#endif // RegistrationParams_h
