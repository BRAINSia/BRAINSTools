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
#ifndef ProbabilityMap_h
#define ProbabilityMap_h
#include "StringValue.h"
#include "FloatValue.h"
#include "CompoundObjectBase.h"

class ProbabilityMap : public CompoundObjectBase
{
public:
  using SuperClass = CompoundObjectBase;
  virtual int
  PrintSelf( std::ostream & os, int indent ) const
  {
    indent += SuperClass::PrintSelf( os, indent );
    os << this->PrintSpaces( indent ) << "=== ProbabilityMap ===" << std::endl;
    return indent + 2;
  }

  ProbabilityMap()
    : CompoundObjectBase( "ProbabilityMap" )
  {
    this->Add( new StringValue( "StructureID", "" ), "StructureID" );
    this->Add( new StringValue( "GenerateVector", "" ), "GenerateVector" );
    this->Add( new FloatValue( "Gaussian", 1.0 ), "Gaussian" );
    this->Add( new StringValue( "Filename", "" ), "Filename" );
    this->Add( new StringValue( "rho", "" ), "rho" );
    this->Add( new StringValue( "phi", "" ), "phi" );
    this->Add( new StringValue( "theta", "" ), "theta" );
  }
};

class ProbabilityMapList : public CompoundObjectBase
{
public:
  ProbabilityMapList()
    : CompoundObjectBase( "ProbabilityMapList" )
  {}
};

#endif
