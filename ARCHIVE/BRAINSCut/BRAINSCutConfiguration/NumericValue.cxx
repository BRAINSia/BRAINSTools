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
#include "FloatValue.h"
#include "IntValue.h"
#include <cstdlib>
#include "BRAINSCutExceptionStringHandler.h"

//
// for each type define Verify, and SetValue (from string)

// eventually set up domain/range checking?
bool
IntValue::Verify() const
{
  return true;
}

bool
FloatValue::Verify() const
{
  return true;
}

void
IntValue::SetValue(const std::string & stringval)
{
  char *test;
  long  val = strtol(stringval.c_str(), &test, 10);

  if( test == stringval.c_str() )
    {
    std::string msg("Can't convert *");
    msg += stringval;
    msg += ") to integer";
    throw BRAINSCutExceptionStringHandler(msg);
    }
  this->m_Value = val;
}

void
FloatValue::SetValue(const std::string & stringval)
{
  char * test;
  double val = strtod(stringval.c_str(), &test);

  if( test == stringval.c_str() )
    {
    std::string msg("Can't convert *");
    msg += stringval;
    msg += ") to float";
    throw BRAINSCutExceptionStringHandler(msg);
    }
  this->NumericValue<double>::SetValue(val);
}
