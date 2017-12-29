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
#ifndef TrainingParameters_h
#define TrainingParameters_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "BooleanValue.h"
#include "ElementParser.h"

#include <cstdlib> // For EXIT_SUCCESS

class TrainingParameters : public ElementParser
{
public:
  typedef ElementParser SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const override
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== TrainingParameters ===" << std::endl;
    return indent + 2;
  }

  TrainingParameters(std::string method) : ElementParser(method.c_str() )
  {
    if( method == "ANNParameters" )
      {
      this->Add(new IntValue("Iterations", "20"), "Iterations");
      this->Add(new IntValue("MaximumVectorsPerEpoch", "2000"), "MaximumVectorsPerEpoch");
      this->Add(new IntValue("EpochIterations", 100), "EpochIterations");
      this->Add(new IntValue("ErrorInterval", 5), "ErrorInterval");
      this->Add(new FloatValue("DesiredError", 1.0), "DesiredError");
      this->Add(new FloatValue("ActivationSlope", 0.001), "ActivationSlope");
      this->Add(new FloatValue("ActivationMinMax", 1.0), "ActivationMinMax");
      this->Add(new IntValue("NumberOfHiddenNodes", 0), "NumberOfHiddenNodes");
      }
    else if( method == "RandomForestParameters" )
      {
      // TODO :: ADD proper initialization
      this->Add(new IntValue("MaxDepth", "5" ), "MaxDepth");
      this->Add(new IntValue("MinSampleCount", "5" ), "MinSampleCount");
      this->Add(new BooleanValue("UseSurrogates", false), "UseSurrogates");
      this->Add(new BooleanValue("CalcVarImportance", false), "CalcVarImportance");
      this->Add(new IntValue("MaxTreeCount", "5" ), "MaxTreeCount");
      }
    else
      {
      std::cout << "ERROR::: the training method `"
                << method
                << " does not supported "
                << std::endl;
      exit(EXIT_FAILURE);
      }
  }
};

#endif // TrainingParameters_h
