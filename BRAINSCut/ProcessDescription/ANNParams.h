#ifndef ANNParams_h
#define ANNParams_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "CompoundObjectBase.h"

class ANNParams : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== ANNParams ===" << std::endl;
    return indent + 2;
  }

  ANNParams() : CompoundObjectBase("ANNParams")
  {
    // this->Add(new IntValue("VectorSize", "0"), "VectorSize");
    this->Add(new IntValue("Iterations", "20"), "Iterations");
    this->Add(new IntValue("MaximumVectorsPerEpoch", "2000"), "MaximumVectorsPerEpoch");
    this->Add(new IntValue("EpochIterations", 100), "EpochIterations");
    this->Add(new IntValue("ErrorInterval", 5), "ErrorInterval");
    this->Add(new FloatValue("DesiredError", 1.0), "DesiredError");
    this->Add(new FloatValue("ActivationSlope", 0.001), "ActivationSlope");
    this->Add(new FloatValue("ActivationMinMax", 1.0), "ActivationMinMax");
    this->Add(new IntValue("NumberOfHiddenNodes", 0), "NumberOfHiddenNodes");
    // this->Add(new FloatValue("LearningRate",0.3),"LearningRate");
    // this->Add(new FloatValue("MomentumRate",0.15),"MomentumRate");
  }
};

#endif // ANNParams_h
