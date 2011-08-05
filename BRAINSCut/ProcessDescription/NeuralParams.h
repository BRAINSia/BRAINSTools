#ifndef NeuralParams_h
#define NeuralParams_h
#include "StringValue.h"
#include "IntValue.h"
#include "FloatValue.h"
#include "CompoundObjectBase.h"

class NeuralParams : public CompoundObjectBase
{
public:
  typedef CompoundObjectBase SuperClass;
  virtual int PrintSelf(std::ostream & os, int indent) const
  {
    indent += SuperClass::PrintSelf(os, indent);
    os << this->PrintSpaces(indent) << "=== NeuralParams ===" << std::endl;
    return indent + 2;
  }

  NeuralParams() : CompoundObjectBase("NeuralParams")
  {
    this->Add(new FloatValue("MaskSmoothingValue", 0.0), "MaskSmoothingValue");
    // this->Add(new IntValue("GaussianSize",1),"GaussianSize");
    this->Add(new IntValue("GradientProfileSize", 1), "GradientProfileSize");
    this->Add(new IntValue("IrisSize", 1), "IrisSize");
    this->Add(new StringValue("TrainingVectorFilename", ""), "TrainingVectorFilename");
    this->Add(new StringValue("TestVectorFilename", ""), "TestVectorFilename");
    this->Add(new StringValue("TrainingModelFilename", ""), "TrainingModelFilename");
  }
};

#endif // NeuralParams_h
