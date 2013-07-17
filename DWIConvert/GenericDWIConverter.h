#ifndef __GenericDWIConverter_h
#define __GenericDWIConverter_h
#include "DWIConverter.h"

class GenericDWIConverter : public DWIConverter
{
public:
  GenericDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections) :  DWIConverter(allHeaders,
                                                                         inputFileNames,
                                                                         useBMatrixGradientDirections)
    {
    }
  virtual ~GenericDWIConverter() {}
protected:
  virtual void ExtractDWIData()
    {
      // throw; // don't call
    }
  virtual void AddFlagsToDictionary()
    {
      // throw; // don't call
    }
};

#endif // __GenericDWIConverter_h
