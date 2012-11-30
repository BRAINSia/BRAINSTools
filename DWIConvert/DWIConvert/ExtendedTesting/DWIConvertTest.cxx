
#include <iostream>
#include "itkTestMain.h"
// #include "DWIConvert.cxx"

/*
extern "C" MODULE_IMPORT int ModuleEntryPoint(int, char* []);

void RegisterTests()
{
  StringToTestFunctionMap["ModuleEntryPoint"] = ModuleEntryPoint;
}
*/

void RegisterTests()
{
  REGISTER_TEST(DWIConvertTest);
}

#undef main
#define main DWIConvertTest
#include "../DWIConvert.cxx"
