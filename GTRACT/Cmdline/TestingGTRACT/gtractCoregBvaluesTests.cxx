
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(gtractCoregBvaluesTest);
}

#undef main
#define main gtractCoregBvaluesTest
#include "../gtractCoregBvalues.cxx"
