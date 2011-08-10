
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(gtractCoRegAnatomyTest);
}

#undef main
#define main gtractCoRegAnatomyTest
#include "../gtractCoRegAnatomy.cxx"
