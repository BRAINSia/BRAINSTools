#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(BRAINSResampleTest);
}

#undef main
#define main BRAINSResampleTest
#include "../BRAINSResample.cxx"
