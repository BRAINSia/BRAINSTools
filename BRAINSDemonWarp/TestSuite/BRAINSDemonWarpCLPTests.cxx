// this file defines the BRAINSDemonWarp-Tests for the test driver
// and all it expects is that you have a function called RegisterTests
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(BRAINSDemonWarpTest);
}

#undef main
#define main BRAINSDemonWarpTest
#include "../BRAINSDemonWarp.cxx"
