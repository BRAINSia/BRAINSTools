// this file defines the BRAINSDemonWarp-Tests for the test driver
// and all it expects is that you have a function called RegisterTests
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(VBRAINSDemonWarpTest);
}

#undef main
#define main VBRAINSDemonWarpTest
#include "../VBRAINSDemonWarp.cxx"
