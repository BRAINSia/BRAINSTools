#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(ApplyWarpTest);
}

#undef main
#define main ApplyWarpTest
#include "ApplyWarp.cxx"
