// this file defines the ImageCompare-Tests for the test driver
// and all it expects is that you have a function called RegisterTests
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include <iostream>
#include "itkTestMain.h"

// A dummy method that does nothing. We intend doing a comparison which is done
// by the itkTestMain.h. (It requires a main program, prior to the compare)
int Dummy(int, char *[])
{
  return EXIT_SUCCESS;
}

void RegisterTests()
{
  REGISTER_TEST(Dummy);
}
