
#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST(gtractCoregBvaluesTest);
}

#undef main
#define main gtractCoregBvaluesTest
#include "../gtractCoregBvalues.cxx"
