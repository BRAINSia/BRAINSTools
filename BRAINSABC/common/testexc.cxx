#include "muException.h"

int main()
{
  try
    {
    muExceptionMacro(<< "FOO");
    }
  catch( mu::Exception & e )
    {
    std::cerr << e << std::endl;
    return -1;
    }

  return 0;
}
