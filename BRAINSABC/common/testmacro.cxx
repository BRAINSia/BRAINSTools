#include <iostream>
#include <sstream>
#include <string>

#define muDisplayMacro(varname)          \
  std::cout << #varname << " = " << varname << std::endl;

#define muStringMacro(strname, s)    \
  std::string strname;        \
    {            \
    std::ostringstream outss;      \
    outss << "" s << std::ends;      \
    strname = outss.str();      \
    }

int main()
{
  int         x = 10;
  std::string s = "abcfoo";

  muDisplayMacro(x);
  muDisplayMacro(s);

  return 0;
}
