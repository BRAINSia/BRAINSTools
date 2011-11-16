#ifndef BRAINSCutExceptionHandler_h
#define BRAINSCutExceptionHandler_h
#include <string>
#include <ostream>

class BRAINSCutExceptionStringHandler
{
public:
  BRAINSCutExceptionStringHandler(const std::string & errorString);
  BRAINSCutExceptionStringHandler(const char *errorString);
  const std::string & Error() const;

  friend std::ostream & operator<<(std::ostream& stream, BRAINSCutExceptionStringHandler ob);

private:
  std::string m_ErrorString;
};

#endif
