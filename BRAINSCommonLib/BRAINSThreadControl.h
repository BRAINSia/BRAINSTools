#ifndef BRAINSThreadControl_h
#define BRAINSThreadControl_h

#include <itksys/SystemTools.hxx>
#include <sstream>
#include "itkMultiThreader.h"

namespace BRAINSUtils
{
inline
void SetThreadCount(int desiredCount)
{
  int threadCount(-1);

  if( desiredCount > 0 )  // NOTE: Default is -1, which then uses the
  // ITK default.
    {
    threadCount = desiredCount;
    }
  else
    {
    std::string numThreads;
    if( itksys::SystemTools::GetEnv("NSLOTS", numThreads) )
      {
      std::istringstream s(numThreads, std::istringstream::in);
      s >> threadCount;
      }
    }
  if( threadCount > 0 )
    {
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threadCount);
    }
}
}

#endif // BRAINSThreadControl_h
