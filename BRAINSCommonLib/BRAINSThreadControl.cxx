
#include "BRAINSThreadControl.h"

namespace BRAINSUtils
{
StackPushITKDefaultNumberOfThreads::StackPushITKDefaultNumberOfThreads(const int desiredCount) :
  m_originalThreadValue( itk::MultiThreader::GetGlobalDefaultNumberOfThreads() )
{
  int threadCount(-1);

  if( desiredCount > 0 )  // NOTE: Default is -1, which then uses the ITK default.
    {
    threadCount = desiredCount;
    }
  else
    {                                          //
    threadCount = this->m_originalThreadValue; // This is the old default.
    }
    { // Process the NSLOTS environmental varialble set by the SGE batch processing system
    int         NSLOTSThreadCount(-1);
    std::string numThreads;
    if( itksys::SystemTools::GetEnv("NSLOTS", numThreads) )
      {
      std::istringstream s(numThreads, std::istringstream::in);
      s >> NSLOTSThreadCount;
      }
    if( NSLOTSThreadCount > threadCount )
      {
      threadCount = NSLOTSThreadCount;
      }
    }
  if( threadCount > 0 )
    {
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threadCount);
    }
}

StackPushITKDefaultNumberOfThreads::~StackPushITKDefaultNumberOfThreads()
{
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(this->m_originalThreadValue);
}
}
