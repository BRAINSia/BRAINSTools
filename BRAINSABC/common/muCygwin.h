// Include this file before including stuff from mu

#ifndef __muCygwin_h
#define __muCygwin_h

#include "itkMacro.h"
#include <exception>
#include <iostream>

// Cygwin exception handling work-around
#undef itkExceptionMacro
#define itkExceptionMacro(x)                                    \
    {                                                             \
    std::cerr << "Exception: " x << std::endl;                  \
    std::cerr << "Possibly crashing about now..." << std::endl; \
    throw "exc";                                                \
    }

// TODO wrap main so that uncaught exception does not crash program
#define MU_DEFINE_MAIN                                          \
  int                                                           \
  main(int argc, char * *argv)                                   \
    {                                                             \
    int r = 0;                                                  \
    try                                                         \
      {                                                         \
      r = _mu_main(argc, argv);                                 \
      }                                                         \
    catch( itk::ExceptionObject & e )                          \
      {                                                         \
      std::cerr << e << std::endl;                              \
      return -1;                                                \
      }                                                         \
    catch( std::exception & e )                                \
      {                                                         \
      std::cerr << "Exception: " << e.what() << std::endl;      \
      return -1;                                                \
      }                                                         \
    catch( char *s )                                           \
      {                                                         \
      std::cerr << "Exception: " << s << std::endl;             \
      return -1;                                                \
      }                                                         \
    return r;                                                   \
    }

#endif
