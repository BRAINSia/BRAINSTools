#ifndef _DEBUG_h
#define _DEBUG_h

// Hints from
//
// http://www.hausmilbe.net/Blog/text_20090224_Debugging_Matlab_Mex_files_with_Valgrind_and_DDD.php
/* macros for debugging */
#define xxDEBUG
#ifdef xxDEBUG

#include "mex.h"
#define MPRINT_MSG(_msg)  mexPrintf("[%s|%s, %3d]: %s\n", __FILE__, __func__, __LINE__, _msg);
#define MPRINT_INT(_val)  mexPrintf("[%s|%s, %3d]: %s = %d\n", __FILE__, __func__, __LINE__, #_val, _val);
#define MPRINT_DBL(_val)  mexPrintf("[%s|%s, %3d]: %s = %f\n", __FILE__, __func__, __LINE__, #_val, _val);
#define MPRINT_PTR(_val)  mexPrintf("[%s|%s, %3d]: %s = %p\n", __FILE__, __func__, __LINE__, #_val, _val);

#else       /* #ifdef DEBUG */

/* do not print messages and values */
#define MPRINT_MSG(_msg)
#define MPRINT_INT(_val)
#define MPRINT_DBL(_val)
#define MPRINT_PTR(_val)

/* assert needs that when no debugging is desired and therefore, assertion is
  ignored */
#ifndef NDEBUG
#define NDEBUG 1
#endif

#endif       /* #ifdef DEBUG */

#include "itkMetaDataDictionary.h"
#include <cstdio>
#include <iostream>
#include <sstream>
enum newlines { none, nl };

//#define __DEBUG 1

#if defined(__DEBUG)

inline unsigned char
space_or_newline(const unsigned index, const unsigned limit)
{
  return (index < limit - 1) ? ' ' : '\n';
}

template <typename T>
inline
void DEBUG_PRINT(const T & var, newlines n = none)
{
  std::stringstream ss;

  ss << var;
  if( n == nl )
    {
    ss << std::endl;
    }
  fputs(ss.str().c_str(), stdout);
  mexWarnMsgTxt(ss.str().c_str());

  fflush(stdout);
}

template <typename T>
inline
void DEBUG_PRINT(const T & var, unsigned index, unsigned limit)
{
  std::stringstream ss;

  ss << var << space_or_newline(index, limit);
  fputs(ss.str().c_str(), stdout);
  mexWarnMsgTxt(ss.str().c_str());
  fflush(stdout);
}

template <typename TImage>
inline
void DEBUG_PRINTDICT(const TImage *im)
{
  std::stringstream               ss;
  const itk::MetaDataDictionary & thisDic = im->GetMetaDataDictionary();

  for( itk::MetaDataDictionary::ConstIterator it = thisDic.Begin();
       it != thisDic.End(); ++it )
    {
    ss << it->first << " " << typeid(it->second).name() << " ";
    it->second->Print(ss);
    }
  fputs(ss.str().c_str(), stdout);
  mexWarnMsgTxt(ss.str().c_str());
  fflush(stdout);
}

#else
template <typename T>
inline
void DEBUG_PRINT(const T &, newlines)
{
}

template <typename T>
inline
void DEBUG_PRINT(const T & )
{
}

template <typename T>
inline
void DEBUG_PRINT(const T &, unsigned, unsigned)
{
}

template <typename TImage>
inline
void DEBUG_PRINTDICT(const TImage *)
{
}

#endif
#endif // _DEBUG_h
