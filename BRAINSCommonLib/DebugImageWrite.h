#ifndef __DebugImageWrite_h
#define __DebugImageWrite_h
#include "BRAINSCommonLib.h"
#if defined(BRAINS_DEBUG_IMAGE_WRITE)
#include "itkIO.h"
#include "itksys/SystemTools.hxx"

#define DEFINE_DEBUG_IMAGE_COUNTER int fileSequenceNumber = 0
extern int fileSequenceNumber;
inline std::string twodigits(unsigned int x)
{
  std::string rval;

  rval = '0' + static_cast<char>(x / 10);
  rval += '0' + static_cast<char>(x % 10);
  return rval;
}

extern int fileSequenceNumber;

#define DebugOutput(imageType, image)                           \
    {                                                            \
    typename imageType::Pointer im = image;                    \
    std::string fname(__FILE__);                               \
    fname = itksys::SystemTools::GetFilenameName(fname);       \
    std::stringstream filename;                                \
    filename << "DBG_";                                        \
    filename << twodigits(fileSequenceNumber);                 \
    filename << "_"                                            \
             << #image << "_" << fname                         \
             << "_" << __LINE__ << ".nii.gz";                  \
    std::cerr << "Writing " << filename.str();                 \
    itkUtil::WriteImage<imageType>(im, filename.str() );         \
    ++fileSequenceNumber;                                      \
    }

#define DebugOutputN(imageType, image, N, name)                      \
    {                                                               \
    typename imageType::Pointer im = image;                       \
    std::string fname(__FILE__);                                  \
    fname = itksys::SystemTools::GetFilenameName(fname);          \
    std::stringstream filename;                                   \
    filename << "DBG_" << twodigits(fileSequenceNumber) << "_"    \
             << #name << "_" << N << "_" << fname                 \
             << "_" << __LINE__ << ".nii.gz";                     \
    std::cerr << "Writing " << filename.str();                    \
    itkUtil::WriteImage<imageType>(im, filename.str() );            \
    ++fileSequenceNumber;                                         \
    }

#define DebugOutputWName(imageType, image, name)                    \
    {                                                               \
    typename imageType::Pointer im = image;                       \
    std::string fname(__FILE__);                                  \
    fname = itksys::SystemTools::GetFilenameName(fname);          \
    std::stringstream filename;                                   \
    filename << "DBG_" << twodigits(fileSequenceNumber) << "_"    \
             << #name "_" << fname                                \
             << "_" << __LINE__ << ".nii.gz";                     \
    std::cerr << "Writing " << filename.str();                    \
    itkUtil::WriteImage<imageType>(im, filename.str() );            \
    ++fileSequenceNumber;                                         \
    }
#else
#define DEFINE_DEBUG_IMAGE_COUNTER
#define DebugOutput(imageType, image)
#define DebugOutputN(imageType, image, N, name)
#define DebugOutputWName(imageType, image, name)
#endif // BRAINS_DEBUG_IMAGE_WRITE
#endif // __DebugImageWrite_h
