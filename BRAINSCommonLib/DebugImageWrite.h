#ifndef __DebugImageWrite_h
#define __DebugImageWrite_h
#include "BRAINSCommonLib.h"
#if defined(BRAINS_DEBUG_IMAGE_WRITE)
#include "itkIO.h"
#include "itksys/SystemTools.hxx"

#define DEFINE_DEBUG_IMAGE_COUNTER \
  namespace DebugImageWrite        \
  {                                \
  int fileSequenceNumber = 0;    \
  }

namespace DebugImageWrite
{
extern int fileSequenceNumber;
inline std::string twodigits(unsigned int x)
{
  std::string rval;

  rval = '0' + static_cast<char>(x / 10);
  rval += '0' + static_cast<char>(x % 10);
  return rval;
}

extern int fileSequenceNumber;

template <typename ImageType>
void DebugOutput(int LINE,
                 const char *FILE,
                 typename ImageType::Pointer img,
                 int imageIndex = -1,
                 const char *name = 0)
{
  std::string fname(FILE);

  fname = itksys::SystemTools::GetFilenameName(fname);
  std::stringstream filename;
  filename << "DBG_";
  filename << twodigits(fileSequenceNumber);
  filename << "_";
  filename << name;
  if( imageIndex != -1 )
    {
    filename << twodigits(fileSequenceNumber);
    }
  filename << "_" << fname
           << "_" << LINE << ".nii.gz";
  std::cerr << "Writing " << filename.str() << " "
            << img.GetPointer() << std::endl;
  itkUtil::WriteImage<ImageType>(img, filename.str() );
  ++fileSequenceNumber;
}
}

#define DebugOutput(imageType, image) \
  DebugImageWrite::DebugOutput<imageType>(__LINE__, __FILE__, image, -1,#image)
#define DebugOutputN(imageType, image, N, name) \
  DebugImageWrite::DebugOutput<imageType>(__LINE__, __FILE__, image, N, #name)
#define DebugOutputWName(imageType, image, name) \
  DebugImageWrite::DebugOutput<imageType>(__LINE__, __FILE__, image, -1,#name)

#else
#define DEFINE_DEBUG_IMAGE_COUNTER
#define DebugOutput(imageType, image)
#define DebugOutputN(imageType, image, N, name)
#define DebugOutputWName(imageType, image, name)
#endif // BRAINS_DEBUG_IMAGE_WRITE
#endif // __DebugImageWrite_h
