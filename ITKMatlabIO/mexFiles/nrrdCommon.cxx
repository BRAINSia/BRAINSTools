#include "nrrdCommon.h"
#include "mex.h"

bool isGradientAxis(const unsigned int kindFlag)
{
  /* kindFlag == 1 indicates a domain -> Physical space
   * kindFlag == 2 indicates a space  -> Physical space
   * kindFlag == 4 indicates a list   -> Gradient Vector
   * kindFlag == 6 indicates a vector -> Gradient Vector */
  return kindFlag == 6 || kindFlag == 4;
}

mxClassID typeNtoM(const int ntype)
{
  mxClassID mtype;

  switch( ntype )
    {
    case nrrdTypeChar:
      mtype = mxINT8_CLASS;
      break;
    case nrrdTypeUChar:
      mtype = mxUINT8_CLASS;
      break;
    case nrrdTypeShort:
      mtype = mxINT16_CLASS;
      break;
    case nrrdTypeUShort:
      mtype = mxUINT16_CLASS;
      break;
    case nrrdTypeInt:
      mtype = mxINT32_CLASS;
      break;
    case nrrdTypeUInt:
      mtype = mxUINT32_CLASS;
      break;
    case nrrdTypeLLong:
      mtype = mxINT64_CLASS;
      break;
    case nrrdTypeULLong:
      mtype = mxUINT64_CLASS;
      break;
    case nrrdTypeFloat:
      mtype = mxSINGLE_CLASS;
      break;
    case nrrdTypeDouble:
      mtype = mxDOUBLE_CLASS;
      break;
    default:
      mtype = mxUNKNOWN_CLASS;
      break;
    }
  return mtype;
}
