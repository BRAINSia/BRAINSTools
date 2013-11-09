#include <BRAINSCommonLib.h>


// If ITK supplied MGHIO, then register it for inclusion in BRAINSTools for use.
#if defined(ITK_HAS_MGHIO)
#include "itkMGHImageIOFactory.h"
void BRAINSRegisterAlternateIO(void)
{
  itk::MGHImageIOFactory::RegisterOneFactory();
}
#else
void BRAINSRegisterAlternateIO(void)
{
}
#endif
