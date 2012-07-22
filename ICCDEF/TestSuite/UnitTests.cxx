#include "vnl/vnl_sample.h"
#include "itkTestMain.h"
#include "itkConfigure.h"

void RegisterTests()
{
  REGISTER_TEST(itkVectorFFTWTest);
  REGISTER_TEST(itkIterativeInverseDisplacementFieldFilterTest);
  REGISTER_TEST(itkIterativeInverseDisplacementFieldFilterTest2);
}
