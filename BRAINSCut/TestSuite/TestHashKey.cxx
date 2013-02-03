#include "FeatureInputVector.h"

int
main(int /*argc*/, char * [] /*argv*/)
{
  FeatureInputVector testFeatureInputVector;

  bool allTestPass = true;

  if( !testFeatureInputVector.DoScalingUnitTests() )
    {
    allTestPass = false;
    }

  return testFeatureInputVector.DoUnitTests();
}
