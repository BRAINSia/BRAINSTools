#include "WriteIdentityTransformCLP.h"
#include "itkTransformFileWriter.h"
#include "itkAffineTransform.h"

int
main(int argc, char * *argv)
{
  PARSE_ARGS;
  typedef itk::AffineTransform<double, 3> TransformationType;

  TransformationType::Pointer transform = TransformationType::New();

  transform->SetIdentity();

  typedef itk::TransformFileWriter WriterType;
  WriterType::Pointer writer =  WriterType::New();

  writer->SetFileName("./identityTransformation.mat");
  writer->SetInput( transform);
  writer->Update();
}
