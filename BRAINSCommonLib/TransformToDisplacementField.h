#ifndef __TransformToDisplacementField.h
#define __TransformToDisplacementField .h
#include "itkIO.h"
#include "CrossOverAffineSystem.h"
#include <itkTransformToDeformationFieldSource.h>

/**
  * Go from any subclass of Transform, to the corresponding deformation field
  */
template <typename DisplacementFieldPointerType, typename TransformPointerType>
DisplacementFieldPointerType
TransformToDeformationField(itk::ImageBase<DisplacementFieldPointerType::ObjectType::ImageDimension> *templateImage,
                            TransformPointerType xfrm)
{
  typedef typename DisplacementFieldPointerType::ObjectType                   OutputType;
  typedef typename itk::TransformToDeformationFieldSource<OutputType, double> TodefType;
  typename TodefType::Pointer todef( TodefType::New() );
  todef->SetOutputParametersFromImage(templateImage);
  todef->SetTransform(xfrm);
  try
    {
    todef->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    throw err; // pass the buck up.
    }
  return todef->GetOutput();
}

#endif // TransformToDisplacementField.h
