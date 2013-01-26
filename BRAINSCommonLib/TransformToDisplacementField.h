#ifndef __TransformToDisplacementField_h
#define __TransformToDisplacementField_h

#include "itkIO.h"
#include "CrossOverAffineSystem.h"

#if (ITK_VERSION_MAJOR < 4)
#include <itkTransformToDeformationFieldSource.h>
#else
#include <itkTransformToDisplacementFieldSource.h>
#endif

/**
  * Go from any subclass of Transform, to the corresponding deformation field
  */
template <typename DisplacementFieldPointerType, typename TransformPointerType>
DisplacementFieldPointerType
TransformToDisplacementField(itk::ImageBase<DisplacementFieldPointerType::ObjectType::ImageDimension> *templateImage,
                             TransformPointerType xfrm)
{
  typedef typename DisplacementFieldPointerType::ObjectType OutputType;
#if (ITK_VERSION_MAJOR < 4)
  typedef typename itk::TransformToDeformationFieldSource<OutputType, double> TodefType;
#else
  typedef typename itk::TransformToDisplacementFieldSource<OutputType, double> TodefType;
#endif
  typename TodefType::Pointer todef( TodefType::New() );
  todef->SetOutputParametersFromImage(templateImage);
  todef->SetTransform(xfrm);
  try
    {
    todef->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    throw; // pass the buck up.
    }
  return todef->GetOutput();
}

#endif // TransformToDisplacementField_h
