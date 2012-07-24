#include "BRAINSTransformConvertCLP.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageFileReader.h"
#include "itkBSplineDeformableTransform.h"
#include "itkIO.h"
#include "itkImageRegionIterator.h"
#include "GenericTransformImage.h"
#include "itkTranslationTransform.h"

//
// transform ranking,
// meaning a lower ranked transform can be
// converted to a higher ranked transform
// VersorRigid3D  = 1
// ScaleVersor3D = 2
// ScaleSkewVersor3D 3
// Affine 4
// BSpline 5
// BSplineROI 5

inline
bool
IsSameClass(const GenericTransformType *result,
            const GenericTransformType *source)
{
  return strcmp(result->GetNameOfClass(), source->GetNameOfClass() ) == 0;
}

inline
bool
IsClass(const GenericTransformType *xfrm, const char *className)
{
  return strcmp(xfrm->GetNameOfClass(), className) == 0;
}

//
// Convert from any type derived from MatrixOffsetTransformType to
// AffineTransform.
bool
ExtractTransform(AffineTransformType::Pointer & result,
                 const GenericTransformType *source)
{
  result->SetIdentity();
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }

  typedef AffineTransformType::Superclass MatrixOffsetTransformType;
  const MatrixOffsetTransformType *matBasePtr =
    dynamic_cast<const MatrixOffsetTransformType *>(source);
  if( matBasePtr == 0 )
    {
    return false;
    }

  result->SetCenter(matBasePtr->GetCenter() );
  result->SetMatrix(matBasePtr->GetMatrix() );
  result->SetTranslation(matBasePtr->GetTranslation() );
  return true;
}

//
// versor rigid 3d case.
bool
ExtractTransform(VersorRigid3DTransformType::Pointer & result,
                 const GenericTransformType *source)
{
  result->SetIdentity();
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }
  // this looks like it should be a convertible transform but
  // I'm not sure.
  typedef itk::TranslationTransform<double, 3> TransTransformType;
  if( IsClass(source, "TranslationTransform") )
    {
    const TransTransformType *translationXfrm =
      dynamic_cast<const TransTransformType *>(source);
    TransTransformType::OutputVectorType offset = translationXfrm->GetOffset();
    result->SetOffset(offset);
    return true;
    }
  // versor == rotation only
  if( IsClass(source, "VersorTransform") )
    {
    typedef itk::VersorTransform<double> VersorTransformType;
    const VersorTransformType *versorXfrm = dynamic_cast<const VersorTransformType *>(source);
    result->SetRotation(versorXfrm->GetVersor() );
    result->SetCenter(versorXfrm->GetCenter() );
    return true;
    }
  return false;
}

//
// scale versor case
bool
ExtractTransform(ScaleVersor3DTransformType::Pointer & result,
                 const GenericTransformType *source)
{
  result->SetIdentity();
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }
  if( IsClass(source, "VersorRigid3DTransform") )
    {
    const VersorRigid3DTransformType *versorRigidXfrm =
      dynamic_cast<const VersorRigid3DTransformType *>(source);
    result->SetRotation(versorRigidXfrm->GetVersor() );
    result->SetTranslation(versorRigidXfrm->GetTranslation() );
    result->SetCenter(versorRigidXfrm->GetCenter() );
    return true;
    }
  // otherwise try VersorRigidTransform
  VersorRigid3DTransformType::Pointer vrx = VersorRigid3DTransformType::New();
  if( ExtractTransform(vrx, source) ) // of VersorRigid3D conversion
                                      // works
    {
    // recurse to do this conversion
    return ExtractTransform(result, vrx.GetPointer() );
    }
  return false;
}

//
// scale skew versor case
bool
ExtractTransform(ScaleSkewVersor3DTransformType::Pointer & result,
                 const GenericTransformType *source)
{
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }
  // is it the parent?
  if( IsClass(source, "ScaleVersor3DTransform") )
    {
    const ScaleVersor3DTransformType *scaleVersorXfrm =
      dynamic_cast<const ScaleVersor3DTransformType *>(source);
    result->SetRotation(scaleVersorXfrm->GetVersor() );
    result->SetTranslation(scaleVersorXfrm->GetTranslation() );
    result->SetCenter(scaleVersorXfrm->GetCenter() );
    result->SetScale(scaleVersorXfrm->GetScale() );
    return true;
    }
  // otherwise try ScaleVersor conversion
  ScaleVersor3DTransformType::Pointer svx = ScaleVersor3DTransformType::New();
  if( ExtractTransform(svx, source) ) // of VersorRigid3D conversion
                                      // works
    {
    // recurse to do this conversion
    return ExtractTransform(result, svx.GetPointer() );
    }
  return false;
}

//
// scale skew versor case
bool
ExtractTransform(BSplineTransformType::Pointer & result,
                 const GenericTransformType *source)
{
  if( IsSameClass(result.GetPointer(), source) )
    {
    const BSplineTransformType *sourceBSpline =
      dynamic_cast<const BSplineTransformType *>(source);
    result->SetFixedParameters(sourceBSpline->GetFixedParameters() );
    result->SetIdentity();
    result->SetParametersByValue(sourceBSpline->GetParameters() );
    result->SetBulkTransform(sourceBSpline->GetBulkTransform() );
    return true;
    }
  else
    {
    result->SetIdentity();
    }
  typedef BSplineTransformType::BulkTransformType BulkXfrmType;

  const BulkXfrmType *bulkXfrm =
    dynamic_cast<const BulkXfrmType *>(source);
  if( bulkXfrm == 0 )
    {
    return false;
    }
  result->SetBulkTransform(bulkXfrm);
  return true;
}

void
TransformConvertError(GenericTransformType *inputXfrm,
                      const std::string & targetClassName)
{
  std::cerr << "Can't convert transform of type "
            << inputXfrm->GetTransformTypeAsString()
            << " to "
            << targetClassName
            << std::endl;
}

#define CHECK_PARAMETER_IS_SET(parameter, message) \
  if( parameter == "" )                             \
    {                                             \
    std::cerr << message << std::endl;            \
    return EXIT_FAILURE;                          \
    }

int main(int argc, char *argv[])
{
  PARSE_ARGS;

  CHECK_PARAMETER_IS_SET(inputTransform,
                         "Missing inputTransform parameter");
  CHECK_PARAMETER_IS_SET(outputTransformType,
                         "Missing outpuTransformType");

  // read the input transform
  itk::TransformFileReader::Pointer reader =
    itk::TransformFileReader::New();
  reader->SetFileName(inputTransform.c_str() );
  reader->Update();
  itk::TransformFileReader::TransformListType *transformList =
    reader->GetTransformList();
  GenericTransformType::Pointer inputXfrm = dynamic_cast<GenericTransformType *>(transformList->front().GetPointer() );

  // Handle BSpline type
  BSplineTransformType::Pointer bsplineInputXfrm =
    dynamic_cast<BSplineTransformType *>(inputXfrm.GetPointer() );
  if( bsplineInputXfrm.IsNotNull() )
    {
    transformList->pop_front();
    if( transformList->size() == 0 )
      {
      std::cerr << "Error, the second transform needed for BSplineDeformableTransform is missing." << std::endl;
      return EXIT_FAILURE;
      }
    BSplineTransformType::BulkTransformType::Pointer bulkXfrm =
      dynamic_cast<BSplineTransformType::BulkTransformType *>(transformList->front().GetPointer() );
    if( bulkXfrm.IsNull() )
      {
      std::cerr << "Error, the second transform is not a bulk transform" << std::endl;
      }
    bsplineInputXfrm->SetBulkTransform(bulkXfrm);
    inputXfrm = bsplineInputXfrm.GetPointer();
    }

  if( outputTransformType == "DisplacementField" )
    {
    CHECK_PARAMETER_IS_SET(referenceVolume,
                           "Missing referenceVolume needed for Displacement Field output");
    CHECK_PARAMETER_IS_SET(displacementVolume,
                           "Missing displacementVolume needed for Displacement Field output");

    typedef itk::Image<short, 3> ReferenceImageType;
    ReferenceImageType::Pointer referenceImage = itkUtil::ReadImage<ReferenceImageType>(referenceVolume);
    if( referenceImage.IsNull() )
      {
      std::cerr << "Can't read Reference Volume " << referenceVolume << std::endl;
      return EXIT_FAILURE;
      }
    // Allocate Displacement Field
    typedef itk::Vector<float, 3>     VectorType;
    typedef itk::Image<VectorType, 3> DisplacementFieldType;
    DisplacementFieldType::Pointer displacementField =
      itkUtil::AllocateImageFromExample<ReferenceImageType, DisplacementFieldType>(referenceImage);

    typedef itk::ImageRegionIterator<DisplacementFieldType> DisplacementIteratorType;
    for( DisplacementIteratorType it(displacementField, displacementField->GetLargestPossibleRegion() );
         !it.IsAtEnd(); ++it )
      {
      DisplacementFieldType::IndexType dispIndex = it.GetIndex();
      DisplacementFieldType::PointType fixedPoint, movingPoint;
      displacementField->TransformIndexToPhysicalPoint(dispIndex, fixedPoint);
      movingPoint = inputXfrm->TransformPoint(fixedPoint);
      VectorType displacement = movingPoint - fixedPoint;
      it.Set(displacement);
      }

    try
      {
      itkUtil::WriteImage<DisplacementFieldType>(displacementField, displacementVolume);
      }
    catch( ... )
      {
      std::cerr << "Error writing displacement field " << displacementVolume << std::endl;
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }

  CHECK_PARAMETER_IS_SET(outputTransform,
                         "Missing outputTransform parameter");

  GenericTransformType::Pointer outputXfrm;

  if( outputTransformType == "Affine" )
    {
    AffineTransformType::Pointer affineXfrm = AffineTransformType::New();
    if( ExtractTransform(affineXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError(inputXfrm, "Affine Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = affineXfrm.GetPointer();
    }
  else if( outputTransformType == "VersorRigid" )
    {
    VersorRigid3DTransformType::Pointer versorRigidXfrm =
      VersorRigid3DTransformType::New();
    if( ExtractTransform(versorRigidXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError(inputXfrm, "VersorRigid3D Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = versorRigidXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleVersor" )
    {
    ScaleVersor3DTransformType::Pointer scaleVersorXfrm =
      ScaleVersor3DTransformType::New();
    if( ExtractTransform(scaleVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError(inputXfrm, "ScaleVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleVersorXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleSkewVersor" )
    {
    ScaleSkewVersor3DTransformType::Pointer scaleSkewVersorXfrm =
      ScaleSkewVersor3DTransformType::New();
    if( ExtractTransform(scaleSkewVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError(inputXfrm, "ScaleSkewVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleSkewVersorXfrm.GetPointer();
    }
#if 0
  else if( outputTransformType == "BSplineDeformable" )
    {
    BSplineTransformType::Pointer bsplineXfrm =
      BSplineTransformType::New();
    if( ExtractTransform(bsplineXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError(inputXfrm, "BSplineDeformable Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = bsplineXfrm.GetPointer();
    }
#endif
  if( outputTransformType == "Same" )
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName(outputTransform);
    for( itk::TransformFileReader::TransformListType::iterator it = transformList->begin();
         it != transformList->end(); ++it )
      {
      transformWriter->AddTransform( (*it).GetPointer() );
      }

    try
      {
      transformWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Can't write " << outputTransform << excp.GetDescription() << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    // write the resulting transform.
    itk::WriteTransformToDisk(outputXfrm.GetPointer(), outputTransform);
    }
  return EXIT_SUCCESS;
}
