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

template<class TScalarType>
inline
bool
IsSameClass(const itk::Transform< TScalarType, 3, 3 > *result,
            const itk::Transform< TScalarType, 3, 3 > *source)
{
  return strcmp(result->GetNameOfClass(), source->GetNameOfClass() ) == 0;
}

template<class TScalarType>
inline
bool
IsClass(const itk::Transform< TScalarType, 3, 3 > *xfrm, const char *className)
{
  return strcmp(xfrm->GetNameOfClass(), className) == 0;
}

template<class TScalarType>
void
TransformConvertError(const itk::Transform< TScalarType, 3, 3 > *inputXfrm,
                      const std::string & targetClassName)
{
  std::cerr << "Can't convert transform of type "
  << inputXfrm->GetTransformTypeAsString()
  << " to "
  << targetClassName
  << std::endl;
}

//
// Convert from any type derived from MatrixOffsetTransformType to
// AffineTransform.
template<class TScalarType>
bool
ExtractTransform(typename itk::AffineTransform< TScalarType, 3 >::Pointer &result,
                 const itk::Transform< TScalarType, 3, 3 > *source)
{
  result->SetIdentity();
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }

  typedef itk::AffineTransform< TScalarType, 3 > AffineTransformType;
  typedef typename AffineTransformType::Superclass MatrixOffsetTransformType;
  const MatrixOffsetTransformType *matBasePtr = dynamic_cast<const MatrixOffsetTransformType *>(source);
  if( matBasePtr == 0 )
    {
    return false;
    }

  result->SetCenter( matBasePtr->GetCenter() );
  result->SetMatrix( matBasePtr->GetMatrix() );
  result->SetTranslation( matBasePtr->GetTranslation() );
  return true;
}

//
// versor rigid 3d case.
template<class TScalarType>
bool
ExtractTransform(typename itk::VersorRigid3DTransform<TScalarType>::Pointer & result,
                 const itk::Transform< TScalarType, 3, 3 > *source)
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
  typedef itk::TranslationTransform<TScalarType, 3> TransTransformType;
  if( IsClass(source, "TranslationTransform") )
    {
    const TransTransformType *translationXfrm = dynamic_cast<const TransTransformType *>(source);
    typename TransTransformType::OutputVectorType offset = translationXfrm->GetOffset();
    result->SetOffset( offset );
    return true;
    }
  // versor == rotation only
  if( IsClass(source, "VersorTransform") )
    {
    typedef itk::VersorTransform<TScalarType> VersorTransformType;
    const VersorTransformType *versorXfrm = dynamic_cast<const VersorTransformType *>(source);

    result->SetRotation( versorXfrm->GetVersor() );
    result->SetCenter( versorXfrm->GetCenter() );
    return true;
    }
  return false;
}

//
// scale versor case
template<class TScalarType>
bool
ExtractTransform(typename itk::ScaleVersor3DTransform<TScalarType>::Pointer & result,
                 const itk::Transform< TScalarType, 3, 3 > *source)
{
  result->SetIdentity();
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }

  typedef itk::VersorRigid3DTransform<TScalarType> VersorRigid3DTransformType;
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
  typename VersorRigid3DTransformType::Pointer vrx = VersorRigid3DTransformType::New();
  if( ExtractTransform<TScalarType>(vrx, source) ) // of VersorRigid3D conversion
                                      // works
    {
    // recurse to do this conversion
    return ExtractTransform<TScalarType>(result, vrx.GetPointer() );
    }
  return false;
}

//
// scale skew versor case
template<class TScalarType>
bool
ExtractTransform(typename itk::ScaleSkewVersor3DTransform< TScalarType >::Pointer & result,
                 const itk::Transform< TScalarType, 3, 3 > *source)
{
  // always able to convert to same type
  if( IsSameClass(result.GetPointer(), source) )
    {
    result->SetParameters( source->GetParameters() );
    result->SetFixedParameters( source->GetFixedParameters() );
    return true;
    }

  // is it the parent?
  typedef itk::ScaleVersor3DTransform<TScalarType> ScaleVersor3DTransformType;
  if( IsClass(source, "ScaleVersor3DTransform") )
    {
    const ScaleVersor3DTransformType *scaleVersorXfrm = dynamic_cast<const ScaleVersor3DTransformType *>(source);
    result->SetRotation(scaleVersorXfrm->GetVersor() );
    result->SetTranslation(scaleVersorXfrm->GetTranslation() );
    result->SetCenter(scaleVersorXfrm->GetCenter() );
    result->SetScale(scaleVersorXfrm->GetScale() );
    return true;
    }
  // otherwise try ScaleVersor conversion
  typename ScaleVersor3DTransformType::Pointer svx = ScaleVersor3DTransformType::New();
  if( ExtractTransform<TScalarType>(svx, source) ) // of VersorRigid3D conversion
                                                                               // works
    {
    // recurse to do this conversion
    return ExtractTransform<TScalarType>(result, svx.GetPointer() );
    }
  return false;
}

#define CHECK_PARAMETER_IS_SET(parameter, message) \
  if( parameter == "" )                             \
    {                                             \
    std::cerr << message << std::endl;            \
    return EXIT_FAILURE;                          \
    }

template<class TScalarType>
int
DoConversion( int argc, char *argv[] )
{
  PARSE_ARGS;

  typedef itk::Transform< TScalarType, 3, 3 >                                 GenericTransformType;
  typedef itk::BSplineDeformableTransform< TScalarType,
                                           GenericTransformImageNS::SpaceDimension,
                                           GenericTransformImageNS::SplineOrder>    BSplineTransformType;

  typedef itk::AffineTransform< TScalarType, 3 >            AffineTransformType;
  typedef itk::VersorRigid3DTransform< TScalarType >        VersorRigid3DTransformType;
  typedef itk::ScaleVersor3DTransform< TScalarType >        ScaleVersor3DTransformType;
  typedef itk::ScaleSkewVersor3DTransform< TScalarType >    ScaleSkewVersor3DTransformType;

  // read the input transform
  typedef itk::TransformFileReaderTemplate<TScalarType>  TransformFileReaderType;
  typename TransformFileReaderType::Pointer reader = TransformFileReaderType::New();
  reader->SetFileName(inputTransform.c_str() );
  reader->Update();
  typename TransformFileReaderType::TransformListType *transformList = reader->GetTransformList();
  typename GenericTransformType::Pointer inputXfrm = dynamic_cast<GenericTransformType *>( transformList->front().GetPointer() );

  std::cout << "------------------------ " << std::endl;
  std::cout << "Input Transform Type Saved on Memory ==> " << inputXfrm->GetTransformTypeAsString() << std::endl;
  std::cout << "* Input transform parameters: " << inputXfrm->GetParameters() << std::endl;
  std::cout << "* Input transform fixed parameters: " << inputXfrm->GetFixedParameters() << std::endl;
  std::cout << "------------------------ " << std::endl;

  // Handle BSpline type
  typename BSplineTransformType::Pointer bsplineInputXfrm = dynamic_cast<BSplineTransformType *>( inputXfrm.GetPointer() );
  if( bsplineInputXfrm.IsNotNull() )
    {
    transformList->pop_front();
    if( transformList->size() == 0 )
      {
      std::cerr << "Error, the second transform needed for BSplineDeformableTransform is missing." << std::endl;
      return EXIT_FAILURE;
      }
    typename BSplineTransformType::BulkTransformType::Pointer bulkXfrm =
      dynamic_cast<typename BSplineTransformType::BulkTransformType *>(transformList->front().GetPointer() );
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

  //output transform processing
  typename GenericTransformType::Pointer outputXfrm;

  if( outputTransformType == "Affine" )
    {
    typename AffineTransformType::Pointer affineXfrm = AffineTransformType::New();
    if( ExtractTransform<TScalarType>(affineXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "Affine Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = affineXfrm.GetPointer();
    }
  else if( outputTransformType == "VersorRigid" )
    {
    typename VersorRigid3DTransformType::Pointer versorRigidXfrm = VersorRigid3DTransformType::New();
    if( ExtractTransform<TScalarType>(versorRigidXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "VersorRigid3D Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = versorRigidXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleVersor" )
    {
    typename ScaleVersor3DTransformType::Pointer scaleVersorXfrm = ScaleVersor3DTransformType::New();
    if( ExtractTransform<TScalarType>( scaleVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "ScaleVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleVersorXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleSkewVersor" )
    {
    typename ScaleSkewVersor3DTransformType::Pointer scaleSkewVersorXfrm = ScaleSkewVersor3DTransformType::New();
    if( ExtractTransform<TScalarType>( scaleSkewVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "ScaleSkewVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleSkewVersorXfrm.GetPointer();
    }

  if( outputTransformType == "Same" )
    {
    typedef typename itk::TransformFileWriterTemplate<TScalarType> TransformWriterType;
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName(outputTransform);
    for( typename itk::TransformFileReaderTemplate<TScalarType>::TransformListType::iterator it = transformList->begin();
         it != transformList->end(); ++it )
      {
      typename GenericTransformType::Pointer outXfrm = dynamic_cast<GenericTransformType *>( (*it).GetPointer() );
      transformWriter->AddTransform( outXfrm );
      //
      std::cout << "Output Transform Type Written to the Disk ==> " << outXfrm->GetTransformTypeAsString() << std::endl;
      std::cout << "* Output transform parameters: " << outXfrm->GetParameters() << std::endl;
      std::cout << "* Output transform fixed parameters: " << outXfrm->GetFixedParameters() << std::endl;
      std::cout << "------------------------ " << std::endl;
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
    std::cout << "Output Transform Type Written to the Disk ==> " << outputXfrm->GetTransformTypeAsString() << std::endl;
    std::cout << "* Output transform parameters: " << outputXfrm->GetParameters() << std::endl;
    std::cout << "* Output transform fixed parameters: " << outputXfrm->GetFixedParameters() << std::endl;
    std::cout << "------------------------ " << std::endl;
    //
    itk::WriteTransformToDisk<TScalarType>(outputXfrm.GetPointer(), outputTransform);
    }
  return EXIT_SUCCESS;
}


int main(int argc, char *argv[])
{
  PARSE_ARGS;

  CHECK_PARAMETER_IS_SET(inputTransform,
                         "Missing inputTransform parameter");
  CHECK_PARAMETER_IS_SET(outputTransform,
                         "Missing outputTransform parameter");
  CHECK_PARAMETER_IS_SET(outputTransformType,
                         "Missing outpuTransformType");
  CHECK_PARAMETER_IS_SET(outputPrecisionType,
                         "Missing outputPrecisionType");

  if( outputPrecisionType == "double" )
    {
    return DoConversion<double>( argc, argv );
    }
  else if( outputPrecisionType == "float" )
    {
    return DoConversion<float>( argc, argv );
    }
  else
    {
    std::cerr << "Error: Invalid parameter for output precision type." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
