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

template<class TInputScalarType, class TOutputScalarType>
inline
bool
IsSameClass(const itk::Transform< TOutputScalarType, 3, 3 > *result,
            const itk::Transform< TInputScalarType, 3, 3 > *source)
{
  return strcmp(result->GetNameOfClass(), source->GetNameOfClass() ) == 0;
}

template<class TInputScalarType>
inline
bool
IsClass(const itk::Transform< TInputScalarType, 3, 3 > *xfrm, const char *className)
{
  return strcmp(xfrm->GetNameOfClass(), className) == 0;
}

template<class TInputScalarType>
void
TransformConvertError(const itk::Transform< TInputScalarType, 3, 3 > *inputXfrm,
                      const std::string & targetClassName)
{
  std::cerr << "Can't convert transform of type "
  << inputXfrm->GetTransformTypeAsString()
  << " to "
  << targetClassName
  << std::endl;
}

//
// Type Conversions //
//
// converts double precision type parameters to single precision type
template<class TInputScalarType, class TOutputScalarType>
itk::OptimizerParameters< TOutputScalarType >
ParametersTypeConvertor(const itk::OptimizerParameters< TInputScalarType >  &sourceParams)
{
  itk::OptimizerParameters< TOutputScalarType > outputParams;
  outputParams.SetSize( sourceParams.GetSize() );
  outputParams.Fill(0);
  for( itk::SizeValueType i = 0; i < sourceParams.GetSize(); ++i )
    {
    outputParams[i] = (TOutputScalarType)( sourceParams[i] );
    }
  return outputParams;
}

// two specializations when both in/output types are the same. No casting needed.
template<>
itk::OptimizerParameters<double>
ParametersTypeConvertor<double, double>(const itk::OptimizerParameters<double> &sourceParams)
{
  return sourceParams;
}

template<>
itk::OptimizerParameters<float>
ParametersTypeConvertor<float, float>(const itk::OptimizerParameters<float> &sourceParams)
{
  return sourceParams;
}

// converts double precision type points to single precision type
template<class TInputScalarType, class TOutputScalarType>
itk::Point<TOutputScalarType, 3>
PointTypeConvertor(const itk::Point<TInputScalarType, 3> &sourcePoint)
{
  itk::Point<TOutputScalarType, 3> outputPoint;
  outputPoint[0] = (TOutputScalarType)( sourcePoint[0] );
  outputPoint[1] = (TOutputScalarType)( sourcePoint[1] );
  outputPoint[2] = (TOutputScalarType)( sourcePoint[2] );
  return outputPoint;
}

template<>
itk::Point<double, 3>
PointTypeConvertor<double, double>(const itk::Point<double, 3> &sourcePoint)
{
  return sourcePoint;
}

template<>
itk::Point<float, 3>
PointTypeConvertor<float, float>(const itk::Point<float, 3> &sourcePoint)
{
  return sourcePoint;
}

// converts double precision type vectors to single precision type
template<class TInputScalarType, class TOutputScalarType>
itk::Vector<TOutputScalarType, 3>
VectorTypeConvertor(const itk::Vector<TInputScalarType, 3> &sourceVector)
{
  itk::Vector<TOutputScalarType, 3> OutputVector;
  OutputVector[0] = (TOutputScalarType)( sourceVector[0] );
  OutputVector[1] = (TOutputScalarType)( sourceVector[1] );
  OutputVector[2] = (TOutputScalarType)( sourceVector[2] );
  return OutputVector;
}

template<>
itk::Vector<double, 3>
VectorTypeConvertor<double, double>(const itk::Vector<double, 3> &sourceVector)
{
  return sourceVector;
}

template<>
itk::Vector<float, 3>
VectorTypeConvertor<float, float>(const itk::Vector<float, 3> &sourceVector)
{
  return sourceVector;
}

// converts double precision type versors to single precision type
template<class TInputScalarType, class TOutputScalarType>
itk::Versor<TOutputScalarType>
VersorTypeConvertor(const itk::Versor<TInputScalarType> &sourceVersor)
{
  itk::Versor<TOutputScalarType> OutputVersor;
  TOutputScalarType x_value, y_value, z_value, w_value;
  x_value = (TOutputScalarType)( sourceVersor.GetX() );
  y_value = (TOutputScalarType)( sourceVersor.GetY() );
  z_value = (TOutputScalarType)( sourceVersor.GetZ() );
  w_value = (TOutputScalarType)( sourceVersor.GetW() );
  OutputVersor.Set(x_value, y_value, z_value, w_value);
  return OutputVersor;
}

template<>
itk::Versor<double>
VersorTypeConvertor<double, double>(const itk::Versor<double> &sourceVersor)
{
  return sourceVersor;
}

template<>
itk::Versor<float>
VersorTypeConvertor<float, float>(const itk::Versor<float> &sourceVersor)
{
  return sourceVersor;
}

// converts double precision type matrices to single precision type
template<class TInputScalarType, class TOutputScalarType>
itk::Matrix<TOutputScalarType, 3, 3>
MatrixTypeConvertor(const itk::Matrix<TInputScalarType, 3, 3> &sourceMatix)
{
  itk::Matrix<TOutputScalarType, 3, 3> outputMatrix;
  for( itk::SizeValueType i=0; i<3; ++i )
    {
    for( itk::SizeValueType j=0; j<3; ++j )
      {
      outputMatrix(i,j) = (TOutputScalarType)( sourceMatix(i,j) );
      }
    }
  return outputMatrix;
}

template<>
itk::Matrix<double, 3, 3>
MatrixTypeConvertor<double, double>(const itk::Matrix<double, 3, 3> &sourceMatix)
{
  return sourceMatix;
}

template<>
itk::Matrix<float, 3, 3>
MatrixTypeConvertor<float, float>(const itk::Matrix<float, 3, 3> &sourceMatix)
{
  return sourceMatix;
}

// conversion between the precision type for two transforms that are in the same class.
template<class TInputScalarType, class TOutputScalarType>
bool
PrecisionConvertor( typename itk::Transform< TOutputScalarType, 3, 3 >::Pointer &outXfrm,
                   const itk::Transform< TInputScalarType, 3, 3 > *source )
{
  if( IsClass(source, "TranslationTransform") )
    {
    typedef itk::TranslationTransform<TOutputScalarType, 3> TransTransformType;
    typename TransTransformType::Pointer transXfrm = TransTransformType::New();
    outXfrm = transXfrm.GetPointer();
    }
  else if( IsClass(source, "AffineTransform") )
    {
    typedef itk::AffineTransform<TOutputScalarType, 3>  AffineTransformType;
    typename AffineTransformType::Pointer affineXfrm = AffineTransformType::New();
    outXfrm = affineXfrm.GetPointer();
    }
  else if( IsClass(source, "VersorTransform") )
    {
    typedef itk::VersorTransform<TOutputScalarType> VersorTransformType;
    typename VersorTransformType::Pointer versorXfrm = VersorTransformType::New();
    outXfrm = versorXfrm.GetPointer();
    }
  else if( IsClass(source, "VersorRigid3DTransform") )
    {
    typedef itk::VersorRigid3DTransform<TOutputScalarType>  VersorRigid3DTransformType;
    typename VersorRigid3DTransformType::Pointer versorRigid3DXfrm = VersorRigid3DTransformType::New();
    outXfrm = versorRigid3DXfrm.GetPointer();
    }
  else if( IsClass(source, "ScaleVersor3DTransform") )
    {
    typedef itk::ScaleVersor3DTransform<TOutputScalarType>  ScaleVersor3DTransformType;
    typename ScaleVersor3DTransformType::Pointer scaleVersor3DXfrm = ScaleVersor3DTransformType::New();
    outXfrm = scaleVersor3DXfrm.GetPointer();
    }
  else if( IsClass(source, "ScaleSkewVersor3DTransform") )
    {
    typedef itk::ScaleSkewVersor3DTransform<TOutputScalarType> ScaleVersor3DTransformType;
    typename ScaleVersor3DTransformType::Pointer scaleSkewVersor3DXfrm = ScaleVersor3DTransformType::New();
    outXfrm = scaleSkewVersor3DXfrm.GetPointer();
    }
  else
    {
    TransformConvertError<TInputScalarType>( source, "Different Precision Type");
    return false;
    }
  outXfrm->SetParameters( ParametersTypeConvertor<TInputScalarType, TOutputScalarType>( source->GetParameters() ) );
  outXfrm->SetFixedParameters( ParametersTypeConvertor<TInputScalarType, TOutputScalarType>( source->GetFixedParameters() ) );
  return true;
}

//
// Convert from any type derived from MatrixOffsetTransformType to
// AffineTransform.
template<class TInputScalarType, class TOutputScalarType>
bool
ExtractTransform(typename itk::AffineTransform< TOutputScalarType, 3 >::Pointer &result,
                 const itk::Transform< TInputScalarType, 3, 3 > *source)
{
  result->SetIdentity();
  // always able to convert to same type
  typedef itk::Transform< TOutputScalarType, 3, 3 > OutputGenericTransformType;
  typename OutputGenericTransformType::Pointer resultXfrm = dynamic_cast<OutputGenericTransformType *>( result.GetPointer() );
  if( IsSameClass(result.GetPointer(), source) )
    {
     if( PrecisionConvertor<TInputScalarType, TOutputScalarType>(resultXfrm, source) )
       {
       result = dynamic_cast<typename itk::AffineTransform< TOutputScalarType, 3 > *>(resultXfrm.GetPointer());
       return true;
       }
     else
       {
       std::cerr << "Can't convert the precision type of the input transform to the output precision type." << std::endl;
       return false;
       }
    }

  typedef itk::AffineTransform< TInputScalarType, 3 > InputAffineTransformType;
  typedef typename InputAffineTransformType::Superclass MatrixOffsetTransformType;
  const MatrixOffsetTransformType *matBasePtr = dynamic_cast<const MatrixOffsetTransformType *>(source);
  if( matBasePtr == 0 )
    {
    return false;
    }

  result->SetCenter( PointTypeConvertor<TInputScalarType, TOutputScalarType>( matBasePtr->GetCenter() ) );
  result->SetMatrix( MatrixTypeConvertor<TInputScalarType, TOutputScalarType>( matBasePtr->GetMatrix() ) );
  result->SetTranslation( VectorTypeConvertor<TInputScalarType, TOutputScalarType>( matBasePtr->GetTranslation() ) );
  return true;
}

//
// versor rigid 3d case.
template<class TInputScalarType, class TOutputScalarType>
bool
ExtractTransform(typename itk::VersorRigid3DTransform<TOutputScalarType>::Pointer & result,
                 const itk::Transform< TInputScalarType, 3, 3 > *source)
{
  result->SetIdentity();
  // always able to convert to same type
  typedef itk::Transform< TOutputScalarType, 3, 3 > OutputGenericTransformType;
  typename OutputGenericTransformType::Pointer resultXfrm = dynamic_cast<OutputGenericTransformType *>( result.GetPointer() );
  if( IsSameClass(result.GetPointer(), source) )
    {
    if( PrecisionConvertor<TInputScalarType, TOutputScalarType>(resultXfrm, source) )
      {
      result = dynamic_cast<typename itk::VersorRigid3DTransform< TOutputScalarType > *>(resultXfrm.GetPointer());
      return true;
      }
    else
      {
      std::cerr << "Can't convert the precision type of the input transform to the output precision type." << std::endl;
      return false;
      }
    }

  // this looks like it should be a convertible transform but
  // I'm not sure.
  typedef itk::TranslationTransform<TInputScalarType, 3> TransTransformType;
  if( IsClass(source, "TranslationTransform") )
    {
    const TransTransformType *translationXfrm = dynamic_cast<const TransTransformType *>(source);
    typename TransTransformType::OutputVectorType offset = translationXfrm->GetOffset();

    result->SetOffset( VectorTypeConvertor<TInputScalarType, TOutputScalarType>(offset) );
    return true;
    }
  // versor == rotation only
  if( IsClass(source, "VersorTransform") )
    {
    typedef itk::VersorTransform<TInputScalarType> VersorTransformType;
    const VersorTransformType *versorXfrm = dynamic_cast<const VersorTransformType *>(source);

    result->SetRotation( VersorTypeConvertor<TInputScalarType, TOutputScalarType>( versorXfrm->GetVersor() ) );
    result->SetCenter( PointTypeConvertor<TInputScalarType, TOutputScalarType>( versorXfrm->GetCenter() ) );
    return true;
    }
  return false;
}

//
// scale versor case
template<class TInputScalarType, class TOutputScalarType>
bool
ExtractTransform(typename itk::ScaleVersor3DTransform<TOutputScalarType>::Pointer & result,
                 const itk::Transform< TInputScalarType, 3, 3 > *source)
{
  result->SetIdentity();
  // always able to convert to same type
  typedef itk::Transform< TInputScalarType, 3, 3 >  InputGenericTransformType;
  typedef itk::Transform< TOutputScalarType, 3, 3 > OutputGenericTransformType;
  typename OutputGenericTransformType::Pointer resultXfrm = dynamic_cast<OutputGenericTransformType *>( result.GetPointer() );
  if( IsSameClass(result.GetPointer(), source) )
    {
    if( PrecisionConvertor<TInputScalarType, TOutputScalarType>(resultXfrm, source) )
      {
      result = dynamic_cast<typename itk::ScaleVersor3DTransform<TOutputScalarType> *>(resultXfrm.GetPointer());
      return true;
      }
    else
      {
      std::cerr << "Can't convert the precision type of the input transform to the output precision type." << std::endl;
      return false;
      }
    }

  typedef itk::VersorRigid3DTransform<TInputScalarType>  InputVersorRigid3DTransformType;
  typedef itk::VersorRigid3DTransform<TOutputScalarType> OutputVersorRigid3DTransformType;
  if( IsClass(source, "VersorRigid3DTransform") )
    {
    const InputVersorRigid3DTransformType *versorRigidXfrm = dynamic_cast<const InputVersorRigid3DTransformType *>(source);

    result->SetRotation( VersorTypeConvertor<TInputScalarType, TOutputScalarType>( versorRigidXfrm->GetVersor() ) );
    result->SetTranslation( VectorTypeConvertor<TInputScalarType, TOutputScalarType>( versorRigidXfrm->GetTranslation() ) );
    result->SetCenter( PointTypeConvertor<TInputScalarType, TOutputScalarType>( versorRigidXfrm->GetCenter() ) );
    return true;
    }
  // otherwise try VersorRigidTransform
  typename OutputVersorRigid3DTransformType::Pointer vrx_out = OutputVersorRigid3DTransformType::New();
  if( ExtractTransform<TInputScalarType, TOutputScalarType>(vrx_out, source) ) // of VersorRigid3D conversion
                                                                               // works
    {
    // Now "vrx_out" has the output precision type. It should be converted to input precision type (vrx_in) before
    // passing that again to "ExtractTransform" function.
    typename InputVersorRigid3DTransformType::Pointer vrx_in = InputVersorRigid3DTransformType::New();
    typename InputGenericTransformType::Pointer vrxInXfrm = dynamic_cast<InputGenericTransformType *>( vrx_in.GetPointer() );
    if( PrecisionConvertor<TOutputScalarType, TInputScalarType>( vrxInXfrm, vrx_out.GetPointer() ) ) // precision type conversion
      {
      vrx_in = dynamic_cast< InputVersorRigid3DTransformType *>(vrxInXfrm.GetPointer());
        // recurse to do this conversion
      return ExtractTransform<TInputScalarType, TOutputScalarType>( result, vrx_in.GetPointer() );
      }
    }
  return false;
}

//
// scale skew versor case
template<class TInputScalarType, class TOutputScalarType>
bool
ExtractTransform(typename itk::ScaleSkewVersor3DTransform< TOutputScalarType >::Pointer & result,
                 const itk::Transform< TInputScalarType, 3, 3 > *source)
{
  // always able to convert to same type
  typedef itk::Transform< TInputScalarType, 3, 3 >  InputGenericTransformType;
  typedef itk::Transform< TOutputScalarType, 3, 3 > OutputGenericTransformType;
  typename OutputGenericTransformType::Pointer resultXfrm = dynamic_cast<OutputGenericTransformType *>( result.GetPointer() );
  if( IsSameClass(result.GetPointer(), source) )
    {
    if( PrecisionConvertor<TInputScalarType, TOutputScalarType>(resultXfrm, source) )
      {
      result = dynamic_cast< typename itk::ScaleSkewVersor3DTransform< TOutputScalarType > *>(resultXfrm.GetPointer());
      return true;
      }
    else
      {
      std::cerr << "Can't convert the precision type of the input transform to the output precision type." << std::endl;
      return false;
      }
    }

  // is it the parent?
  typedef itk::ScaleVersor3DTransform<TInputScalarType>  InputScaleVersor3DTransformType;
  typedef itk::ScaleVersor3DTransform<TOutputScalarType> OutputScaleVersor3DTransformType;
  if( IsClass(source, "ScaleVersor3DTransform") )
    {
    const InputScaleVersor3DTransformType *scaleVersorXfrm = dynamic_cast<const InputScaleVersor3DTransformType *>(source);

    result->SetRotation( VersorTypeConvertor<TInputScalarType, TOutputScalarType>( scaleVersorXfrm->GetVersor() ) );
    result->SetTranslation( VectorTypeConvertor<TInputScalarType, TOutputScalarType>( scaleVersorXfrm->GetTranslation() ) );
    result->SetCenter( PointTypeConvertor<TInputScalarType, TOutputScalarType>( scaleVersorXfrm->GetCenter() ) );
    result->SetScale( VectorTypeConvertor<TInputScalarType, TOutputScalarType>( scaleVersorXfrm->GetScale() ) );
    return true;
    }
  // otherwise try ScaleVersor conversion
  typename OutputScaleVersor3DTransformType::Pointer svx_out = OutputScaleVersor3DTransformType::New();
  if( ExtractTransform<TInputScalarType, TOutputScalarType>(svx_out, source) ) // of VersorRigid3D conversion
                                                                               // works
    {
    // Now "svx_out" has the output precision type. It should be converted to input precision type (svx_in) before
    // passing that again to "ExtractTransform" function.
    typename InputScaleVersor3DTransformType::Pointer svx_in = InputScaleVersor3DTransformType::New();
    typename InputGenericTransformType::Pointer svxInXfrm = dynamic_cast<InputGenericTransformType *>( svx_in.GetPointer() );
    if( PrecisionConvertor<TOutputScalarType, TInputScalarType>( svxInXfrm, svx_out.GetPointer() ) )
      {
      svx_in = dynamic_cast< InputScaleVersor3DTransformType *>(svxInXfrm.GetPointer());
      // recurse to do this conversion
      return ExtractTransform<TInputScalarType, TOutputScalarType>( result, svx_in.GetPointer() );
      }
    }
  return false;
}

#define CHECK_PARAMETER_IS_SET(parameter, message) \
  if( parameter == "" )                             \
    {                                             \
    std::cerr << message << std::endl;            \
    return EXIT_FAILURE;                          \
    }

template<class TInputScalarType, class TOutputScalarType>
int
DoConversion( int argc, char *argv[] )
{
  PARSE_ARGS;

  typedef itk::Transform< TInputScalarType, 3, 3 >                                   InputGenericTransformType;
  typedef itk::BSplineDeformableTransform< TInputScalarType,
                                           GenericTransformImageNS::SpaceDimension,
                                           GenericTransformImageNS::SplineOrder>     BSplineTransformType;

  typedef itk::Transform< TOutputScalarType, 3, 3 >               OutputGenericTransformType;
  typedef itk::AffineTransform< TOutputScalarType, 3 >            AffineTransformType;
  typedef itk::VersorRigid3DTransform< TOutputScalarType >        VersorRigid3DTransformType;
  typedef itk::ScaleVersor3DTransform< TOutputScalarType >        ScaleVersor3DTransformType;
  typedef itk::ScaleSkewVersor3DTransform< TOutputScalarType >    ScaleSkewVersor3DTransformType;

  // read the input transform
  typedef itk::TransformFileReaderTemplate<TInputScalarType>  TransformFileReaderType;
  typename TransformFileReaderType::Pointer reader = TransformFileReaderType::New();
  reader->SetFileName(inputTransform.c_str() );
  reader->Update();
  typename TransformFileReaderType::TransformListType *transformList = reader->GetTransformList();
  typename InputGenericTransformType::Pointer inputXfrm = dynamic_cast<InputGenericTransformType *>( transformList->front().GetPointer() );

  std::cout << "------------------------ " << std::endl;
  std::cout << "Input Transform Type ==> " << inputXfrm->GetTransformTypeAsString() << std::endl;
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
  typename OutputGenericTransformType::Pointer outputXfrm;

  if( outputTransformType == "Affine" )
    {
    typename AffineTransformType::Pointer affineXfrm = AffineTransformType::New();
    if( ExtractTransform<TInputScalarType, TOutputScalarType>(affineXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TInputScalarType>(inputXfrm, "Affine Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = affineXfrm.GetPointer();
    }
  else if( outputTransformType == "VersorRigid" )
    {
    typename VersorRigid3DTransformType::Pointer versorRigidXfrm = VersorRigid3DTransformType::New();
    if( ExtractTransform<TInputScalarType, TOutputScalarType>(versorRigidXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TInputScalarType>(inputXfrm, "VersorRigid3D Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = versorRigidXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleVersor" )
    {
    typename ScaleVersor3DTransformType::Pointer scaleVersorXfrm = ScaleVersor3DTransformType::New();
    if( ExtractTransform<TInputScalarType, TOutputScalarType>( scaleVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TInputScalarType>(inputXfrm, "ScaleVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleVersorXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleSkewVersor" )
    {
    typename ScaleSkewVersor3DTransformType::Pointer scaleSkewVersorXfrm = ScaleSkewVersor3DTransformType::New();
    if( ExtractTransform<TInputScalarType, TOutputScalarType>( scaleSkewVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TInputScalarType>(inputXfrm, "ScaleSkewVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleSkewVersorXfrm.GetPointer();
    }

  if( outputTransformType == "Same" )
    {
    typedef typename itk::TransformFileWriterTemplate<TOutputScalarType> TransformWriterType;
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName(outputTransform);
    for( typename itk::TransformFileReaderTemplate<TInputScalarType>::TransformListType::iterator it = transformList->begin();
         it != transformList->end(); ++it )
      {
      typename InputGenericTransformType::Pointer inXfrm = dynamic_cast<InputGenericTransformType *>( (*it).GetPointer() );
      typename OutputGenericTransformType::Pointer outXfrm;
      if( PrecisionConvertor<TInputScalarType, TOutputScalarType>(outXfrm, inXfrm) )
        {
        transformWriter->AddTransform( outXfrm );
        //
        std::cout << "Output Transform Type ==> " << outXfrm->GetTransformTypeAsString() << std::endl;
        std::cout << "* Output transform parameters: " << outXfrm->GetParameters() << std::endl;
        std::cout << "* Output transform fixed parameters: " << outXfrm->GetFixedParameters() << std::endl;
        std::cout << "------------------------ " << std::endl;
        }
      else
        {
        std::cerr << "Can't convert the input transform to the same transform but different precision type." << std::endl;
        return EXIT_FAILURE;
        }
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
    std::cout << "Output Transform Type ==> " << outputXfrm->GetTransformTypeAsString() << std::endl;
    std::cout << "* Output transform parameters: " << outputXfrm->GetParameters() << std::endl;
    std::cout << "* Output transform fixed parameters: " << outputXfrm->GetFixedParameters() << std::endl;
    std::cout << "------------------------ " << std::endl;
    //
    itk::WriteTransformToDisk<TOutputScalarType>(outputXfrm.GetPointer(), outputTransform);
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
  CHECK_PARAMETER_IS_SET(inputPrecisionType,
                         "Missing inputPrecisionType");
  CHECK_PARAMETER_IS_SET(outputPrecisionType,
                         "Missing outputPrecisionType");

  if( inputPrecisionType == "double" && outputPrecisionType == "double" )
    {
    return DoConversion<double, double>( argc, argv );
    }
  else if( inputPrecisionType == "double" && outputPrecisionType == "float" )
    {
    return DoConversion<double, float>( argc, argv );
    }
  else if( inputPrecisionType == "float" && outputPrecisionType == "double" )
    {
    return DoConversion<float, double>( argc, argv );
    }
  else if( inputPrecisionType == "float" && outputPrecisionType == "float" )
    {
    return DoConversion<float, float>( argc, argv );
    }
  else
    {
    std::cerr << "Error: Invalid parameters for input and output precision type." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
