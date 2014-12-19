/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "BRAINSTransformConvertCLP.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageFileReader.h"
#include "itkBSplineDeformableTransform.h"
#include "itkIO.h"
#include "itkImageRegionIterator.h"
#include "GenericTransformImage.h"
#include "itkTranslationTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkComposeDisplacementFieldsImageFilter.h"

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

  typedef itk::AffineTransform< TScalarType, 3 > LocalAffineTransformType;
  typedef typename LocalAffineTransformType::Superclass MatrixOffsetTransformType;
  const MatrixOffsetTransformType *matBasePtr = dynamic_cast<const MatrixOffsetTransformType *>(source);
  if( matBasePtr == ITK_NULLPTR )
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

  typedef itk::VersorRigid3DTransform<TScalarType> LocalVersorRigid3DTransformType;
  if( IsClass(source, "VersorRigid3DTransform") )
    {
    const LocalVersorRigid3DTransformType *versorRigidXfrm =
      dynamic_cast<const LocalVersorRigid3DTransformType *>(source);
    result->SetCenter(versorRigidXfrm->GetCenter() );
    result->SetRotation(versorRigidXfrm->GetVersor() );
    result->SetTranslation(versorRigidXfrm->GetTranslation() );
    return true;
    }
  // otherwise try VersorRigidTransform
  typename LocalVersorRigid3DTransformType::Pointer vrx = LocalVersorRigid3DTransformType::New();
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
  typedef itk::ScaleVersor3DTransform<TScalarType> LocalScaleVersor3DTransformType;
  if( IsClass(source, "ScaleVersor3DTransform") )
    {
    const LocalScaleVersor3DTransformType *scaleVersorXfrm =
      dynamic_cast<const LocalScaleVersor3DTransformType *>(source);
    result->SetCenter(scaleVersorXfrm->GetCenter() );
    result->SetRotation(scaleVersorXfrm->GetVersor() );
    result->SetTranslation(scaleVersorXfrm->GetTranslation() );
    result->SetScale(scaleVersorXfrm->GetScale() );
    return true;
    }
  // otherwise try ScaleVersor conversion
  typename LocalScaleVersor3DTransformType::Pointer svx = LocalScaleVersor3DTransformType::New();
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
  BRAINSRegisterAlternateIO();

  typedef itk::Transform< TScalarType, 3, 3 >                       GenericTransformType;
  typedef itk::BSplineDeformableTransform < TScalarType, 3, 3>      BSplineTransformType;

  typedef itk::AffineTransform< TScalarType, 3 >            LocalAffineTransformTYpe;
  typedef itk::VersorRigid3DTransform< TScalarType >        LocalVersorRigid3DTransformType;
  typedef itk::ScaleVersor3DTransform< TScalarType >        LocalScaleVersor3DTransformType;
  typedef itk::ScaleSkewVersor3DTransform< TScalarType >    LocalScaleSkewVersor3DTransformType;

  // read the input transform
  typedef itk::TransformFileReaderTemplate<TScalarType>  TransformFileReaderType;
  typename TransformFileReaderType::Pointer reader = TransformFileReaderType::New();
  reader->SetFileName(inputTransform.c_str() );
  reader->Update();
  typename TransformFileReaderType::TransformListType *transformList = reader->GetTransformList();
  typename GenericTransformType::Pointer inputXfrm = dynamic_cast<GenericTransformType *>( transformList->front().GetPointer() );

  const std::string inputTransformTypeName = inputXfrm->GetTransformTypeAsString();
  std::cout << "------------------------ " << std::endl;
  std::cout << "Input Transform Type Saved in Memory ==> " << inputTransformTypeName << std::endl;
  if( inputTransformTypeName.find("CompositeTransform") == std::string::npos )
    {
    std::cout << "* Input transform parameters: " << inputXfrm->GetParameters() << std::endl;
    std::cout << "* Input transform fixed parameters: " << inputXfrm->GetFixedParameters() << std::endl;
    }
  std::cout << "------------------------ " << std::endl;

  // Handle BSpline type
  typename BSplineTransformType::Pointer bsplineInputXfrm =
    dynamic_cast<BSplineTransformType *>( inputXfrm.GetPointer() );
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

  if( inputTransformTypeName.find("CompositeTransform") != std::string::npos )
    {
    if( outputTransformType == "Same" )
      {
      typedef itk::CompositeTransform<TScalarType,3>                                      CompositeTransformType;
      typedef itk::DisplacementFieldTransform<TScalarType, 3>                             DisplacementFieldTransformType;
      typedef typename DisplacementFieldTransformType::DisplacementFieldType              DisplacementFieldType;
      typedef typename itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType,
                                                                 DisplacementFieldType>   ComposerType;

      typename CompositeTransformType::Pointer compToWrite;

      typename CompositeTransformType::ConstPointer compXfrm =
        dynamic_cast<const CompositeTransformType *>( inputXfrm.GetPointer() );
      if( compXfrm.IsNotNull() )
        {
        compToWrite = compXfrm->Clone();

        // If the last four transforms are displacementFieldType, we assume that they are
        // forward and inverse displacement fields of SyN internal transforms (FixedToMiddle and MovingToMiddle transforms),
        // so we compose them into one SyN displacementFieldTransformType.
        //
        unsigned int numOfTransforms = compToWrite->GetNumberOfTransforms();
        if( (compToWrite->GetNthTransform( numOfTransforms-1 )->GetTransformCategory() == GenericTransformType::DisplacementField)
           && (compToWrite->GetNthTransform( numOfTransforms-2 )->GetTransformCategory() == GenericTransformType::DisplacementField)
           && (compToWrite->GetNthTransform( numOfTransforms-3 )->GetTransformCategory() == GenericTransformType::DisplacementField)
           && (compToWrite->GetNthTransform( numOfTransforms-4 )->GetTransformCategory() == GenericTransformType::DisplacementField) )
          {
          typename DisplacementFieldTransformType::Pointer fixedToMiddleForwardTx =
            dynamic_cast<DisplacementFieldTransformType *>( compToWrite->GetNthTransform( numOfTransforms-4 ).GetPointer() );
          typename DisplacementFieldTransformType::Pointer fixedToMiddleInverseTx =
            dynamic_cast<DisplacementFieldTransformType *>( compToWrite->GetNthTransform( numOfTransforms-3 ).GetPointer() );
          typename DisplacementFieldTransformType::Pointer movingToMiddleForwardTx =
            dynamic_cast<DisplacementFieldTransformType *>( compToWrite->GetNthTransform( numOfTransforms-2 ).GetPointer() );
          typename DisplacementFieldTransformType::Pointer movingToMiddleInverseTx =
            dynamic_cast<DisplacementFieldTransformType *>( compToWrite->GetNthTransform( numOfTransforms-1 ).GetPointer() );

          typename DisplacementFieldTransformType::Pointer fixedToMiddleTransform = DisplacementFieldTransformType::New();
          fixedToMiddleTransform->SetDisplacementField( fixedToMiddleForwardTx->GetDisplacementField() );
          fixedToMiddleTransform->SetInverseDisplacementField( fixedToMiddleInverseTx->GetDisplacementField() );

          typename DisplacementFieldTransformType::Pointer movingToMiddleTransform = DisplacementFieldTransformType::New();
          movingToMiddleTransform->SetDisplacementField( movingToMiddleForwardTx->GetDisplacementField() );
          movingToMiddleTransform->SetInverseDisplacementField( movingToMiddleInverseTx->GetDisplacementField() );

          typename DisplacementFieldTransformType::Pointer resultSyNTransform = DisplacementFieldTransformType::New();

          typename ComposerType::Pointer composer = ComposerType::New();
          composer->SetDisplacementField( movingToMiddleTransform->GetInverseDisplacementField() );
          composer->SetWarpingField( fixedToMiddleTransform->GetDisplacementField() );
          composer->Update();

          typename ComposerType::Pointer inverseComposer = ComposerType::New();
          inverseComposer->SetDisplacementField( fixedToMiddleTransform->GetInverseDisplacementField() );
          inverseComposer->SetWarpingField( movingToMiddleTransform->GetDisplacementField() );
          inverseComposer->Update();

          resultSyNTransform->SetDisplacementField( composer->GetOutput() );
          resultSyNTransform->SetInverseDisplacementField( inverseComposer->GetOutput() );

          // First remove the last four displacement field transform related to the internal states
          compToWrite->RemoveTransform();
          compToWrite->RemoveTransform();
          compToWrite->RemoveTransform();
          compToWrite->RemoveTransform();
          // Then add the restored SyN transform
          compToWrite->AddTransform( resultSyNTransform );

          std::string transfromToWriteName = outputTransform + "Composite.h5";

          // This function writes both forward and inverse composite transforms to the disk
          std::cout << "Converted Transforms are Written to the Disk ==> " << compToWrite->GetTransformTypeAsString() << std::endl;
          itk::WriteTransformToDisk<TScalarType>( compToWrite.GetPointer(), transfromToWriteName );
          }
        }
      else
        {
        return EXIT_FAILURE;
        }
      }
    else
      {
      std::cerr << "Input transform is a CompositeTransform! Output transform must be the \"Same\" type!" << std::endl;
      return EXIT_FAILURE;
      }
    }

  //
  // if no transform name given, don't write transform
  if(outputTransform.size() == 0)
    {
    return EXIT_SUCCESS;
    }

  //output transform processing
  typename GenericTransformType::Pointer outputXfrm;

  if( outputTransformType == "Affine" )
    {
    typename LocalAffineTransformTYpe::Pointer affineXfrm = LocalAffineTransformTYpe::New();
    if( ExtractTransform<TScalarType>(affineXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "Affine Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = affineXfrm.GetPointer();
    }
  else if( outputTransformType == "VersorRigid" )
    {
    typename LocalVersorRigid3DTransformType::Pointer versorRigidXfrm = LocalVersorRigid3DTransformType::New();
    if( ExtractTransform<TScalarType>(versorRigidXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "VersorRigid3D Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = versorRigidXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleVersor" )
    {
    typename LocalScaleVersor3DTransformType::Pointer scaleVersorXfrm = LocalScaleVersor3DTransformType::New();
    if( ExtractTransform<TScalarType>( scaleVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "ScaleVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleVersorXfrm.GetPointer();
    }
  else if( outputTransformType == "ScaleSkewVersor" )
    {
    typename LocalScaleSkewVersor3DTransformType::Pointer scaleSkewVersorXfrm = LocalScaleSkewVersor3DTransformType::New();
    if( ExtractTransform<TScalarType>( scaleSkewVersorXfrm, inputXfrm.GetPointer() ) == false )
      {
      TransformConvertError<TScalarType>(inputXfrm, "ScaleSkewVersor Transform");
      return EXIT_FAILURE;
      }
    outputXfrm = scaleSkewVersorXfrm.GetPointer();
    }

  if( inputTransformTypeName.find("CompositeTransform") == std::string::npos ) // Input composite is assumed to be a state file
                                                                               // not a transform, so it is converted differently
                                                                               // and is already written to the disk.
    {
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
      //
      itk::WriteTransformToDisk<TScalarType>(outputXfrm.GetPointer(), outputTransform);
      }
    }

  return EXIT_SUCCESS;
}


int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  CHECK_PARAMETER_IS_SET(inputTransform,
                         "Missing inputTransform parameter");
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
