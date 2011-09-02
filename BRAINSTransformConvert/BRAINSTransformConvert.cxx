#include "BRAINSTransformConvertCLP.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageFileReader.h"
#include "itkBSPlineDeformableTransform.h"
#include "itkIO.h"
#include "itkImageRegionIterator.h"

#define CHECK_PARAMETER_IS_SET(parameter, message) \
  if( parameter == "" )                             \
    {                                             \
    std::cerr << message << std::endl;            \
    return EXIT_FAILURE;                          \
    }

int main(int argc, char *argv[])
{
  typedef itk::Transform<double>                        TransformType;
  typedef itk::BSplineDeformableTransform<double, 3, 3> BSplineDeformableTransformType;

  PARSE_ARGS;

  CHECK_PARAMETER_IS_SET(inputTransform,
                         "Missing inputTransform parameters");
  CHECK_PARAMETER_IS_SET(outputTransformType,
                         "Missing outpuTransformType");

  // read the input transform
  itk::TransformFileReader::Pointer reader =
    itk::TransformFileReader::New();
  reader->SetFileName(inputTransform.c_str() );
  reader->Update();

  if( outputTransformType == "DisplacementField" )
    {
    CHECK_PARAMETER_IS_SET(referenceVolume,
                           "Missing referenceVolume needed for Displacement Field output");
    CHECK_PARAMETER_IS_SET(displacementVolume,
                           "Missing displacementVolume needed for Displacement Field output");

    itk::TransformFileReader::TransformListType *transformList =
      reader->GetTransformList();
    TransformType::Pointer xfrm = dynamic_cast<TransformType *>(transformList->front().GetPointer() );
    // Handle BSpline type
    BSplineDeformableTransformType::Pointer bsplineXfrm =
      dynamic_cast<BSplineDeformableTransformType *>(xfrm.GetPointer() );

    if( bsplineXfrm.IsNotNull() )
      {
      transformList->pop_front();
      if( transformList->size() == 0 )
        {
        std::cerr << "Error, the second transform needed for BSplineDeformableTransform is missing." << std::endl;
        return EXIT_FAILURE;
        }
      BSplineDeformableTransformType::BulkTransformType::Pointer bulkXfrm =
        dynamic_cast<BSplineDeformableTransformType::BulkTransformType *>(transformList->front().GetPointer() );
      if( bulkXfrm.IsNull() )
        {
        std::cerr << "Error, the second transform is not a bulk transform" << std::endl;
        }
      bsplineXfrm->SetBulkTransform(bulkXfrm);
      xfrm = bsplineXfrm.GetPointer();
      }

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
      movingPoint = xfrm->TransformPoint(fixedPoint);
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
    ;
    }
  else if( outputTransformType == "Affine" )
    {
    }
  else if( outputTransformType == "VersorRigid" )
    {
    }
  else if( outputTransformType == "ScaleVersor" )
    {
    }
  else if( outputTransformType == "ScaleSkewVersor" )
    {
    }
  else if( outputTransformType == "BSPlineDeformable" )
    {
    }
  return EXIT_SUCCESS;
}
