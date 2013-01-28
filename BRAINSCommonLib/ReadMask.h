#ifndef __ReadMask_h
#define __ReadMask_h

#include "itkIO.h"
#include "itkImageMaskSpatialObject.h"

template <class MaskType, unsigned VDimension>
typename MaskType::Pointer
ReadImageMask(const std::string & filename,
              typename itk::ImageBase<VDimension> * /*referenceImage*/)
{
  typedef unsigned char                MaskPixelType;
  typedef itk::Image<MaskPixelType, 3> ReadMaskImageType;

  typename ReadMaskImageType::Pointer OrientedMaskImage = itkUtil::ReadImage<ReadMaskImageType>(filename);
  // TODO:  May want to check that physical spaces overlap?

  // convert mask image to mask
  typedef itk::ImageMaskSpatialObject<ReadMaskImageType::ImageDimension> ReadImageMaskSpatialObjectType;
  typename ReadImageMaskSpatialObjectType::Pointer mask = ReadImageMaskSpatialObjectType::New();
  mask->SetImage(OrientedMaskImage);
  //
  mask->ComputeObjectToWorldTransform();
  // return pointer to mask
  typename MaskType::Pointer p = dynamic_cast<MaskType *>( mask.GetPointer() );
  if( p.IsNull() )
    {
    itkGenericExceptionMacro(<< "Failed conversion to Mask");
    }
  return p;
}

#endif // LoadMask_h
