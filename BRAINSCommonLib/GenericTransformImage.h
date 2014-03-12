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
#ifndef _GenericTransformImage_H_
#define _GenericTransformImage_H_

#include "BRAINSCommonLibWin32Header.h"
#include <iostream>
#include "itkMacro.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkScaleVersor3DTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkAffineTransform.h"
#include <itkBSplineDeformableTransform.h>
#include <itkThinPlateR2LogRSplineKernelTransform.h>
#include "itkVersorRigid3DTransform.h"
#if (ITK_VERSION_MAJOR > 3)
#include "itkCompositeTransform.h"
#endif
#include "ConvertToRigidAffine.h"
#include "itkResampleImageFilter.h"
#include "itkImageDuplicator.h"
#include "Imgmath.h"

#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
//
// TODO:  The next two should be hidden in the cxx files again.
typedef itk::TransformFileReader                    TransformReaderType;
typedef itk::TransformFileReader::TransformListType TransformListType;

// TODO:  This should really be taken as a typedef from the BSpline class.
// TODO:  These should be hidden in the BRAINSFit namespace.
typedef itk::Transform<double, 3, 3> GenericTransformType;

namespace GenericTransformImageNS
{
static const unsigned int SpaceDimension = 3;
static const unsigned int SplineOrder = 3;
}

typedef double CoordinateRepType;
typedef itk::BSplineDeformableTransform<
    CoordinateRepType,
    GenericTransformImageNS::SpaceDimension,
    GenericTransformImageNS::SplineOrder> BSplineTransformType;

typedef itk::AffineTransform<double, 3>                      AffineTransformType;
typedef itk::VersorRigid3DTransform<double>                  VersorRigid3DTransformType;
typedef itk::ScaleVersor3DTransform<double>                  ScaleVersor3DTransformType;
typedef itk::ScaleSkewVersor3DTransform<double>              ScaleSkewVersor3DTransformType;
typedef itk::ThinPlateR2LogRSplineKernelTransform<double, 3> ThinPlateSpline3DTransformType;

#if (ITK_VERSION_MAJOR > 3)
typedef itk::CompositeTransform<double, 3> BRAINSCompositeTransformType;
#endif

namespace itk
{
/**
  * \author Hans J. Johnson
  * \brief A utility function to write ITK compliant transforms to disk in a way
  *that is compliant with the ReadTransformFromDisk
  * \param genericTransformToWrite A pointer to baseclass
  *itk::Transform<double,3,3> that is
  * polymorphically cast to a real class like (i.e. itk::VersorRigid,
  *itk::Affine, itk::BSpline, or similar)
  * \param outputTransform the filename of the output transform.
  * \code
  * //To convert from non-const smart pointers ploymorphically to the smart
  *ConstPointer up the class tree, use the GetPointer
  * AffineTransformType::Pointer myAffine=AffineTransformType::New(); //NOTE:
  * This is not a const smart pointer
  * WriteTransformToDisk<TScalarType>(myAffine.GetPointer(), "myAffineFile.mat");
  * \endcode
  */
template<class TScalarType>
extern void WriteTransformToDisk( itk::Transform<TScalarType, 3, 3> const *const genericTransformToWrite,
                                 const std::string & outputTransform);

/**
  * \author Hans J. Johnson
  * \brief A utility function to read ITK compliant transforms to disk in a way
  *that is compliant with the WriteTransformFromDisk
  * \param outputTransform the filename of the output transform.
  * \return A pointer to baseclass itk::Transform<double,3,3> that is
  * polymorphically convertable to a real class like (i.e. itk::VersorRigid,
  *itk::Affine, itk::BSpline, or similar)
  * \code
  * //To convert from non-const smart pointers ploymorphically to the smart
  *ConstPointer up the class tree, use the GetPointer
  * GenericTransformType::Pointer
  *myGenericTransform=ReadTransformFromDisk(myAffine.GetPointer(),
  *"myAffineFile.mat");
  *
  * VersorRigid3DTransformType::Pointer myVersorRigid3D=NULL;
  * {
  * const std::string transformFileType = myGenericTransform->GetNameOfClass();
  * if ( transformFileType == "VersorRigid3DTransform" )
  *   {
  *   myVersorRigid3D->SetParameters( versorRigid->GetParameters() );
  *   myVersorRigid3D->SetFixedParameters( versorRigid->GetFixedParameters() );
  *   }
  *   NOTE: It is more safe to copy parameters into the concrete class rather
  *than attempting to dynamically
  *   cast the base classes.  The reason is that the smart pointer management
  *becomes very unweildy and
  *   is hard to keep straight between the pointer that may delete the base
  *class, and the pointer that
  *   is the derived class.
  * }
  * \endcode
  */
extern GenericTransformType::Pointer ReadTransformFromDisk(const std::string & initialTransform);

/**
  * \author Hans J. Johnson
  * \brief A utility function to write ITK compliant transforms to disk in a way
  *that is compliant with the ReadTransformFromDisk
  * \param genericTransformToWrite A pointer to baseclass
  *itk::Transform<double,3,3> that is
  * polymorphically cast to a real class like (i.e. itk::VersorRigid,
  *itk::Affine, itk::BSpline, or similar)
  * \param outputTransform the filename of the output transform.
  * \code
  * //To convert from non-const smart pointers ploymorphically to the smart
  *ConstPointer up the class tree, use the GetPointer
  * AffineTransformType::Pointer myAffine=AffineTransformType::New(); //NOTE:
  * This is not a const smart pointer
  * WriteTransformToDisk<TScalarType>(myAffine.GetPointer(), "myAffineFile.mat");
  * \endcode
  */
extern VersorRigid3DTransformType::Pointer ComputeRigidTransformFromGeneric(
  const GenericTransformType::ConstPointer genericTransformToWrite);

/**
  * \author Hans J. Johnson
  * \brief Special purpose convenience function -- should not have a public
  *interface.
  */
extern int WriteBothTransformsToDisk(const GenericTransformType::ConstPointer genericTransformToWrite,
                                     const std::string & outputTransform, const std::string & strippedOutputTransform);

/**
  * \author Hans J. Johnson
  * \brief Special purpose convenience function -- should not have a public
  *interface.
  */
extern int WriteStrippedRigidTransformToDisk(const GenericTransformType::ConstPointer genericTransformToWrite,
                                             const std::string & strippedOutputTransform);
}

/**
  * \author Hans J. Johnson
  * \brief A class to transform images
  */
template <class InputImageType, class OutputImageType>
typename OutputImageType::Pointer
TransformResample(
  typename InputImageType::ConstPointer inputImage,
  typename itk::ImageBase<InputImageType::ImageDimension>::ConstPointer ReferenceImage,
  const typename InputImageType::PixelType defaultValue,
  typename itk::InterpolateImageFunction<InputImageType,
  typename itk::NumericTraits<typename InputImageType::PixelType>::RealType>::Pointer interp,
  typename GenericTransformType::ConstPointer transform);

/**
  * \author Hans J. Johnson
  * \brief A class to transform images
  */
template <class InputImageType, class OutputImageType, class DisplacementImageType>
typename OutputImageType::Pointer
  TransformWarp(
  InputImageType const *const inputImage,
  const itk::ImageBase<InputImageType::ImageDimension> *ReferenceImage,
  typename InputImageType::PixelType defaultValue,
  typename itk::InterpolateImageFunction<InputImageType,
                                         typename itk::NumericTraits<typename InputImageType::PixelType>::RealType>
  ::Pointer interp,
  typename DisplacementImageType::Pointer displacementField);

/**
  * \author Hans J. Johnson
  * \brief A class to transform images.  Only one of genericTransform or
  *DisplacementField can be non-null.
  */
template <typename InputImageType, class OutputImageType, typename DisplacementImageType>
typename OutputImageType::Pointer GenericTransformImage(
  InputImageType const *const OperandImage,
  const itk::ImageBase<InputImageType::ImageDimension> *ReferenceImage,
  // typename DisplacementImageType::Pointer DisplacementField,
  typename GenericTransformType::Pointer genericTransform,
  typename InputImageType::PixelType suggestedDefaultValue, // NOTE:  This is
                                                            // ignored in the
                                                            // case of binary
                                                            // image!
  const std::string & interpolationMode,
  const bool binaryFlag);

#ifndef ITK_MANUAL_INSTANTIATION
#include "GenericTransformImage.hxx"
#endif

#endif
