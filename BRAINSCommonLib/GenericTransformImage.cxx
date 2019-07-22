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
#include <iostream>
#include "GenericTransformImage.h"

#include "itkAffineTransform.h"
#include "itkCenteredAffineTransform.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include "itkCenteredSimilarity2DTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkFixedCenterOfRotationAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkRigid2DTransform.h"
#include "itkRigid3DPerspectiveTransform.h"
#include "itkScalableAffineTransform.h"
#include "itkScaleLogarithmicTransform.h"
#include "itkScaleSkewVersor3DTransform.h"
#include "itkScaleTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkVersorTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkBSplineTransform.h"
#include "itkTransformFactory.h"
#include <itksys/SystemTools.hxx>

// #include "itkSimilarity2DTransfor3DPerspectiveTransform.h"

namespace itk
{
itk::VersorRigid3DTransform< double >::Pointer
ComputeRigidTransformFromGeneric( const itk::Transform< double, 3, 3 >::ConstPointer genericTransformToWrite )
{
  using VersorRigid3DTransformType = itk::VersorRigid3DTransform< double >;
  using ScaleVersor3DTransformType = itk::ScaleVersor3DTransform< double >;
  using ScaleSkewVersor3DTransformType = itk::ScaleSkewVersor3DTransform< double >;

  VersorRigid3DTransformType::Pointer versorRigid = VersorRigid3DTransformType::New();
  versorRigid->SetIdentity();
  // //////////////////////////////////////////////////////////////////////////
  // ConvertTransforms
  if ( genericTransformToWrite.IsNotNull() )
  {
    try
    {
      const std::string transformFileType = genericTransformToWrite->GetNameOfClass();
      if ( transformFileType == "VersorRigid3DTransform" )
      {
        const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< VersorRigid3DTransformType const * >( genericTransformToWrite.GetPointer() );
        AssignRigid::ExtractVersorRigid3DTransform( versorRigid, tempInitializerITKTransform );
      }
      else if ( transformFileType == "ScaleVersor3DTransform" )
      {
        const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< ScaleVersor3DTransformType const * >( genericTransformToWrite.GetPointer() );
        AssignRigid::ExtractVersorRigid3DTransform( versorRigid, tempInitializerITKTransform );
      }
      else if ( transformFileType == "ScaleSkewVersor3DTransform" )
      {
        const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< ScaleSkewVersor3DTransformType const * >( genericTransformToWrite.GetPointer() );
        AssignRigid::ExtractVersorRigid3DTransform( versorRigid, tempInitializerITKTransform );
      }
      else if ( transformFileType == "AffineTransform" )
      {
        using AffineTransformType = itk::AffineTransform< double, 3 >;
        const AffineTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< AffineTransformType const * >( genericTransformToWrite.GetPointer() );
        AssignRigid::ExtractVersorRigid3DTransform( versorRigid, tempInitializerITKTransform );
      }
      else //  NO SUCH CASE!!
      {
        std::cout << "Compute Rigid transform from generic: "
                  << "Unsupported input transform file -- first transform typestring, " << transformFileType
                  << " not equal to any recognized type VersorRigid3DTransform OR "
                  << " ScaleVersor3DTransform OR ScaleSkewVersor3DTransform OR AffineTransform" << std::endl;
        return nullptr;
      }
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cout << "[FAILED]" << std::endl;
      std::cerr << "Error while converting the genericTransformToWrite to Rigid" << std::endl;
      std::cerr << excp << std::endl;
      return nullptr;
    }
  }
  return versorRigid;
}

template < typename TInputScalarType, typename TWriteScalarType >
int
WriteBothTransformsToDisk(
  const typename itk::Transform< TInputScalarType, 3, 3 >::ConstPointer genericTransformToWrite,
  const std::string & outputTransform, const std::string & strippedOutputTransform )
{
  using VersorRigid3DTransformType = itk::VersorRigid3DTransform< TInputScalarType >;
  using CompositeTransformType = itk::CompositeTransform< TInputScalarType, 3 >;
  // //////////////////////////////////////////////////////////////////////////
  // Write output transforms of BRAINSFit.
  /*
   * the input of genericTransformToWrite from the BRIANSFit is a composite transform.
   * If this composite transform has just one component, we write the included transform.
   * If the composite transform has more than one component i.e {a linear transform; a BSpline or SyN transform},
   * we write the composite transform itself.
   */

  if ( genericTransformToWrite.IsNull() )
  {
    return 0;
  }
  if ( std::string( genericTransformToWrite->GetNameOfClass() ) != "CompositeTransform" )
  {
    itkGenericExceptionMacro( << "Error in type conversion" );
  }
  const typename CompositeTransformType::ConstPointer genericCompositeTransform =
    static_cast< const CompositeTransformType * >( genericTransformToWrite.GetPointer() );
  if ( genericCompositeTransform->GetNumberOfTransforms() > 1 )
  {
    std::cout << "Write the output composite transform to the disk ..." << std::endl;
    const std::string extension = itksys::SystemTools::GetFilenameLastExtension( outputTransform );
    std::string       compositeTransformName( outputTransform );
    if ( extension != ".h5" && extension != ".hd5" )
    {
      compositeTransformName = itksys::SystemTools::GetFilenameWithoutExtension( outputTransform ) + "Composite.h5";
      std::cerr << "Warning: Composite transforms should always be HDF files" << std::endl
                << "Changing filename from " << outputTransform << " to " << compositeTransformName << std::endl;
    }
    itk::WriteTransformToDisk< TInputScalarType, TWriteScalarType >( genericCompositeTransform.GetPointer(),
                                                                     compositeTransformName.c_str() );
  }
  else
  {
    typename CompositeTransformType::TransformTypePointer genericComponent =
      genericCompositeTransform->GetNthTransform( 0 );
    try
    {
      const std::string transformFileType = genericComponent->GetNameOfClass();
      if ( transformFileType == "VersorRigid3DTransform" || transformFileType == "ScaleVersor3DTransform" ||
           transformFileType == "ScaleSkewVersor3DTransform" || transformFileType == "AffineTransform" )
      {
        if ( outputTransform.size() > 0 ) // Write out the transform
        {
          itk::WriteTransformToDisk< TInputScalarType, TWriteScalarType >( genericComponent, outputTransform );
        }
      }
      else if ( transformFileType == "BSplineTransform" )
      {
        using BSplineTransformType = itk::BSplineTransform< TInputScalarType, 3, 3 >;

        const typename BSplineTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< BSplineTransformType const * >( genericComponent.GetPointer() );
        if ( strippedOutputTransform.size() > 0 )
        {
          std::cout << "ERROR:  The rigid component of a BSpline transform is not supported." << std::endl;
        }
        if ( outputTransform.size() > 0 )
        {
          itk::WriteTransformToDisk< TInputScalarType, TWriteScalarType >( tempInitializerITKTransform,
                                                                           outputTransform );
        }
      }
      else if ( transformFileType == "displacementFieldTransform" )
      {
        if ( outputTransform.size() > 0 ) // Write out the transform
        {
          itk::WriteTransformToDisk< TInputScalarType, TWriteScalarType >( genericComponent, outputTransform );
        }
      }
      else //  NO SUCH CASE!!
      {
        std::cout << "Unsupported transform file -- " << transformFileType
                  << " not equal to any recognized type VersorRigid3DTransform OR"
                  << " ScaleVersor3DTransform OR ScaleSkewVersor3DTransform OR AffineTransform OR BSplineTransform."
                  << std::endl;
        return -1;
      }
      // Should just write out the rigid transform here.
      if ( strippedOutputTransform.size() > 0 )
      {
        typename VersorRigid3DTransformType::Pointer versorRigid =
          itk::ComputeRigidTransformFromGeneric( genericComponent.GetPointer() );
        if ( versorRigid.IsNotNull() )
        {
          itk::WriteTransformToDisk< TInputScalarType, TWriteScalarType >( versorRigid.GetPointer(),
                                                                           strippedOutputTransform );
        }
      }
    }
    catch ( itk::ExceptionObject & excp )
    {
      throw excp; // rethrow exception, handle in some other scope.
    }
  }
  return 0;
}

template int
WriteBothTransformsToDisk< double, double >( const itk::Transform< double, 3, 3 >::ConstPointer genericTransformToWrite,
                                             const std::string &                                outputTransform,
                                             const std::string & strippedOutputTransform );
template int
WriteBothTransformsToDisk< double, float >( const itk::Transform< double, 3, 3 >::ConstPointer genericTransformToWrite,
                                            const std::string &                                outputTransform,
                                            const std::string & strippedOutputTransform );
/*
// INFO: to make this function work for input float type, we need to change "ComputeRigidTransformFromGeneric" function
to be a template over scalarType
//
template int WriteBothTransformsToDisk<float,double>(const itk::Transform<float, 3, 3>::ConstPointer
genericTransformToWrite, const std::string & outputTransform, const std::string & strippedOutputTransform); template int
WriteBothTransformsToDisk<float,float>(const itk::Transform<float, 3, 3>::ConstPointer genericTransformToWrite, const
std::string & outputTransform, const std::string & strippedOutputTransform);
*/

template < typename TInputScalarType, typename TWriteScalarType >
int
WriteStrippedRigidTransformToDisk(
  const typename itk::Transform< TInputScalarType, 3, 3 >::ConstPointer genericTransformToWrite,
  const std::string &                                                   strippedOutputTransform )
{
  return WriteBothTransformsToDisk< TInputScalarType, TWriteScalarType >(
    genericTransformToWrite, std::string( "" ), strippedOutputTransform );
}

template int
WriteStrippedRigidTransformToDisk< double, double >(
  const itk::Transform< double, 3, 3 >::ConstPointer genericTransformToWrite,
  const std::string &                                strippedOutputTransform );
template int
WriteStrippedRigidTransformToDisk< double, float >(
  const itk::Transform< double, 3, 3 >::ConstPointer genericTransformToWrite,
  const std::string &                                strippedOutputTransform );

template < typename TScalarType >
typename itk::Transform< TScalarType, 3, 3 >::Pointer
ReadTransformFromDisk( const std::string & initialTransform )
{
  typename itk::Transform< TScalarType, 3, 3 >::Pointer genericTransform = nullptr;

  using ThinPlateSpline3DTransformType = typename itk::ThinPlateR2LogRSplineKernelTransform< TScalarType, 3 >;
  using ScaleVersor3DTransformType = typename itk::ScaleVersor3DTransform< TScalarType >;
  using ScaleSkewVersor3DTransformType = typename itk::ScaleSkewVersor3DTransform< TScalarType >;
  using VersorRigid3DTransformType = typename itk::VersorRigid3DTransform< TScalarType >;
  using AffineTransformType = typename itk::AffineTransform< TScalarType, 3 >;
  using BSplineTransformType = typename itk::BSplineTransform< TScalarType, 3, 3 >;
  using BRAINSCompositeTransformType = typename itk::CompositeTransform< TScalarType, 3 >;
  using TransformFileReaderType = typename itk::TransformFileReaderTemplate< TScalarType >;

  typename TransformFileReaderType::Pointer transformListReader = TransformFileReaderType::New();

  using TransformListType = typename TransformFileReaderType::TransformListType;

  TransformListType currentTransformList;

  typename TransformFileReaderType::TransformListType::const_iterator currentTransformListIterator;

  try
  {
    if ( initialTransform.size() > 0 )
    {
      std::cout << "Read ITK transform from file: " << initialTransform << std::endl;

      transformListReader->SetFileName( initialTransform.c_str() );
      transformListReader->Update();

      currentTransformList = *( transformListReader->GetTransformList() );
      if ( currentTransformList.size() == 0 )
      {
        itkGenericExceptionMacro( << "Number of currentTransformList = " << currentTransformList.size()
                                  << "FATAL ERROR: Failed to read transform" << initialTransform );
      }
      unsigned i = 0;
      for ( typename TransformListType::const_iterator it = currentTransformList.begin();
            it != currentTransformList.end();
            ++it, ++i )
      {
        std::cout << "HACK: " << i << "  " << ( *( it ) )->GetNameOfClass() << std::endl;
      }
    }
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "[FAILED]" << std::endl;
    std::cerr << "Error while reading the the file " << initialTransform << std::endl;
    std::cerr << excp << std::endl;
    throw;
  }

  if ( currentTransformList.size() == 1 ) // Most simple transform types
  {
    // NOTE:  The dynamic casting here circumvents the standard smart pointer
    // behavior.  It is important that
    // by making a new copy and transfering the parameters it is more safe.  Now
    // we only need to ensure
    // that currentTransformList.begin() smart pointer is stable (i.e. not
    // deleted) while the variable
    // temp has a reference to it's internal structure.
    const std::string transformFileType = ( *( currentTransformList.begin() ) )->GetNameOfClass();
    if ( transformFileType == "VersorRigid3DTransform" )
    {
      const typename VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< VersorRigid3DTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename VersorRigid3DTransformType::Pointer tempCopy = VersorRigid3DTransformType::New();
      AssignRigid::AssignConvertedTransform( tempCopy, tempInitializerITKTransform );
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "ScaleVersor3DTransform" )
    {
      const typename ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< ScaleVersor3DTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename ScaleVersor3DTransformType::Pointer tempCopy = ScaleVersor3DTransformType::New();
      AssignRigid::AssignConvertedTransform( tempCopy, tempInitializerITKTransform );
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "ScaleSkewVersor3DTransform" )
    {
      const typename ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< ScaleSkewVersor3DTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename ScaleSkewVersor3DTransformType::Pointer tempCopy = ScaleSkewVersor3DTransformType::New();
      AssignRigid::AssignConvertedTransform( tempCopy, tempInitializerITKTransform );
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "AffineTransform" )
    {
      const typename AffineTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< AffineTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename AffineTransformType::Pointer tempCopy = AffineTransformType::New();
      AssignRigid::AssignConvertedTransform( tempCopy, tempInitializerITKTransform );
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "ThinPlateR2LogRSplineKernelTransform" )
    {
      const typename ThinPlateSpline3DTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< ThinPlateSpline3DTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename ThinPlateSpline3DTransformType::Pointer tempCopy = ThinPlateSpline3DTransformType::New();
      tempCopy->SetFixedParameters( tempInitializerITKTransform->GetFixedParameters() );
      tempCopy->SetParametersByValue( tempInitializerITKTransform->GetParameters() );
      tempCopy->ComputeWMatrix();
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "BSplineTransform" )
    {
      const typename BSplineTransformType::ConstPointer tempInitializerITKTransform =
        static_cast< BSplineTransformType const * >( ( *( currentTransformList.begin() ) ).GetPointer() );
      typename BSplineTransformType::Pointer tempCopy = BSplineTransformType::New();
      tempCopy->SetFixedParameters( tempInitializerITKTransform->GetFixedParameters() );
      tempCopy->SetParametersByValue( tempInitializerITKTransform->GetParameters() );
      genericTransform = tempCopy.GetPointer();
    }
    else if ( transformFileType == "CompositeTransform" )
    {
      try
      {
        const typename BRAINSCompositeTransformType::ConstPointer tempInitializerITKTransform =
          static_cast< const BRAINSCompositeTransformType * >( ( *( currentTransformList.begin() ) ).GetPointer() );
        typename BRAINSCompositeTransformType::Pointer tempCopy = BRAINSCompositeTransformType::New();
        const typename BRAINSCompositeTransformType::TransformQueueType & transformQueue =
          tempInitializerITKTransform->GetTransformQueue();
        for ( unsigned int i = 0; i < transformQueue.size(); ++i )
        {
          tempCopy->AddTransform( tempInitializerITKTransform->GetNthTransform( i ) );
        }
        genericTransform = tempCopy.GetPointer();
      }
      catch ( itk::ExceptionObject & excp )
      {
        std::cerr << "[FAILED]" << std::endl;
        std::cerr << "Error while reading the the file " << initialTransform << std::endl;
        std::cerr << excp << std::endl;
        throw excp;
      }
    }
    else
    {
      std::cerr << "ERROR:  Invalid type (" << transformFileType << ") " << __FILE__ << " " << __LINE__ << std::endl;
    }
  }
  else if ( currentTransformList.size() > 1 )
  {
    std::cout << "Adding all transforms in the list to a composite transform file..." << std::endl;
    try
    {
      typename BRAINSCompositeTransformType::Pointer tempCopy = BRAINSCompositeTransformType::New();

      for ( typename TransformListType::const_iterator it = currentTransformList.begin();
            it != currentTransformList.end();
            ++it )
      {
        tempCopy->AddTransform( static_cast< itk::Transform< TScalarType, 3, 3 > * >( ( *it ).GetPointer() ) );
      }
      genericTransform = tempCopy.GetPointer();
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "[FAILED]" << std::endl;
      std::cerr << "Error while adding all input transforms to a composite transform file " << initialTransform
                << std::endl;
      std::cerr << excp << std::endl;
      throw excp;
    }
  }
  return genericTransform;
}

template itk::Transform< double, 3, 3 >::Pointer
ReadTransformFromDisk< double >( const std::string & initialTransform );
//
// INFO: to make ReadTransformFromDisk function work with float type, we need to change "AssignConvertedTransform"
// functions
//       to be templates over ScalarType
// template itk::Transform<float, 3, 3>::Pointer ReadTransformFromDisk<float>(const std::string & initialTransform);

itk::Transform< double, 3, 3 >::Pointer
ReadTransformFromDisk( const std::string & initialTransform )
{
  return ReadTransformFromDisk< double >( initialTransform );
}

template < typename TInputScalarType, typename TWriteScalarType >
void WriteTransformToDisk( itk::Transform< TInputScalarType, 3, 3 > const * const MyTransform,
                           const std::string &                                    TransformFilename )
{
  /*
   *  Convert the transform to the appropriate assumptions and write it out as requested.
   */

  std::string transformFileType = MyTransform->GetNameOfClass();
  if ( !MyTransform )
  {
    transformFileType = MyTransform->GetNameOfClass();
  }

  // First check if the input transform is a displacementField transform
  using InputDisplacementFieldTransformType = itk::DisplacementFieldTransform< TInputScalarType, 3 >;
  if ( transformFileType ==
       "DisplacementFieldTransform" ) // if it's a displacement field transform, we write that as a float displacement
  {
    const InputDisplacementFieldTransformType * dispXfrm =
      static_cast< const InputDisplacementFieldTransformType * >( MyTransform );
    using InputDisplacementFieldType = typename InputDisplacementFieldTransformType::DisplacementFieldType;
    typename InputDisplacementFieldType::ConstPointer inputDispField = dispXfrm->GetDisplacementField();

    // Define a float type displacement field
    using VectorType = itk::Vector< TWriteScalarType, 3 >;
    using OutputDisplacementFieldType = itk::Image< VectorType, 3 >;
    typename OutputDisplacementFieldType::Pointer outputDispField = OutputDisplacementFieldType::New();
    outputDispField->CopyInformation( inputDispField );
    outputDispField->SetRegions( inputDispField->GetLargestPossibleRegion() );
    outputDispField->Allocate();

    itk::ImageRegionConstIterator< InputDisplacementFieldType > inIt( inputDispField,
                                                                      inputDispField->GetLargestPossibleRegion() );
    inIt.GoToBegin();
    itk::ImageRegionIterator< OutputDisplacementFieldType > outIt( outputDispField,
                                                                   outputDispField->GetLargestPossibleRegion() );
    outIt.GoToBegin();
    while ( !outIt.IsAtEnd() )
    {
      outIt.Set( static_cast< VectorType >( inIt.Get() ) );
      ++inIt;
      ++outIt;
    }

    using DisplacementFieldWriter = itk::ImageFileWriter< OutputDisplacementFieldType >;
    typename DisplacementFieldWriter::Pointer dispWriter = DisplacementFieldWriter::New();
    dispWriter->SetInput( outputDispField );
    dispWriter->SetFileName( TransformFilename.c_str() );
    try
    {
      dispWriter->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
      std::cerr << "Can't write the displacement field transform file " << TransformFilename << std::endl;
      std::cerr << "Exception Object caught: " << std::endl;
      std::cerr << err << std::endl;
    }
  }
  else // other transforms (not displacementField transforms)
  {
    using TransformWriterType = itk::TransformFileWriterTemplate< TWriteScalarType >;
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName( TransformFilename.c_str() );
    transformWriter->SetUseCompression( true );

    const std::string extension = itksys::SystemTools::GetFilenameLastExtension( TransformFilename );
    std::string       inverseTransformFileName( TransformFilename );
    inverseTransformFileName.replace(
      inverseTransformFileName.end() - extension.size(), inverseTransformFileName.end(), "_Inverse.h5" );
    typename TransformWriterType::Pointer inverseTransformWriter = TransformWriterType::New();
    inverseTransformWriter->SetFileName( inverseTransformFileName.c_str() );
    inverseTransformWriter->SetUseCompression( true );
    bool inverseTransformExists = true;

    // if the transform to write is not displacementField transform.
    {
      transformWriter->AddTransform( MyTransform );
      if ( MyTransform->GetInverseTransform().IsNull() )
      {
        inverseTransformExists = false;
      }
      else
      {
        /*
         First, we need to check whether the inverse transform has the type of "MatrixOffsetTransformBase".
         In such cases, the inverse transform cannot be written in a different precision type (e.g. float)
         rather than its original precision (e.g. double) because the TransformFileWriter will return an error:
         ------
         itk::ERROR: Could not create an instance of MatrixOffsetTransformBase_float_3_3
         ------
         Therefore, we cast the inverse transform to the AffineTransform type before it is passed to
         TransformFileWriter.
         */
        using GenericTransformType = typename itk::MatrixOffsetTransformBase< TInputScalarType, 3, 3 >;
        // Check if MyTransform is a GenericTransformType
        if ( transformFileType.find( "AffineTransform" ) != std::string::npos ||
             transformFileType == "MatrixOffsetTransformBase" || transformFileType == "Rigid3DTransform" ||
             transformFileType == "Euler3DTransform" || transformFileType == "CenteredEuler3DTransform" ||
             transformFileType == "QuaternionRigidTransform" || transformFileType == "VersorTransform" ||
             transformFileType == "VersorRigid3DTransform" || transformFileType == "ScaleSkewVersor3DTransform" ||
             transformFileType == "ScaleVersor3DTransform" || transformFileType == "Similarity3DTransform" ||
             transformFileType == "ScaleTransform" || transformFileType == "ScaleLogarithmicTransform" )
        {
          typename GenericTransformType::Pointer genericTransform =
            static_cast< GenericTransformType * >( MyTransform->GetInverseTransform().GetPointer() );
          using InverseTransformType = itk::AffineTransform< TInputScalarType, 3 >;
          const typename InverseTransformType::Pointer tempInverseTransform = InverseTransformType::New();
          tempInverseTransform->SetMatrix( genericTransform->GetMatrix() );
          tempInverseTransform->SetOffset( genericTransform->GetOffset() );

          inverseTransformWriter->AddTransform( tempInverseTransform );
        }
        else
        {
          inverseTransformWriter->AddTransform( MyTransform->GetInverseTransform() );
        }
      }
    }
    try
    {
      transformWriter->Update();
    }
    catch ( itk::ExceptionObject & excp )
    {
      throw excp;
    }
    if ( inverseTransformExists )
    {
      try
      {
        inverseTransformWriter->Update();
      }
      catch ( itk::ExceptionObject & excp )
      {
        throw excp;
        // Writing the inverseTransform is optional,  throw excp;
      }
    }
  }
  // Test if the forward file exists.
  if ( !itksys::SystemTools::FileExists( TransformFilename.c_str() ) )
  {
    itk::ExceptionObject e( __FILE__, __LINE__, "Failed to write file", "WriteTransformToDisk" );
    std::ostringstream   msg;
    msg << "The file was not successfully created. " << std::endl << "Filename = " << TransformFilename << std::endl;
    e.SetDescription( msg.str().c_str() );
    throw e;
  }
}

template void
WriteTransformToDisk< double, double >( itk::Transform< double, 3, 3 > const * const MyTransform,
                                        const std::string &                          TransformFilename );
template void
WriteTransformToDisk< double, float >( itk::Transform< double, 3, 3 > const * const MyTransform,
                                       const std::string &                          TransformFilename );
template void
WriteTransformToDisk< float, double >( itk::Transform< float, 3, 3 > const * const MyTransform,
                                       const std::string &                         TransformFilename );
template void
WriteTransformToDisk< float, float >( itk::Transform< float, 3, 3 > const * const MyTransform,
                                      const std::string &                         TransformFilename );

template < typename TScalarType >
void WriteTransformToDisk( itk::Transform< TScalarType, 3, 3 > const * const MyTransform,
                           const std::string &                               TransformFilename )
{
  WriteTransformToDisk< TScalarType, TScalarType >( MyTransform, TransformFilename );
}

template void
WriteTransformToDisk< double >( itk::Transform< double, 3, 3 > const * const MyTransform,
                                const std::string &                          TransformFilename );
template void
WriteTransformToDisk< float >( itk::Transform< float, 3, 3 > const * const MyTransform,
                               const std::string &                         TransformFilename );

} // end namespace itk
