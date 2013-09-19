#include <iostream>
#include "GenericTransformImage.h"

#include "itkAffineTransform.h"
#include "itkBSplineDeformableTransform.h"
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
#include "itkTransformFactory.h"
#include <itksys/SystemTools.hxx>

// #include "itkSimilarity2DTransfor3DPerspectiveTransform.h"

namespace itk
{
VersorRigid3DTransformType::Pointer ComputeRigidTransformFromGeneric(
  const GenericTransformType::ConstPointer genericTransformToWrite)
{
  typedef VersorRigid3DTransformType VersorRigidTransformType;
  VersorRigidTransformType::Pointer versorRigid = VersorRigidTransformType::New();
  versorRigid->SetIdentity();
  // //////////////////////////////////////////////////////////////////////////
  // ConvertTransforms
  if( genericTransformToWrite.IsNotNull() )
    {
    try
      {
      const std::string transformFileType = genericTransformToWrite->GetNameOfClass();
      if( transformFileType == "VersorRigid3DTransform" )
        {
        const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
          dynamic_cast<VersorRigid3DTransformType const *>( genericTransformToWrite.GetPointer() );
        if( tempInitializerITKTransform.IsNull() )
          {
          itkGenericExceptionMacro(<< "Error in type conversion");
          }
        AssignRigid::ExtractVersorRigid3DTransform(versorRigid, tempInitializerITKTransform);
        }
      else if( transformFileType == "ScaleVersor3DTransform" )
        {
        const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
          dynamic_cast<ScaleVersor3DTransformType const *>( genericTransformToWrite.GetPointer() );
        if( tempInitializerITKTransform.IsNull() )
          {
          itkGenericExceptionMacro(<< "Error in type conversion");
          }
        AssignRigid::ExtractVersorRigid3DTransform(versorRigid, tempInitializerITKTransform);
        }
      else if( transformFileType == "ScaleSkewVersor3DTransform" )
        {
        const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
          dynamic_cast<ScaleSkewVersor3DTransformType const *>( genericTransformToWrite.GetPointer() );
        if( tempInitializerITKTransform.IsNull() )
          {
          itkGenericExceptionMacro(<< "Error in type conversion");
          }
        AssignRigid::ExtractVersorRigid3DTransform(versorRigid, tempInitializerITKTransform);
        }
      else if( transformFileType == "AffineTransform" )
        {
        const AffineTransformType::ConstPointer tempInitializerITKTransform =
          dynamic_cast<AffineTransformType const *>( genericTransformToWrite.GetPointer() );
        if( tempInitializerITKTransform.IsNull() )
          {
          itkGenericExceptionMacro(<< "Error in type conversion");
          }
        AssignRigid::ExtractVersorRigid3DTransform(versorRigid, tempInitializerITKTransform);
        }
      /*
        * else if(transformFileType == "BSplineDeformableTransform")
        * {
        * BSplineTransformType::Pointer tempInitializerITKTransform
        * = dynamic_cast<BSplineTransformType *>(
        *    genericTransformToWrite.GetPointer() );
                if ( tempInitializerITKTransform.IsNull() )
                  {
itkGenericExceptionMacro(<< "Error in type conversion");
                  }
        * //AssignRigid::ExtractVersorRigid3DTransform(versorRigid,
        *    tempInitializerITKTransform->GetBulkTransform());
        * versorRigid=NULL; //NOT: Perhaps it makes sense to extract the rigid
        *    part of the bulk transform.  But that is pretty obscure case.
        * return NULL;
        * }
        */
      else      //  NO SUCH CASE!!
        {
        std::cout
          << "Unsupported initial transform file -- TransformBase first transform typestring, "
          << transformFileType
          << " not equal to any recognized type VersorRigid3DTransform OR "
          << " ScaleVersor3DTransform OR ScaleSkewVersor3DTransform OR AffineTransform"
          << std::endl;
        return NULL;
        }
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cout << "[FAILED]" << std::endl;
      std::cerr
        << "Error while converting the genericTransformToWrite to Rigid"
        << std::endl;
      std::cerr << excp << std::endl;
      return NULL;
      }
    }
  return versorRigid;
}

int WriteBothTransformsToDisk(const GenericTransformType::ConstPointer genericTransformToWrite,
                              const std::string & outputTransform,
                              const std::string & strippedOutputTransform)
{
  // //////////////////////////////////////////////////////////////////////////
  // Write out tranfoms.
  if( genericTransformToWrite.IsNull() )
    {
    return 0;
    }
  try
    {
    const std::string transformFileType = genericTransformToWrite->GetNameOfClass();
    if( transformFileType == "VersorRigid3DTransform" )
      {
      if( outputTransform.size() > 0 )  // Write out the transform
        {
        itk::WriteTransformToDisk<double>(genericTransformToWrite, outputTransform);
        }
      }
    else if( transformFileType == "ScaleVersor3DTransform" )
      {
      if( outputTransform.size() > 0 )  // Write out the transform
        {
        itk::WriteTransformToDisk<double>(genericTransformToWrite, outputTransform);
        }
      }
    else if( transformFileType == "ScaleSkewVersor3DTransform" )
      {
      if( outputTransform.size() > 0 )  // Write out the transform
        {
        itk::WriteTransformToDisk<double>(genericTransformToWrite, outputTransform);
        }
      }
    else if( transformFileType == "AffineTransform" )
      {
      if( outputTransform.size() > 0 )  // Write out the transform
        {
        itk::WriteTransformToDisk<double>(genericTransformToWrite, outputTransform);
        }
      }
    else if( transformFileType == "BSplineDeformableTransform" )
      {
      const BSplineTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<BSplineTransformType  const *>( genericTransformToWrite.GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      if( strippedOutputTransform.size() > 0 )
        {
        std::cout << "ERROR:  The rigid component of a BSpline transform is not supported." << std::endl;
        }
      if( outputTransform.size() > 0 )
        {
        itk::WriteTransformToDisk<double>(genericTransformToWrite, outputTransform);
        }
      }
    else      //  NO SUCH CASE!!
      {
      std::cout << "Unsupported initial transform file -- TransformBase first transform typestring, "
                << transformFileType
                << " not equal to any recognized type VersorRigid3DTransform OR"
                << " ScaleVersor3DTransform OR ScaleSkewVersor3DTransform OR AffineTransform"
                << std::endl;
      return -1;
      }
    // Should just write out the rigid transform here.
    if( strippedOutputTransform.size() > 0  )
      {
      typedef VersorRigid3DTransformType VersorRigidTransformType;
      VersorRigidTransformType::Pointer versorRigid = itk::ComputeRigidTransformFromGeneric(genericTransformToWrite);
      if( versorRigid.IsNotNull() )
        {
        itk::WriteTransformToDisk<double>(versorRigid.GetPointer(), strippedOutputTransform);
        }
      }
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cout << "Exception: " << excp << std::endl;
    throw excp; // reohrow exception, handle in some other scope.
    }
  return 0;
}

int WriteStrippedRigidTransformToDisk(const GenericTransformType::ConstPointer genericTransformToWrite,
                                      const std::string & strippedOutputTransform)
{
  return WriteBothTransformsToDisk(genericTransformToWrite, std::string(""), strippedOutputTransform);
}

GenericTransformType::Pointer ReadTransformFromDisk(const std::string & initialTransform)
{
  GenericTransformType::Pointer genericTransform = NULL;
  // read in the initial ITKTransform
  TransformReaderType::Pointer      transformListReader =  TransformReaderType::New();
  TransformListType                 currentTransformList;
  TransformListType::const_iterator currentTransformListIterator;

  try
    {
    if( initialTransform.size() > 0 )
      {
      std::cout << "Read ITK transform from file: " << initialTransform << std::endl;

      transformListReader->SetFileName( initialTransform.c_str() );
      transformListReader->Update();

      currentTransformList = *( transformListReader->GetTransformList() );
      if( currentTransformList.size() == 0 )
        {
        itkGenericExceptionMacro( << "Number of currentTransformList = " << currentTransformList.size()
                                  << "FATAL ERROR: Failed to read transform" << initialTransform);
        }
      std::cout << "HACK: " << currentTransformList.size() << "  "
                << ( *( currentTransformList.begin() ) )->GetNameOfClass() << std::endl;
      }
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "[FAILED]" << std::endl;
    std::cerr << "Error while reading the the file " << initialTransform << std::endl;
    std::cerr << excp << std::endl;
    throw;
    }

  if( currentTransformList.size() == 1 )  // Most simple transform types
    {
    // NOTE:  The dynamic casting here circumvents the standard smart pointer
    // behavior.  It is important that
    // by making a new copy and transfering the parameters it is more safe.  Now
    // we only need to ensure
    // that currentTransformList.begin() smart pointer is stable (i.e. not
    // deleted) while the variable
    // temp has a reference to it's internal structure.
    const std::string transformFileType = ( *( currentTransformList.begin() ) )->GetNameOfClass();
    if( transformFileType == "VersorRigid3DTransform" )
      {
      const VersorRigid3DTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<VersorRigid3DTransformType const *>( ( *( currentTransformList.begin() ) ).GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      VersorRigid3DTransformType::Pointer tempCopy = VersorRigid3DTransformType::New();
      AssignRigid::AssignConvertedTransform(tempCopy,
                                            tempInitializerITKTransform);
      genericTransform = tempCopy.GetPointer();
      }
    else if( transformFileType == "ScaleVersor3DTransform" )
      {
      const ScaleVersor3DTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<ScaleVersor3DTransformType const *>( ( *( currentTransformList.begin() ) ).GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      ScaleVersor3DTransformType::Pointer tempCopy = ScaleVersor3DTransformType::New();
      AssignRigid::AssignConvertedTransform(tempCopy,
                                            tempInitializerITKTransform);
      genericTransform = tempCopy.GetPointer();
      }
    else if( transformFileType == "ScaleSkewVersor3DTransform" )
      {
      const ScaleSkewVersor3DTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<ScaleSkewVersor3DTransformType const *>( ( *( currentTransformList.begin() ) ).GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      ScaleSkewVersor3DTransformType::Pointer tempCopy = ScaleSkewVersor3DTransformType::New();
      AssignRigid::AssignConvertedTransform(tempCopy,
                                            tempInitializerITKTransform);
      genericTransform = tempCopy.GetPointer();
      }
    else if( transformFileType == "AffineTransform" )
      {
      const AffineTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<AffineTransformType const *>( ( *( currentTransformList.begin() ) ).GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      AffineTransformType::Pointer tempCopy = AffineTransformType::New();
      AssignRigid::AssignConvertedTransform(tempCopy,
                                            tempInitializerITKTransform);
      genericTransform = tempCopy.GetPointer();
      }
    else if( transformFileType == "ThinPlateR2LogRSplineKernelTransform" )
      {
      const ThinPlateSpline3DTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<ThinPlateSpline3DTransformType const *>( ( *( currentTransformList.begin() ) ).GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      ThinPlateSpline3DTransformType::Pointer tempCopy = ThinPlateSpline3DTransformType::New();
      tempCopy->SetFixedParameters( tempInitializerITKTransform->GetFixedParameters() );
      tempCopy->SetParametersByValue( tempInitializerITKTransform->GetParameters() );
      tempCopy->ComputeWMatrix();
      genericTransform = tempCopy.GetPointer();
      }
#if (ITK_VERSION_MAJOR > 3)
    else if( transformFileType == "CompositeTransform" )
      {
      try
        {
        const CompositeTransformType::ConstPointer tempInitializerITKTransform =
          dynamic_cast<const CompositeTransformType *>( ( *( currentTransformList.begin() ) ).GetPointer() );
        if( tempInitializerITKTransform.IsNull() )
          {
          itkGenericExceptionMacro(<< "Error in type conversion");
          }
        CompositeTransformType::Pointer                    tempCopy = CompositeTransformType::New();
        const CompositeTransformType::TransformQueueType & transformQueue =
          tempInitializerITKTransform->GetTransformQueue();
        for( unsigned int i = 0; i < transformQueue.size(); ++i )
          {
          tempCopy->AddTransform(tempInitializerITKTransform->GetNthTransform(i) );
          }
        genericTransform = tempCopy.GetPointer();
        }
      catch( itk::ExceptionObject & excp )
        {
        std::cerr << "[FAILED]" << std::endl;
        std::cerr << "Error while reading the the file " << initialTransform << std::endl;
        std::cerr << excp << std::endl;
        throw excp;
        }
      }
#endif
    else
      {
      std::cerr << "ERROR:  Invalid type (" << transformFileType << ") " << __FILE__ << " " << __LINE__ << std::endl;
      }
    }
  else if( currentTransformList.size() == 2 )  // A special case for
                                               // BSplineTransforms
  // To recombine the bulk and the bSpline transforms.
    {
    // transformListReader->GetTransformList();
    TransformListType::const_iterator initializeTransformsListIterator =
      currentTransformList.begin();

    const GenericTransformType::ConstPointer FirstTransform = dynamic_cast<GenericTransformType const *>
      ( ( *( initializeTransformsListIterator ) ).GetPointer() );
    if( FirstTransform.IsNull() )
      {
      itkGenericExceptionMacro(<< "Error in type conversion");
      }
    const std::string FirstTransformFileType = FirstTransform->GetNameOfClass();

    ++initializeTransformsListIterator; // Increment to next iterator

    const GenericTransformType::ConstPointer SecondTransform = dynamic_cast<GenericTransformType const *>
      ( ( *( initializeTransformsListIterator ) ).GetPointer() );
    if( SecondTransform.IsNull() )
      {
      itkGenericExceptionMacro(<< "Error in type conversion");
      }
    const std::string SecondTransformFileType = SecondTransform->GetNameOfClass();

    BSplineTransformType::Pointer outputBSplineTransform =
      BSplineTransformType::New();
    outputBSplineTransform->SetIdentity();

    // Now get the BSpline information
    if( FirstTransformFileType == "BSplineDeformableTransform" )
      {
      const BSplineTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<BSplineTransformType const *>
        ( FirstTransform.GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      outputBSplineTransform->SetFixedParameters( tempInitializerITKTransform->GetFixedParameters() );
      outputBSplineTransform->SetParametersByValue( tempInitializerITKTransform->GetParameters() );
      outputBSplineTransform->SetBulkTransform(SecondTransform);
      }
    else if( SecondTransformFileType == "BSplineDeformableTransform" )
      {
      const BSplineTransformType::ConstPointer tempInitializerITKTransform =
        dynamic_cast<BSplineTransformType const *>
        ( SecondTransform.GetPointer() );
      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      outputBSplineTransform->SetFixedParameters( tempInitializerITKTransform->GetFixedParameters() );
      outputBSplineTransform->SetParametersByValue( tempInitializerITKTransform->GetParameters() );
      outputBSplineTransform->SetBulkTransform(FirstTransform);
      }
    else
      {
      std::cout << "[FAILED]" << std::endl;
      std::cerr << "Error using the currentTransformList has two elements, but"
                << " neither of them are a BSplineDeformableTransform/"
                << std::endl
                << "There should not be more than two transforms in the transform list."
                << std::endl;
      return NULL;
      }
    genericTransform = outputBSplineTransform.GetPointer();
    }
  else if( currentTransformList.size() > 2 )
    {
    // Error, too many transforms on transform list.
    std::cout << "[FAILED]" << std::endl;
    std::cerr << "Error using the currentTransformList for initializing a"
              << " BSPlineDeformableTransform:"
              << std::endl
              << "There should not be more than two transforms in the transform list."
              << std::endl;
    return NULL;
    }
  return genericTransform;
}

template<class TScalarType>
void WriteTransformToDisk( itk::Transform<TScalarType, 3, 3> const *const MyTransform, const std::string & TransformFilename )
{
  /*
    *  Convert the transform to the appropriate assumptions and write it out as
    *requested.
    */
    {
    typedef itk::TransformFileWriterTemplate<TScalarType> TransformWriterType;
    typename TransformWriterType::Pointer transformWriter =  TransformWriterType::New();
    transformWriter->SetFileName( TransformFilename.c_str() );

    const std::string extension = itksys::SystemTools::GetFilenameLastExtension(TransformFilename);
    std::string       inverseTransformFileName(TransformFilename);
    inverseTransformFileName.replace(inverseTransformFileName.end() - extension.size(),
                                     inverseTransformFileName.end(), "_Inverse.h5");
    typename TransformWriterType::Pointer inverseTransformWriter =  TransformWriterType::New();
    inverseTransformWriter->SetFileName( inverseTransformFileName.c_str() );
    const std::string transformFileType = MyTransform->GetNameOfClass();
    bool              inverseTransformExists = true;
    if( transformFileType == "BSplineDeformableTransform" )
      {
      typedef itk::BSplineDeformableTransform< TScalarType,
                                               GenericTransformImageNS::SpaceDimension,
                                               GenericTransformImageNS::SplineOrder> LocalBSplineTransformType;

      const typename LocalBSplineTransformType::ConstPointer tempInitializerITKTransform
        = dynamic_cast<LocalBSplineTransformType const *>( MyTransform );

      if( tempInitializerITKTransform.IsNull() )
        {
        itkGenericExceptionMacro(<< "Error in type conversion");
        }
      // NOTE: Order was reversed to get BSpline first, then Bulk
      // transform in an attempt to appease Slicer3.
      // Bulk transform is assumed to be second in Slicer3.
      transformWriter->AddTransform(tempInitializerITKTransform);
      transformWriter->AddTransform( tempInitializerITKTransform->GetBulkTransform() );

      if( tempInitializerITKTransform->GetInverseTransform().IsNull() )
        {
        inverseTransformExists = false;
        }
      else
        {
        inverseTransformWriter->AddTransform(tempInitializerITKTransform->GetInverseTransform() );
        }
      if( tempInitializerITKTransform->GetBulkTransform()->GetInverseTransform().IsNull() )
        {
        inverseTransformExists = false;
        }
      else
        {
        inverseTransformWriter->AddTransform( tempInitializerITKTransform->GetBulkTransform()->GetInverseTransform() );
        }
      }
    else
      {
      transformWriter->AddTransform(MyTransform);
      if( MyTransform->GetInverseTransform().IsNull() )
        {
        inverseTransformExists = false;
        }
      else
        {
        inverseTransformWriter->AddTransform(MyTransform->GetInverseTransform() );
        }
      }

    try
      {
      transformWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      throw excp;
      }
    if( inverseTransformExists )
      {
      try
        {
        inverseTransformWriter->Update();
        }
      catch( itk::ExceptionObject & excp )
        {
        throw excp;
        // Writing the inverseTransform is optional,  throw excp;
        }
      }
    // Test if the forward file exists.
    if( !itksys::SystemTools::FileExists( TransformFilename.c_str() ) )
      {
      itk::ExceptionObject e(__FILE__, __LINE__, "Failed to write file", "WriteTransformToDisk");
      std::ostringstream   msg;
      msg << "The file was not successfully created. "
          << std::endl << "Filename = " << TransformFilename
          << std::endl;
      e.SetDescription( msg.str().c_str() );
      throw e;
      }
    }
}

template void WriteTransformToDisk<double>( itk::Transform<double, 3, 3> const *const MyTransform, const std::string & TransformFilename );
template void WriteTransformToDisk<float>( itk::Transform<float, 3, 3> const *const MyTransform, const std::string & TransformFilename );

} // end namespace itk
