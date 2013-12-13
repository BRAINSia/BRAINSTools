#ifndef BRAINSCutGenerateProbability_h
#define BRAINSCutGenerateProbability_h

#include "BRAINSCutDataHandler.h"
#include "BRAINSCutConfiguration.h"
#include "itkDisplacementFieldTransform.h"
#include <itkIO.h>

class BRAINSCutGenerateProbability
{
public:
  BRAINSCutGenerateProbability( BRAINSCutDataHandler& dataHandler);

  void SetTrainingDataSetsList();

  void GenerateProbabilityMaps();

private:
  BRAINSCutDataHandler* myDataHandler;

  /** DataSets */
  std::list<DataSet *> trainingDataSetList;

  void GenerateSymmetricalSphericalCoordinateImage();

  void CreateNewFloatImageFromTemplate( WorkingImageType::Pointer & PointerToOutputImage,
                                        const WorkingImageType::Pointer & PreInitializedImage);

  void XYZToSpherical( const itk::Point<float, 3> & LocationWithOriginAtCenterOfImage, float & rho, float & phi,
                       float & theta);

  template <class WarperImageType>
  typename WarperImageType::Pointer ImageWarper(  const std::string & RegistrationFilename,
                                                  const std::string & ImageName,
                                                  typename WarperImageType::Pointer ReferenceImage  )
  {
    const bool useTransform = ( RegistrationFilename.find(".mat") != std::string::npos ||
                                RegistrationFilename.find(".h5") != std::string::npos ||
                                RegistrationFilename.find(".hdf5") != std::string::npos ||
                                RegistrationFilename.find(".txt") != std::string::npos
                                );

    typename WarperImageType::Pointer PrincipalOperandImage; // One name for the
                                                             // image to be
                                                             // warped.
      {
      typedef typename itk::ImageFileReader<WarperImageType> ReaderType;
      typename ReaderType::Pointer imageReader = ReaderType::New();

      imageReader->SetFileName(ImageName);
      imageReader->Update();
      PrincipalOperandImage = imageReader->GetOutput();
      }

    typedef double                                        VectorComponentType;
    typedef typename itk::Vector<VectorComponentType, 3> VectorPixelType;
    typedef typename itk::Image<VectorPixelType,  3>     LocalDisplacementFieldType;

    // An empty SmartPointer constructor sets up someImage.IsNull() to represent
    // a
    // not-supplied state:
    typename LocalDisplacementFieldType::Pointer DisplacementField;
    // An empty SmartPointer constructor sets up someTransform.IsNull() to
    // represent a not-supplied state:
    typename GenericTransformType::Pointer genericTransform;
    // typename WarperImageType::Pointer ReferenceImage;
    // if there is no *mat file.
    //
    // HACK ALERT HACK ALERT
    // The GenericTransformImage function would apply either a
    // 'generic' transform or a DisplacementField transform. If both
    // were non-null, the DisplacementField took precedence. To get
    // the same behavior with the DisplacementField parameter removed
    // you have to enforce this 'OR NOT AND' behavior at the point of
    // call
    if( !useTransform )  // that is, it's a warp by deformation field:
      {
      typedef typename itk::ImageFileReader<LocalDisplacementFieldType> DefFieldReaderType;
      typename DefFieldReaderType::Pointer fieldImageReader = DefFieldReaderType::New();
      fieldImageReader->SetFileName(RegistrationFilename);
      fieldImageReader->Update();
      DisplacementField = fieldImageReader->GetOutput();

      // Resample deformation filed to reference image size Check that
      // deformation
      // field and reference image have same dimensions.

      // and---  ReferenceImage.IsNull() represents the delayed default
      typedef itk::DisplacementFieldTransform<DeformationScalarType,LocalDisplacementFieldType::ImageDimension>
        DisplacementFieldTransformType;
      typename DisplacementFieldTransformType::Pointer dispXfrm =
        DisplacementFieldTransformType::New();
      dispXfrm->SetDisplacementField(DisplacementField.GetPointer());
      genericTransform = dispXfrm.GetPointer();
      }
    else // there EXIST *mat file.
      {
      std::cout << "!!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!!" << std::endl
                << "* Mat file exists!" << std::endl
                << "!!!!!!!!!!!! CAUTION !!!!!!!!!!!!!!!!!!!" << std::endl;
      genericTransform = itk::ReadTransformFromDisk(RegistrationFilename);
      }
    const double defaultValue = 0;
    const typename std::string interpolationMode = "Linear";
    const typename std::string pixelType = "short";

    typename WarperImageType::Pointer TransformedImage =
      GenericTransformImage<WarperImageType, WarperImageType, LocalDisplacementFieldType>(
        PrincipalOperandImage,
        ReferenceImage,
        genericTransform,
        defaultValue,
        interpolationMode,
        pixelType == "binary");

    return TransformedImage;
  }
};
#endif
