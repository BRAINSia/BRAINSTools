#if !defined( DebugImageViewerClient_h )
#define DebugImageViewerClient_h

#include "BRAINSCommonLib.h"

#ifdef USE_DebugImageViewer
#include <string>
#include <itksys/SystemTools.hxx>
#include <vtkClientSocket.h>
#include <itkImage.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkSpatialOrientation.h>
#include <itkSpatialOrientationAdapter.h>
#include <itkOrientImageFilter.h>
#include <stdio.h>
// #include <itkIO.h>
// #include <itkIO2.h>

namespace DebugImageViewerUtil
{
typedef itk::SpatialOrientationAdapter SOAdapterType;
typedef SOAdapterType::DirectionType   DirectionType;

template <class InputImageType, class OutputImageType>
typename OutputImageType::Pointer
ScaleAndCast(const typename InputImageType::Pointer & image,
             const typename OutputImageType::PixelType OutputMin,
             const typename OutputImageType::PixelType OutputMax)
{
  typedef itk::RescaleIntensityImageFilter<InputImageType,
                                           OutputImageType> R2CRescaleFilterType;
  typename R2CRescaleFilterType::Pointer RealToProbMapCast
    = R2CRescaleFilterType::New();
  RealToProbMapCast->SetOutputMinimum(OutputMin);
  RealToProbMapCast->SetOutputMaximum(OutputMax);
  RealToProbMapCast->SetInput(image);
  try
    {
    RealToProbMapCast->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception in Image cast." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(-1);
    }
  typename OutputImageType::Pointer returnScaled = RealToProbMapCast->GetOutput();
  return returnScaled;
}

/** AllocateImageFromExample creates and allocates an image of the type OutputImageType,
 * using TemplateImageType as the source of size and spacing...
 *
 */
template <class TemplateImageType, class OutputImageType>
typename OutputImageType::Pointer
AllocateImageFromExample(
  const typename TemplateImageType::Pointer & TemplateImage)
{
  typename OutputImageType::Pointer rval = OutputImageType::New();
  rval->CopyInformation(TemplateImage);
  rval->SetRegions( TemplateImage->GetLargestPossibleRegion() );
  rval->Allocate();
  return rval;
}

template <class ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::Pointer & inputImage,
            itk::SpatialOrientation::ValidCoordinateOrientationFlags orient)
{
  typename itk::OrientImageFilter<ImageType, ImageType>::Pointer orienter
    = itk::OrientImageFilter<ImageType, ImageType>::New();

  orienter->SetDesiredCoordinateOrientation(orient);
  orienter->UseImageDirectionOn();
  orienter->SetInput(inputImage);
  orienter->Update();
  typename ImageType::Pointer returnval
    = orienter->GetOutput();
  returnval->DisconnectPipeline();
  orienter->ReleaseDataFlagOn();
  return returnval;
}

template <class ImageType>
typename ImageType::Pointer
OrientImage(typename ImageType::Pointer & inputImage,
            const typename ImageType::DirectionType & dirCosines)
{
  return OrientImage<ImageType>
           ( inputImage,
           SOAdapterType().FromDirectionCosines(
             dirCosines) );
}
}

class DebugImageViewerClient
{
public:
  DebugImageViewerClient() : m_Sock(0), m_Enabled(false), m_PromptUser(false)
  {
  }

  ~DebugImageViewerClient()
  {
    if( m_Sock )
      {
      this->m_Sock->CloseSocket();
      this->m_Sock->Delete();
      }
  }

  void SetPromptUser(bool x)
  {
    m_PromptUser = x;
  }

  /** Send an image to the viewer */
  template <class ImageType>
  void SendImage(const typename ImageType::Pointer & image,
                 unsigned viewIndex = 0)
  {
    this->Send<ImageType>(image, viewIndex);
    if( this->m_PromptUser )
      {
      //
      // make sure we connect to interactive input
      FILE *in = fopen("/dev/tty", "r");
      std::cerr << ">>>>>>>>>Hit enter to continue " << std::flush;
      char buf[256];
      fgets(buf, 255, in);
      fclose(in);
      }
  }

  /** Send one component of a vector image to the viewer */
  template <class ImageType>
  void SendImage(const typename ImageType::Pointer & image,
                 unsigned viewIndex,
                 unsigned vectorIndex)
  {
    this->Send<ImageType>(image, viewIndex, vectorIndex);
    if( this->m_PromptUser )
      {
      //
      // make sure we connect to interactive input
      FILE *in = fopen("/dev/tty", "r");
      std::cerr << ">>>>>>>>>Hit enter to continue " << std::flush;
      char buf[256];
      fgets(buf, 255, in);
      fclose(in);
      }
  }

  /** enable sending of images to the viewer */
  void SetEnabled(bool enabled)
  {
    this->m_Enabled = enabled;
    if( enabled )
      {
      this->_Init();
      }
    else if( this->m_Sock != 0 )
      {
      this->m_Sock->CloseSocket();
      this->m_Sock->Delete();
      this->m_Sock = 0;
      }
  }

  bool Enabled()
  {
    return m_Enabled;
  }

private:
  void _Init()
  {
    this->m_Sock = vtkClientSocket::New();
    this->m_Sock->ConnectToServer("localhost", 19345);
  }

  template <class ImageType>
  void Send(const typename ImageType::Pointer & image, unsigned int viewIndex);

  template <class ImageType>
  void Send(const typename ImageType::Pointer & image, unsigned int viewIndex, unsigned int vectorIndex);

private:
  vtkClientSocket *m_Sock;
  bool             m_Enabled;
  bool             m_PromptUser;
};

template <class ImageType>
void
DebugImageViewerClient::Send(const typename ImageType::Pointer & image, unsigned int viewIndex)
{
  if( !this->Enabled() )
    {
    return;
    }
  typedef itk::Image<unsigned char, 3>   TransferImageType;
  typedef TransferImageType::SizeType    SizeType;
  typedef TransferImageType::SpacingType SpacingType;
  typedef TransferImageType::PointType   PointType;
  //
  // make sure image is in a known image type
  TransferImageType::Pointer xferImage
    = DebugImageViewerUtil::ScaleAndCast<ImageType, TransferImageType>(image, 0, 255);
  typename TransferImageType::DirectionType DesiredDirectionCos;
  DesiredDirectionCos[0][0] = 1; DesiredDirectionCos[0][1] = 0;
  DesiredDirectionCos[0][2] = 0;
  DesiredDirectionCos[1][0] = 0; DesiredDirectionCos[1][1] = 1;
  DesiredDirectionCos[1][2] = 0;
  DesiredDirectionCos[2][0] = 0; DesiredDirectionCos[2][1] = 0;
  DesiredDirectionCos[2][2] = 1;
  xferImage = DebugImageViewerUtil::OrientImage<TransferImageType>(xferImage,
                                                                   DesiredDirectionCos);
  //
  // get size
  SizeType size = xferImage->GetLargestPossibleRegion().GetSize();

  unsigned int bufferSize = size[0] * size[1] * size[2]
    * sizeof( typename TransferImageType::PixelType );

  // get spacing
  SpacingType spacing = xferImage->GetSpacing();
  // get orientation
  itk::SpatialOrientation::ValidCoordinateOrientationFlags orientation
    = itk::SpatialOrientationAdapter().FromDirectionCosines
        ( xferImage->GetDirection() );
  // get origin
  PointType origin = xferImage->GetOrigin();
  for( unsigned int i = 0; i < 3; i++ )
    {
    this->m_Sock->Send( &size[i], sizeof( SizeType::SizeValueType ) );
    }
  for( unsigned int i = 0; i < 3; i++ )
    {
    this->m_Sock->Send( &spacing[i], sizeof( SpacingType::ValueType ) );
    }
  this->m_Sock->Send( &orientation, sizeof( orientation ) );
  // send origin
  for( unsigned int i = 0; i < 3; i++ )
    {
    double x = origin[i];
    this->m_Sock->Send( &x, sizeof( double ) );
    }
  this->m_Sock->Send( &viewIndex, sizeof( viewIndex ) );
  // transfer image pixels
  this->m_Sock->Send(xferImage->GetBufferPointer(),
                     bufferSize);
  //   std::cerr << "DebugImageViewer: size = " << size
  //             << " spacing = " << spacing << std::endl
  //             << "orientation = " << orientation
  //             << " viewIndex = " <<  " bufferSize " << bufferSize
  //             << std::endl;
  //   std::cerr.flush();
}

template <class ImageType>
void
DebugImageViewerClient::Send(const typename ImageType::Pointer & image,
                             unsigned int viewIndex,
                             unsigned int vectorIndex)
{
  if( !this->Enabled() )
    {
    return;
    }
  //
  // assume ImageType is a vector type.
  typedef typename ImageType::PixelType::ComponentType ScalarPixelType;
  typedef typename itk::Image<ScalarPixelType, ImageType::ImageDimension>
    ScalarImageType;
  typename ScalarImageType::Pointer scalarImage
    = DebugImageViewerUtil::AllocateImageFromExample<ImageType, ScalarImageType>
        (image);
  typename itk::ImageRegionConstIterator<ImageType>
  sourceIt( image, image->GetLargestPossibleRegion() );
  typename itk::ImageRegionIterator<ScalarImageType>
  destIt( scalarImage, scalarImage->GetLargestPossibleRegion() );
  for( ; !sourceIt.IsAtEnd(); ++sourceIt, ++destIt )
    {
    destIt.Set(sourceIt.Get()[vectorIndex]);
    }
  this->Send<ScalarImageType>(scalarImage, viewIndex);
}

#endif // USE_DebugImageViewer
#endif // DebugImageViewerClient_h
