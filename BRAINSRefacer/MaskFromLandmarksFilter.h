//
// Created by Leinoff, Alexander on 6/17/16.
//

#ifndef BRAINSTOOLS_MASKFROMLANDMARKS_H
#define BRAINSTOOLS_MASKFROMLANDMARKS_H

#include <Slicer3LandmarkIO.h>
#include <itkMacro.h>
#include <itkImageToImageFilter.h>
#include <itkPoint.h>

template <typename TInputImage, typename TInputMask>
class MaskFromLandmarksFilter
  : public itk::ImageToImageFilter< TInputImage, TInputMask>
{
public:
  typedef MaskFromLandmarksFilter Self;
  typedef itk::SmartPointer<Self> Pointer;

 // const unsigned int Dimension = 3;

  itkNewMacro(Self);
  itkTypeMacro(Self, ImageToImageFilter)

  //this doesn't work?? why not ?? \/\/
  //itkSetMacro(Landmarks, LandmarksMapType)
  itkSetMacro(LandmarksFileName, std::string)
  itkSetMacro(ReverseMask, bool)

protected:
  MaskFromLandmarksFilter()
  {
    this->m_ReverseMask = false;
  }
  ~MaskFromLandmarksFilter() {};

  void GenerateData() ITK_OVERRIDE
  {
    typedef TInputMask ImageMaskType;
    typedef ImageMaskType OutputImageType;

    typedef typename OutputImageType::Pointer      ImagePointer;
    typedef typename TInputImage::ConstPointer ImageConstPointer;

    LandmarksMapType myLandmarks = ReadSlicer3toITKLmk(m_LandmarksFileName);
    typedef itk::Point<double, 3> PointType;

    PointType rightEye   = myLandmarks.find("RE")->second;
    PointType leftEye    = myLandmarks.find("LE")->second;
    PointType dens_axis  = myLandmarks.find("dens_axis")->second;

    std::cout << "rightEye:\t" << rightEye << std::endl;
    std::cout << "leftEye:\t" << leftEye << std::endl;
    std::cout << "dens_axis:\t" << dens_axis << std::endl;


    // find the a,b,c for the plane equation ax + by + c = 0
    // first get two vectors in the plane u and v
    typedef itk::Vector<double,3> VectorType;
    VectorType u = rightEye - leftEye;
    VectorType v = dens_axis - leftEye;

    std::cout << "vector u: \t" << u << std::endl;
    std::cout << "vector v: \t" << v << std::endl;

    VectorType cross = itk::CrossProduct(u, v);
    std::cout << "vector cross: \t" << cross <<std::endl;

    //for ax + by + cz = d plug in cross for abc and one of the points (dens_axis) for xyz
    VectorType leftEyeVector;
    leftEyeVector[0] = leftEye[0];
    leftEyeVector[1] = leftEye[1];
    leftEyeVector[2] = leftEye[2];

    VectorType::ComponentType d = cross * leftEyeVector;
    std::cout << "d:\t" << d << std::endl;

    std::cout << "equation of plane is:" << std::endl;
    std::cout << cross[0] << "x + " << cross[1] << "y + " << cross[2] << "z = " << d << std::endl;

    // make a mask based on the plane:

    // first create mask image
//    ImageMaskType::Pointer maskImageLM = ImageMaskType::New();
    ImageConstPointer inputImage = this->GetInput();
    ImagePointer      outputImage = this->GetOutput();

    //OutputImageRegionType outputRegion = inputImage->GetLargestPossibleRegion();
    outputImage->SetOrigin(inputImage->GetOrigin());
    outputImage->SetSpacing(inputImage->GetSpacing());
    outputImage->SetDirection(inputImage->GetDirection());
    outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outputImage->Allocate();

    typedef itk::ImageRegionConstIterator<TInputImage> ImageConstIterator;
    ImageConstIterator inputIterator(inputImage, inputImage->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<ImageMaskType> MaskIteratorType;
    MaskIteratorType outputIterator(outputImage, outputImage->GetLargestPossibleRegion());

    inputIterator.GoToBegin();
    outputIterator.GoToBegin();

    // go through image to see if pixel is in mask or not??
    unsigned char maskVal = this->m_ReverseMask ? 0 : 1;
    unsigned char unMaskVal = this->m_ReverseMask ? 1 : 0;
    while(!inputIterator.IsAtEnd())
      {
      PointType currentPoint;
      inputImage->TransformIndexToPhysicalPoint(inputIterator.GetIndex(), currentPoint);
      //  std::cout << currentPoint <<std::endl;

      VectorType currentVector;
      currentVector[0] = currentPoint[0];
      currentVector[1] = currentPoint[1];
      currentVector[2] = currentPoint[2];

      if( d > cross * currentVector )
        {
        outputIterator.Set( maskVal ); // normally 1
        }
      else
        {
        outputIterator.Set( unMaskVal ); //normally 0
        }
      ++outputIterator;
      ++inputIterator;
      }
  }

private:
  std::string      m_LandmarksFileName;
  bool             m_ReverseMask;
};


#endif //BRAINSTOOLS_MASKFROMLANDMARKS_H
