//
// Created by Jeffrey Obadal, Alexander Leinoff on 5/31/16.
//
#ifndef BRAINSTOOLS_CREATERANDOMBSPLINE_H
#define BRAINSTOOLS_CREATERANDOMBSPLINE_H

#include <itkBSplineTransform.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <ctime>

template<typename TInputImage,
  typename TPixelType,
    unsigned int NDimension,
      unsigned int NBSplineOrder>
class CreateRandomBSpline
  : public itk::ImageToImageFilter< TInputImage, TInputImage>
{
public:
  typedef CreateRandomBSpline Self;
  typedef itk::SmartPointer <Self> Pointer;

  typedef itk::BSplineTransform<TPixelType,NDimension,NBSplineOrder> BSplineType;
  typedef typename BSplineType::Pointer BSplinePointer;

  itkNewMacro(Self);
  itkTypeMacro(Self, ImageToImageFilter);

  itkGetMacro(BSplineOutput, BSplinePointer)

  itkSetMacro(BSplineControlPoints, unsigned int)
  itkGetMacro(BSplineControlPoints, unsigned int)

  itkSetMacro(RandMin, int)
  itkGetMacro(RandMin, int)

  itkSetMacro(RandMax, int)
  itkGetMacro(RandMax, int)

  itkSetMacro(RandScale, double)
  itkGetMacro(RandScale, double)

  itkSetMacro(Verbose, bool)
  itkSetMacro(Debug, bool)

  typedef TInputImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::PointType ImagePointType;

protected:
  CreateRandomBSpline()
  {
    m_BSplineOutput=BSplineType::New();
    m_BSplineControlPoints = 8; //default value
    m_RandMax = 5; //default value
    m_RandMin = -5; //default value
    m_RandScale = 5;

    this->m_Verbose = false;
    this->m_Debug = false;
  };
  ~CreateRandomBSpline(){};

  void GenerateData() ITK_OVERRIDE
  {
    if( m_Debug )
      {
      std::cout << "File:  " << __FILE__ << std::endl;
      std::cout << "Line:  " << __LINE__ << std::endl;
      std::cout << "In function GenerateData()" << std::endl;
      }

    typedef typename ImageType::RegionType ImageRegionType;

    const unsigned int sizeVal = 400; // The image should be large enough to contain most heads, 40 cm^3 should do it
    if(sizeVal % 2 != 0)
      {
      std::cerr << "File: " << __FILE__ << std::endl;
      std::cerr << "Line: " << __LINE__ << std::endl;
      std::cerr << "BSpline size must be divisible by 0" << std::endl;
      }
    double originVal = sizeVal * -0.5;  //automatically determine origin based on size

    typename ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename ImageType::SizeType size;
    size[0] = sizeVal;
    size[1] = sizeVal;
    size[2] = sizeVal;

    ImageRegionType subjectRegion(start, size);

    typename BSplineType::OriginType  fixedOrigin;
    fixedOrigin[0] = originVal;
    fixedOrigin[1] = originVal;
    fixedOrigin[2] = originVal;

    this->GetBSplineOutput()->SetTransformDomainOrigin(fixedOrigin);

    typename BSplineType::MeshSizeType meshSize;
    meshSize.Fill(this->GetBSplineControlPoints() - NBSplineOrder);
    this->GetBSplineOutput()->SetTransformDomainMeshSize(meshSize);

    typename BSplineType::PhysicalDimensionsType fixedBSplinePhysicalDimensions;
    typename ImageType::SpacingType  spacing;
    spacing[0] = 1;
    spacing[1] = 1;
    spacing[2] = 1;

    for( size_t i = 0; i < NDimension; ++i)
      {
      fixedBSplinePhysicalDimensions[i] = spacing[i] * size[i];
      }
    this->GetBSplineOutput()->SetTransformDomainPhysicalDimensions(fixedBSplinePhysicalDimensions);

    //Get the number of paramaters/nodes required for this BSpline
    const unsigned int numberOfParameters = this->GetBSplineOutput()->GetNumberOfParameters();

    //Setup a paramaters variable for the bspline
    typename BSplineType::ParametersType bSplineParams( numberOfParameters );

    //set up random control point creation
    std::srand(std::time(nullptr));

    ImagePointer coefficientImgLR = this->GetBSplineOutput()->GetCoefficientImages()[0];
    ImagePointer coefficientImgPA = this->GetBSplineOutput()->GetCoefficientImages()[1];
    ImagePointer coefficientImgSI = this->GetBSplineOutput()->GetCoefficientImages()[2];

    typedef typename itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
    IteratorType LRit(coefficientImgLR, coefficientImgLR->GetLargestPossibleRegion());
    IteratorType PAit(coefficientImgPA, coefficientImgPA->GetLargestPossibleRegion());
    IteratorType SIit(coefficientImgSI, coefficientImgSI->GetLargestPossibleRegion());

    LRit.GoToBegin();
    PAit.GoToBegin();
    SIit.GoToBegin();

    if( m_Debug )
      {
      std::cout << "File:  " << __FILE__ << std::endl;
      std::cout << "Line:  " << __LINE__ << std::endl;
      std::cout << "In function GenerateData()" << std::endl;

      std::cout << "inputToBSpline" << std::endl;
      this->GetInput()->Print(std::cout, 0);
      std::cout << "done Input" << std::endl;

      std::cout << "LR coefficient image" << std::endl;
      coefficientImgLR->Print(std::cout, 0);
      std::cout << std::endl << "doneLR" << std::endl;

      std::cout << "PA coefficient image" << std::endl;
      coefficientImgPA->Print(std::cout, 0);
      std::cout << std::endl << "donePA" << std::endl;

      std::cout << std::endl << "SI coefficient image" << std::endl;
      coefficientImgSI->Print(std::cout, 0);
      std::cout << std::endl << "doneSI" << std::endl;
      }

    LRit.GoToBegin();
    PAit.GoToBegin();
    SIit.GoToBegin();
    for( ; !LRit.IsAtEnd(); ++LRit, ++PAit, ++SIit)
      {
      ImagePointType pointLR;
      ImagePointType pointPA;
      ImagePointType pointSI;
      coefficientImgLR->TransformIndexToPhysicalPoint(LRit.GetIndex(), pointLR);
      coefficientImgPA->TransformIndexToPhysicalPoint(PAit.GetIndex(), pointPA);
      coefficientImgSI->TransformIndexToPhysicalPoint(SIit.GetIndex(), pointSI);

      double x = pointLR[0];

      //if( y > 0 ) continue; //only front half of skull
      if( x > 0 )
        {
        LRit.Set( myRandom() );
        PAit.Set( myRandom() );//might have to do these seperateley(different loop) for each one if there is no guarentee of order
        SIit.Set( myRandom() );
        }
      }


    LRit.GoToBegin();
    PAit.GoToBegin();
    SIit.GoToBegin();
    for( ; !LRit.IsAtEnd(); ++LRit, ++PAit, ++SIit)
      {

      ImagePointType pointLR;
      ImagePointType pointPA;
      ImagePointType pointSI;
      coefficientImgLR->TransformIndexToPhysicalPoint(LRit.GetIndex(), pointLR);
      coefficientImgPA->TransformIndexToPhysicalPoint(PAit.GetIndex(), pointPA);
      coefficientImgSI->TransformIndexToPhysicalPoint(SIit.GetIndex(), pointSI);

      double x = pointLR[0];
      double y = pointPA[1];
      double z = pointSI[2];

      //if( y > 0 ) continue; //Only front half of skull
      if( x <= 0 )
        {
        ImagePointType pos;
        pos[0] = x * -1;
        pos[1] = y;
        pos[2] = z;
        typename ImageType::IndexType idx;
        coefficientImgLR->TransformPhysicalPointToIndex(pos, idx);
        //set the weight along LR direction to be anti symmetric and other direction symmetric
        LRit.Set( coefficientImgLR->GetPixel(idx) * -1 );
        PAit.Set( coefficientImgPA->GetPixel(idx) ); //might have to do these seperateley for each one if there is no guarentee of order
        SIit.Set( coefficientImgSI->GetPixel(idx) );
        }
      }

    //We shouldn't need to set the parameters since we set the coefficient images
    typedef typename BSplineType::CoefficientImageArray CoefficeintImages;
    CoefficeintImages images;
    images[0] = coefficientImgLR;
    images[1] = coefficientImgPA;
    images[2] = coefficientImgSI;
    this->GetBSplineOutput()->SetCoefficientImages(images);
  }

private:
  inline double myRandom()
  {
    int min = this->GetRandMin();
    int max = this->GetRandMax();
    int range = max - min;
    return static_cast< double >( rand() % range +1 - max ) * this->GetRandScale();
  }

  BSplinePointer m_BSplineOutput;
  unsigned int m_BSplineControlPoints;
  int m_RandMin;
  int m_RandMax;
  double m_RandScale;
  bool             m_Verbose;
  bool             m_Debug;
};
#endif //BRAINSTOOLS_CREATERANDOMBSPLINE_H
