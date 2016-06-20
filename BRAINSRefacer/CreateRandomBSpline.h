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
  };
  ~CreateRandomBSpline(){};

  void GenerateData() ITK_OVERRIDE
  {
    std::cout << "Hello From CreateRandomBSpline!!!" << std::endl;

#if 0 //old way of initializing BSpline
    typename BSplineType::MeshSizeType meshSize;                  //Setup a mesh that contains the number of controlpoints
    meshSize.Fill(m_BSplineControlPoints - NBSplineOrder);         //Inspired from itk example "BSplineWarping2.cxx"

    this->GetBSplineOutput()->SetTransformDomainMeshSize(meshSize);            //TODO: ask if it is possible to have "uneven" control points.
                                                                // EG. 8 on LR axis 7 on SI 6 on AP. We did this in simple itk
                                                                //so it should be possible in cpp itk

    typedef typename ImageType::RegionType ImageRegionType;
    ImageRegionType subjectRegion = this->GetInput()->GetLargestPossibleRegion();


    this->GetBSplineOutput()->SetTransformDomainOrigin(this->GetInput()->GetOrigin());           //Origin
    this->GetBSplineOutput()->SetTransformDomainDirection(this->GetInput()->GetDirection());     //Direction
    this->GetBSplineOutput()->SetTransformDomainPhysicalDimensions((                    //PhysicalDimensions
                                                    this->GetInput()->GetSpacing()[0]*(subjectRegion.GetSize()[0]-1),         //Should all be set to the same as the subject Image
                                                    this->GetInput()->GetSpacing()[1]*(subjectRegion.GetSize()[1]-1),
                                                    this->GetInput()->GetSpacing()[2]*(subjectRegion.GetSize()[2]-1)
                                                  ));
#endif
    typedef typename ImageType::RegionType ImageRegionType;
    ImageRegionType subjectRegion = this->GetInput()->GetLargestPossibleRegion();

    //one new way of initializing bSPline
    //this method was inspired by itk example for deformable registration15
    typename BSplineType::PhysicalDimensionsType fixedBSplinePhysicalDimensions;
    typename BSplineType::MeshSizeType meshSize;
    typename BSplineType::OriginType  fixedOrigin;

    for( size_t i = 0; i < NDimension; ++i)
      {
      fixedOrigin[i] = this->GetInput()->GetOrigin()[i];
      fixedBSplinePhysicalDimensions[i] = this->GetInput()->GetSpacing()[i] *
        static_cast<double>(subjectRegion.GetSize()[i]-1);

      }

    meshSize.Fill(this->GetBSplineControlPoints() - NBSplineOrder);
    this->GetBSplineOutput()->SetTransformDomainOrigin(fixedOrigin);
    this->GetBSplineOutput()->SetTransformDomainPhysicalDimensions(fixedBSplinePhysicalDimensions);
    this->GetBSplineOutput()->SetTransformDomainMeshSize(meshSize);
    this->GetBSplineOutput()->SetTransformDomainDirection(this->GetInput()->GetDirection() );


    //Get the number of paramaters/nodes required for this BSpline
    const unsigned int numberOfParameters = this->GetBSplineOutput()->GetNumberOfParameters();

    //Setup a paramaters variable for the bspline
    typename BSplineType::ParametersType bSplineParams( numberOfParameters );

    std::srand(time(nullptr));

#if 1
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

    std::cout<< "inputToBSpline" <<std::endl;
    this->GetInput()->Print(std::cout,0);
    std::cout<<"done Input" << std::endl;

    std::cout<<"LR"<<std::endl;
    coefficientImgLR->Print(std::cout,0);
    std::cout<<std::endl<<"doneLR"<<std::endl;

    std::cout<<"PA"<<std::endl;
    coefficientImgPA->Print(std::cout,0);
    std::cout<<std::endl<<"donePA"<<std::endl;

    std::cout<<std::endl<<"SI"<<std::endl;
    coefficientImgSI->Print(std::cout,0);
    std::cout<<std::endl<<"doneSI"<<std::endl;

    // assume spacing, origin, IndexToPointMatrix, PointToIndexMatrix
    // is the same between all coefficient images

    ///std::cout<<"-------------------" << LRit.GetIndex() << std::endl;



    LRit.GoToBegin();
    PAit.GoToBegin();
    SIit.GoToBegin();
    for( ; !LRit.IsAtEnd(); ++LRit, ++PAit, ++SIit)
      {
#if 0 //old method
      ImagePointType point;
      coefficientImgLR->TransformIndexToPhysicalPoint(LRit.GetIndex(), point);
#endif
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

#if 0
      ImagePointType point;
      coefficientImgLR->TransformIndexToPhysicalPoint(LRit.GetIndex(), point);

      double x = point[0];
      double y = point[1];
      double z = point[2];
#endif
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

    //std::cout << this->GetBSplineOutput()->GetParameters() <<std::endl;

    //There should be a better way to do this:

#endif


    //Old method below
#if 0
    for( unsigned int n = 0; n < numberOfNodes; ++ n)
      {
      bSplineParams[n] = myRandom();
      bSplineParams[n + numberOfNodes] = myRandom();      // "y" coord;
      bSplineParams[n + numberOfNodes * 2] = myRandom();  // "z" coord;
      // TODO: x,y,z seem like they are the wrong coordinate system. Get a better model
      }


    this->GetBSplineOutput()->SetParameters(bSplineParams);
    std::cout << this->GetBSplineOutput()->GetParameters() <<std::endl;
#endif
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
};
#endif //BRAINSTOOLS_CREATERANDOMBSPLINE_H
