//
// Created by Jeffrey Obadal, Alexander Leinoff on 5/31/16.
//
#ifndef BRAINSTOOLS_CREATERANDOMBSPLINE_H
#define BRAINSTOOLS_CREATERANDOMBSPLINE_H

#include <itkBSplineTransform.h>

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
  itkGetMacro(BSplineControlPoints,  unsigned int)

  typedef TInputImage ImageType;

protected:
  CreateRandomBSpline()
  {
    m_BSplineOutput=BSplineType::New();
  };
  ~CreateRandomBSpline(){};

  void GenerateData() ITK_OVERRIDE
  {
    std::cout << "Hello From CreateRandomBSpline!!!" << std::endl;

    typename BSplineType::MeshSizeType meshSize;                  //Setup a mesh that contains the number of controlpoints
    meshSize.Fill(m_BSplineControlPoints - NBSplineOrder);         //Inspired from itk example "BSplineWarping2.cxx"

    this->GetBSplineOutput()->SetTransformDomainMeshSize(meshSize);            //TODO: ask if it is possible to have "uneven" control points.
                                                                // EG. 8 on LR axis 7 on SI 6 on AP. We did this in simple itk
                                                                //so it should be possible in cpp itk

    typedef typename ImageType::RegionType ImageRegionType;
    ImageRegionType subjectRegion = this->GetInput()->GetBufferedRegion();

    this->GetBSplineOutput()->SetTransformDomainOrigin(this->GetInput()->GetOrigin());           //Origin
    this->GetBSplineOutput()->SetTransformDomainDirection(this->GetInput()->GetDirection());     //Direction
    this->GetBSplineOutput()->SetTransformDomainPhysicalDimensions((                    //PhysicalDimensions
                                                    this->GetInput()->GetSpacing()[0]*(subjectRegion.GetSize()[0]-1),         //Should all be set to the same as the subject Image
                                                      this->GetInput()->GetSpacing()[1]*(subjectRegion.GetSize()[1]-1),
                                                      this->GetInput()->GetSpacing()[2]*(subjectRegion.GetSize()[2]-1)
                                                  ));

    //Get the number of paramaters/nodes required for this BSpline
    const unsigned int numberOfParameters = this->GetBSplineOutput()->GetNumberOfParameters();
    const unsigned int numberOfNodes = numberOfParameters / NDimension;

    //Setup a paramaters variable for the bspline
    typename BSplineType::ParametersType bSplineParams( numberOfParameters );

    std::srand(time(nullptr));

    for( unsigned int n = 0; n < numberOfNodes; ++ n)
      {
      bSplineParams[n] = myRandom();
      bSplineParams[n + numberOfNodes] = myRandom();      // "y" coord;
      bSplineParams[n + numberOfNodes * 2] = myRandom();  // "z" coord;
      // TODO: x,y,z seem like they are the wrong coordinate system. Get a better model
      }

    this->GetBSplineOutput()->SetParameters(bSplineParams);
  }

private:
  inline double myRandom()
  {
    const int min = -5;
    const int max = 5;
    const int range = max - min;
    return static_cast< double >( rand() % range +1 - max )*0.05;
  }
  BSplinePointer m_BSplineOutput;
  unsigned int m_BSplineControlPoints;
};
#endif //BRAINSTOOLS_CREATERANDOMBSPLINE_H