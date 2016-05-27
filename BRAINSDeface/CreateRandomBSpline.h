//
// Created by Jeffrey Obadal, Alexander Leinoff on 5/31/16.
//

#include <itkBSplineTransform.h>

#ifndef BRAINSTOOLS_CREATERANDOMBSPLINE_H
#define BRAINSTOOLS_CREATERANDOMBSPLINE_H


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
  //itkGetMacro(BSplineOutput, BSplinePointer)


  itkSetMacro(BSplineControlPoints, unsigned int)
  itkGetMacro(BSplineControlPoints,  unsigned int)

  typedef TInputImage ImageType;

protected:
  CreateRandomBSpline()
  {
    m_BSplineOutput=BSplineType::New();
  };
  ~CreateRandomBSpline(){};
  //void PrintSelf( std::ostream & os, itk::Indent indent ) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE
  {
    std::cout << "Hello From ctreateRandomBSpline!!!" << std::endl;
    //BSplinePointer myBSpline = BSplineType::New();

    typename BSplineType::MeshSizeType meshSize;                  //Setup a mesh that contains the number of controlpoints
    meshSize.Fill(m_BSplineControlPoints - NBSplineOrder);         //Inspired from itk example "BSplineWarping2.cxx"

    this->GetBSplineOutput()->SetTransformDomainMeshSize(meshSize);            //TODO: ask if it is possible to have "uneven" control points.
                                                                // EG. 8 on LR axis 7 on SI 6 on AP. We did this in simple itk
                                                                //so it should be possible in cpp itk

    typedef typename ImageType::RegionType ImageRegionType;
//  ImageRegionType subjectRegion = subject->GetLargestPossibleRegion();
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

    //print out
   // std::cout << "Number of params: " << numberOfParameters << std::endl;
   // std::cout << "Number of nodes:  " << numberOfNodes << std::endl;

    //Setup a paramaters variable for the bspline
    typename BSplineType::ParametersType bSplineParams( numberOfParameters );

    //  From ITK Example "BSplineWarping2"
    //  The B-spline grid should now be fed with coeficients at each node. Since
    //  this is a two dimensional grid, each node should receive two coefficients.
    //  Each coefficient pair is representing a displacement vector at this node.
    //  The coefficients can be passed to the B-spline in the form of an array where
    //  the first set of elements are the first component of the displacements for
    //  all the nodes, and the second set of elemets is formed by the second
    //  component of the displacements for all the nodes.

    //  In the ITK Example, the read the points in from a file. Here, they will be
    //  generated randomly. This should put the xyz coordinates in the correct space
    //  a better way would probably be to use the image coefficient array of the bspline,
    //  but for now I will use the method from the ITK example

    std::srand(time(nullptr));

    for( unsigned int n = 0; n < numberOfNodes; ++ n)
      {
      bSplineParams[n] = myRandom();
      bSplineParams[n + numberOfNodes] = myRandom();      // "y" coord;
      bSplineParams[n + numberOfNodes * 2] = myRandom();  // "z" coord;
      // TODO: x,y,z seem like they are the wrong coordinate system. Get a better model
      }

    this->GetBSplineOutput()->SetParameters(bSplineParams);

    //myBSpline->Print(std::cout,3);

    //this->GetBSplineOutput()=myBSpline;
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