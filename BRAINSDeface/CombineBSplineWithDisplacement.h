//
// Created by Alexander Leinoff on 6/1/16.
//

#ifndef BRAINSTOOLS_COMBINEBSPLINEWITHDISPLACEMENT_H
#define BRAINSTOOLS_COMBINEBSPLINEWITHDISPLACEMENT_H

#include <itkMultiplyImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>

template<
  typename TReferenceImageType,
  typename TOutputVectorImageType,
  typename PixelType,
  unsigned int NDimension,
  unsigned int NBSplineOrder>
class CombineBSplineWithDisplacement
  : public itk::ImageToImageFilter<TReferenceImageType, TOutputVectorImageType>
{
public:
  typedef CombineBSplineWithDisplacement Self;
  typedef itk::SmartPointer< Self > Pointer;

  itkNewMacro(Self);
  itkTypeMacro(Self, ImageToImageFilter)

  typedef TReferenceImageType ImageType;
  typedef typename ImageType::Pointer ImageTypePointer;

  typedef itk::BSplineTransform< PixelType,NDimension,NBSplineOrder> BSplineType;
  typedef typename BSplineType::Pointer BSplinePointer;

  typedef itk::Vector<PixelType, NDimension > VectorPixelType;
  typedef itk::Image< VectorPixelType, NDimension> DisplacementFieldImageType;
  typedef typename DisplacementFieldImageType::Pointer DisplacementFieldPointer;



  itkSetMacro(BSplineInput, BSplinePointer)
  itkGetMacro(BSplineInput, BSplinePointer)

  itkSetMacro(DistanceMap, ImageTypePointer )
  itkGetMacro(DistanceMap, ImageTypePointer )

  itkSetMacro(ComposedImage, DisplacementFieldPointer)
  itkGetMacro(ComposedImage, DisplacementFieldPointer)


protected:
  CombineBSplineWithDisplacement(){};
  ~CombineBSplineWithDisplacement(){};

  void GenerateData() ITK_OVERRIDE
  {
    std::cout << "In Generate data method of CombineRandomBSplineWithDisplacement" << std::endl;

    const TReferenceImageType * subject = this->GetInput();
    BSplinePointer bSpline = this->GetBSplineInput();

    typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, PixelType> TransformToDisplacementFilterType;
    typename TransformToDisplacementFilterType::Pointer bSplineDisplacementFieldGenerator = TransformToDisplacementFilterType::New();
    bSplineDisplacementFieldGenerator->UseReferenceImageOn();
    bSplineDisplacementFieldGenerator->SetReferenceImage(subject);
    bSplineDisplacementFieldGenerator->SetTransform(bSpline);

    std::cout<<"Extracting component images from displacement field"<<std::endl;
    //multiply the displacement field by the distance map to get the "smooth displacement that doesn't affect the brain
    //first extrace scalar elemnts from vector image
    typedef typename itk::VectorIndexSelectionCastImageFilter<DisplacementFieldImageType, ImageType> ImageExtractionFilterType;
    typename ImageExtractionFilterType::Pointer xTractDisplacementFilter = ImageExtractionFilterType::New();
    typename ImageExtractionFilterType::Pointer yTractDisplacementFilter = ImageExtractionFilterType::New();
    typename ImageExtractionFilterType::Pointer zTractDisplacementFilter = ImageExtractionFilterType::New();

    xTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
    xTractDisplacementFilter->SetIndex(0);
    typename ImageType::Pointer   xDisplacement = xTractDisplacementFilter->GetOutput();

    yTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
    yTractDisplacementFilter->SetIndex(1);
    typename ImageType::Pointer  yDisplacement = yTractDisplacementFilter->GetOutput();

    zTractDisplacementFilter->SetInput(bSplineDisplacementFieldGenerator->GetOutput());
    zTractDisplacementFilter->SetIndex(2);
    typename ImageType::Pointer  zDisplacement = zTractDisplacementFilter->GetOutput();

    //multiply by distancemap

    std::cout<<"Multiplying bSplineDisplacement by Distance map"<<std::endl;

    typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplyFilterType;
    typename MultiplyFilterType::Pointer xMult = MultiplyFilterType::New();
    typename MultiplyFilterType::Pointer yMult = MultiplyFilterType::New();
    typename MultiplyFilterType::Pointer zMult = MultiplyFilterType::New();

    ImageTypePointer distanceMapImage = this->GetDistanceMap();

    xMult->SetInput1(xDisplacement);
    xMult->SetInput2(distanceMapImage);
    typename ImageType::Pointer xMultImage = xMult->GetOutput();
    xMult->Update();

    yMult->SetInput1(yDisplacement);
    yMult->SetInput2(distanceMapImage);
    typename ImageType::Pointer yMultImage = yMult->GetOutput();
    yMult->Update();

    zMult->SetInput1(zDisplacement);
    zMult->SetInput2(distanceMapImage);
    typename ImageType::Pointer zMultImage = zMult->GetOutput();
    zMult->Update();

    std::cout<<"Composing new image from displacement and bSpline product components"<<std::endl;

    typedef itk::ComposeImageFilter<ImageType, DisplacementFieldImageType> ComposeFilterType;
    typename ComposeFilterType::Pointer composeDisplacements = ComposeFilterType::New();
    composeDisplacements->SetInput(0, xMult->GetOutput());
    composeDisplacements->SetInput(1, yMult->GetOutput());
    composeDisplacements->SetInput(2, zMult->GetOutput());


    this->SetComposedImage(composeDisplacements->GetOutput());
    //typename DisplacementFieldImageType::Pointer bSplineDistanceMapCombination = composeDisplacements->GetOutput();
    composeDisplacements->Update();
/*
   //This is what I want to do but i end up iterating this->SetOutput( composeDisplacements->GetOutput() );
    //There has got to be a better way to do this
    std::cout<<"Finished generating data for combinebsplinewithdisplacement"<<std::endl;
    //typedef itk::ImageIteratorWithIndex<DisplacementFieldImageType> ConstIter;
    typedef itk::ImageIteratorWithIndex<DisplacementFieldImageType> Iter;

    DisplacementFieldPointer composedImage = composeDisplacements->GetOutput();
    typedef typename DisplacementFieldImageType::RegionType OutputRegionType;
    OutputRegionType outputRegion = composedImage->GetLargestPossibleRegion(); //Getbuffered largestpossible
    DisplacementFieldPointer outputImage = this->GetOutput();
    composeDisplacements->Update();

    //outputImage->SetRegions(outputRegion);
    //outputImage->Allocate();

    Iter composeIter( composedImage, outputRegion);
    //Iter outputIter( outputImage, outputImage->GetLargestPossibleRegion());

    composeIter.GoToBegin();
    //outputIter.GoToBegin();

    while (!composeIter.IsAtEnd())
      {
    //  outputIter.Set(composeIter.Get());
      ++composeIter;
   //   ++outputIter;
      }
*/
  }

private:
  BSplinePointer m_BSplineInput;
  ImageTypePointer m_DistanceMap;
  DisplacementFieldPointer m_ComposedImage;

};
#endif //BRAINSTOOLS_COMBINEBSPLINEWITHDISPLACEMENT_H