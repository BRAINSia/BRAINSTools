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
  using Self = CombineBSplineWithDisplacement;
  using Pointer = itk::SmartPointer< Self >;

  itkNewMacro(Self);
  itkTypeMacro(Self, ImageToImageFilter)

  using ImageType = TReferenceImageType;
  using ImageTypePointer = typename ImageType::Pointer;

  using BSplineType = itk::BSplineTransform< PixelType,NDimension,NBSplineOrder>;
  using BSplinePointer = typename BSplineType::Pointer;

  using VectorPixelType = itk::Vector<PixelType, NDimension >;
  using DisplacementFieldImageType = itk::Image< VectorPixelType, NDimension>;
  using DisplacementFieldPointer = typename DisplacementFieldImageType::Pointer;

  itkSetMacro(Verbose, bool)
  itkSetMacro(Debug, bool)

  itkSetMacro(BSplineInput, BSplinePointer)
  itkGetMacro(BSplineInput, BSplinePointer)

  itkSetMacro(DistanceMap, ImageTypePointer )
  itkGetMacro(DistanceMap, ImageTypePointer )

  itkSetMacro(ComposedImage, DisplacementFieldPointer)
  itkGetMacro(ComposedImage, DisplacementFieldPointer)

protected:
  CombineBSplineWithDisplacement()
  {
    this->m_Verbose = false;
    this->m_Debug = false;
  }
  ~CombineBSplineWithDisplacement() override{};

  void GenerateData() override
  {
    if( m_Debug )
      {
        std::cout << "File:  " << __FILE__ << std::endl;
        std::cout << "Line:  " << __LINE__ << std::endl;
        std::cout << "In Generate data method of CombineRandomBSplineWithDisplacement" << std::endl;
      }

    const TReferenceImageType * subject = this->GetInput();
    BSplinePointer bSpline = this->GetBSplineInput();

    using TransformToDisplacementFilterType = itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, PixelType>;
    typename TransformToDisplacementFilterType::Pointer bSplineDisplacementFieldGenerator = TransformToDisplacementFilterType::New();
    bSplineDisplacementFieldGenerator->UseReferenceImageOn();
    bSplineDisplacementFieldGenerator->SetReferenceImage(subject);
    bSplineDisplacementFieldGenerator->SetTransform(bSpline);

    if( m_Debug )
      {
      std::cout << "File:  " << __FILE__ << std::endl;
      std::cout << "Line:  " << __LINE__ << std::endl;
      std::cout << "In function GenerateData()" << std::endl;

      std::cout << "Extracting component images from displacement field" << std::endl;
      }
    //multiply the displacement field by the distance map to get the "smooth displacement that doesn't affect the brain
    //first extrace scalar elemnts from vector image
    using ImageExtractionFilterType = typename itk::VectorIndexSelectionCastImageFilter<DisplacementFieldImageType, ImageType>;
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

    if( m_Verbose || m_Debug)
      {
      std::cout<<"Multiplying bSplineDisplacement by Distance map"<<std::endl;
      }

    using MultiplyFilterType = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>;
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

    if( m_Debug || m_Verbose )
      {
      std::cout<<"Composing new image from displacement and bSpline product components"<<std::endl;
      }

    using ComposeFilterType = itk::ComposeImageFilter<ImageType, DisplacementFieldImageType>;
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
    //using ConstIter = itk::ImageIteratorWithIndex<DisplacementFieldImageType>;
    using Iter = itk::ImageIteratorWithIndex<DisplacementFieldImageType>;

    DisplacementFieldPointer composedImage = composeDisplacements->GetOutput();
    using OutputRegionType = typename DisplacementFieldImageType::RegionType;
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

  bool             m_Verbose;
  bool             m_Debug;

};
#endif //BRAINSTOOLS_COMBINEBSPLINEWITHDISPLACEMENT_H
