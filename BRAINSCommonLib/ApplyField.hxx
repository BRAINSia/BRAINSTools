#ifndef __ApplyField_hxx
#define __ApplyField_hxx
#include "ApplyField.h"
#include "itkImage.h"
#include "itkWarpImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkReinitializeLevelSetImageFilter.h"
#include "itkIO.h"

namespace itk
{
template <typename TDisplacementField, typename TInputImage,
          typename TOutputImage>
ApplyField<TDisplacementField, TInputImage,
           TOutputImage>::ApplyField() :
  m_InputImage(0),
  m_OutputImage(0),
  m_DisplacementField(0),
  m_DefaultPixelValue(0)
{
}

template <typename TDisplacementField, typename TInputImage,
          typename TOutputImage>
ApplyField<TDisplacementField, TInputImage,
           TOutputImage>::~ApplyField()
{
}

template <typename TDisplacementField, typename TInputImage,
          typename TOutputImage>
void ApplyField<TDisplacementField, TInputImage,
                TOutputImage>::Execute()
{
  if( m_InputImage.IsNull() )
    {
    std::cout << "ERROR:  No Input image give.! " << std::endl;
    }

  typedef WarpImageFilter<InputImageType, OutputImageType,
                          TDisplacementField> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput(m_InputImage);
  warper->SetDeformationField(m_DisplacementField);
  warper->SetOutputParametersFromImage(m_DisplacementField);
  warper->SetEdgePaddingValue(m_DefaultPixelValue);
  warper->Update();
  std::cout << "  Registration Applied" << std::endl;
  m_OutputImage = warper->GetOutput();
}

template <typename TDisplacementField, typename TInputImage,
          typename TOutputImage>
void ApplyField<TDisplacementField, TInputImage,
                TOutputImage>::ReleaseDataFlagOn()
{
  // m_InputImage->DisconnectPipeline();
  // m_DisplacementField->DisconnectPipeline();
}
}
#endif
