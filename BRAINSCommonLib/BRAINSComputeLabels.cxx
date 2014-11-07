#include "BRAINSComputeLabels.h"


LabelCountMapType GetMinLabelCount(ByteImageType::Pointer & labelsImage)
{
  // We will choose "KNN_SamplesPerLabel" from each posterior class.
  typedef itk::LabelStatisticsImageFilter<ByteImageType,ByteImageType> LabelStatisticsImageFilterType;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput(labelsImage);
  labelStatisticsImageFilter->SetInput(labelsImage);
  labelStatisticsImageFilter->Update();

  typedef LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;

  LabelCountMapType labelCountMap;

  for(ValidLabelValuesType::const_iterator vIt=labelStatisticsImageFilter->GetValidLabelValues().begin();
    vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
    ++vIt)
    {
    if ( labelStatisticsImageFilter->HasLabel(*vIt) )
      {
      size_t labelValue = static_cast<size_t>(*vIt); //Converting to size_t.  Should be LabelPixelType, but that requires refactoring that I don't have time for
      const size_t currentLabelCount = labelStatisticsImageFilter->GetCount( labelValue );
      std::cout << "label: " << (size_t)labelValue << " count: " << labelStatisticsImageFilter->GetCount( labelValue ) << std::endl;
      labelCountMap[labelValue] = currentLabelCount;
      }
    }
  return labelCountMap;
}
