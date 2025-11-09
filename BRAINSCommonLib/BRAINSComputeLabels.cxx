#include "BRAINSComputeLabels.h"


LabelCountMapType
GetMinLabelCount(ByteImageType::Pointer & labelsImage, const vnl_vector<unsigned int> & PriorLabelCodeVector)
{
  using LabelStatisticsImageFilterType = itk::LabelStatisticsImageFilter<ByteImageType, ByteImageType>;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput(labelsImage);
  labelStatisticsImageFilter->SetInput(labelsImage);
  labelStatisticsImageFilter->Update();

  LabelCountMapType labelCountMap;

  for (unsigned int iclass : PriorLabelCodeVector)
  {
    const auto labelValue = static_cast<size_t>(iclass);
    if (labelStatisticsImageFilter->HasLabel(labelValue))
    {
      const size_t currentLabelCount = labelStatisticsImageFilter->GetCount(labelValue);
      // std::cout << "label: " << (size_t)labelValue << " count: " << currentLabelCount << std::endl;
      labelCountMap[labelValue] = currentLabelCount;
    }
    else
    {
      // std::cout << "label: " << (size_t)labelValue << " count: 0" << std::endl;
      labelCountMap[labelValue] = 0;
    }
  }
  for (auto & it : labelCountMap)
  {
    std::cout << "label: " << static_cast<size_t>(it.first) << " count: " << it.second << std::endl;
  }
  return labelCountMap;
}
