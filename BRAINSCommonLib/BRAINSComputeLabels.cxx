#include "BRAINSComputeLabels.h"


LabelCountMapType GetMinLabelCount(ByteImageType::Pointer & labelsImage,
                                   const vnl_vector<unsigned int> & PriorLabelCodeVector)
{
  typedef itk::LabelStatisticsImageFilter<ByteImageType,ByteImageType> LabelStatisticsImageFilterType;
  LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  labelStatisticsImageFilter->SetLabelInput(labelsImage);
  labelStatisticsImageFilter->SetInput(labelsImage);
  labelStatisticsImageFilter->Update();

  LabelCountMapType labelCountMap;

  for( size_t iclass = 0; iclass < PriorLabelCodeVector.size(); ++iclass )
     {
     size_t labelValue = static_cast<size_t>( PriorLabelCodeVector[iclass] );
     if( labelStatisticsImageFilter->HasLabel(labelValue) )
       {
       const size_t currentLabelCount = labelStatisticsImageFilter->GetCount( labelValue );
       std::cout << "label: " << (size_t)labelValue << " count: " << labelStatisticsImageFilter->GetCount( labelValue ) << std::endl;
       labelCountMap[labelValue] = currentLabelCount;
       }
      else
       {
       std::cout << "label: " << (size_t)labelValue << " count: 0" << std::endl;
       labelCountMap[labelValue] = 0;
       }
     }

  return labelCountMap;
}
