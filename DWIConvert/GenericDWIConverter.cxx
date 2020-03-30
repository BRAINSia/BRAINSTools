//
// Created by Hui Xie on 12/19/16.
//
#include "GenericDWIConverter.h"

GenericDWIConverter::GenericDWIConverter(DWIConverter::FileNamesContainer & inputFileNames)
  : DWIConverter(inputFileNames)
{}

void
GenericDWIConverter::LoadFromDisk()
{
  itkGenericExceptionMacro(<< " LoadFromDisk not relevant" << std::endl);
}

GenericDWIConverter::~GenericDWIConverter() = default;

void
GenericDWIConverter::ExtractDWIData()
{
  itkGenericExceptionMacro(<< " ExtractDWIData not relevant" << std::endl);
}

void
GenericDWIConverter::AddFlagsToDictionary()
{
  itkGenericExceptionMacro(<< " AddFlagsToDictionary not relevant" << std::endl);
}

/**
 * @brief Return common fields.  Does nothing for FSL
 * @return empty map
 */
DWIConverter::CommonDicomFieldMapType
GenericDWIConverter::GetCommonDicomFieldsMap() const
{
  return CommonDicomFieldMapType();
}
