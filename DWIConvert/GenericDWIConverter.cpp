//
// Created by Hui Xie on 12/19/16.
//
#include "GenericDWIConverter.h"

GenericDWIConverter::GenericDWIConverter( DWIConverter::FileNamesContainer &inputFileNames , const bool FSLFileFormatHorizontalBy3Rows)
:  DWIConverter( inputFileNames, FSLFileFormatHorizontalBy3Rows)
{
}

void GenericDWIConverter::LoadFromDisk() ITK_OVERRIDE
{
  itkGenericExceptionMacro(<< " LoadFromDisk not relevant" << std::endl);
}

GenericDWIConverter::~GenericDWIConverter() {}

void GenericDWIConverter::ExtractDWIData() ITK_OVERRIDE
{
  itkGenericExceptionMacro(<< " ExtractDWIData not relevant" << std::endl);
}

void GenericDWIConverter::AddFlagsToDictionary() ITK_OVERRIDE
{
  itkGenericExceptionMacro(<< " AddFlagsToDictionary not relevant" << std::endl);
}

/**
* @brief Return common fields.  Does nothing for FSL
* @return empty map
*/
DWIConverter::CommonDicomFieldMapType GenericDWIConverter::GetCommonDicomFieldsMap() const ITK_OVERRIDE
{
  return CommonDicomFieldMapType();
}