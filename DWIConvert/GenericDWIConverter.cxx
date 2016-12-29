//
// Created by Hui Xie on 12/19/16.
//
#include "GenericDWIConverter.h"

GenericDWIConverter::GenericDWIConverter( DWIConverter::FileNamesContainer &inputFileNames , const bool FSLFileFormatHorizontalBy3Rows)
:  DWIConverter( inputFileNames, FSLFileFormatHorizontalBy3Rows)
{
}

void GenericDWIConverter::LoadFromDisk()
{
  itkGenericExceptionMacro(<< " LoadFromDisk not relevant" << std::endl);
}

GenericDWIConverter::~GenericDWIConverter() {}

void GenericDWIConverter::ExtractDWIData()
{
  itkGenericExceptionMacro(<< " ExtractDWIData not relevant" << std::endl);
}

void GenericDWIConverter::AddFlagsToDictionary()
{
  itkGenericExceptionMacro(<< " AddFlagsToDictionary not relevant" << std::endl);
}

/**
* @brief Return common fields.  Does nothing for FSL
* @return empty map
*/
DWIConverter::CommonDicomFieldMapType GenericDWIConverter::GetCommonDicomFieldsMap() const
{
  return CommonDicomFieldMapType();
}
