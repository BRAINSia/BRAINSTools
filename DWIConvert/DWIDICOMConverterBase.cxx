//
// Created by Johnson, Hans J on 11/24/16.
//

#include "DWIDICOMConverterBase.h"

/**
 * @brief Return common fields.  Does nothing for FSL
 * @return empty map
 */
DWIDICOMConverterBase::CommonDicomFieldMapType
DWIDICOMConverterBase::GetCommonDicomFieldsMap() const
{
return this->m_CommonDicomFieldsMap;
}
