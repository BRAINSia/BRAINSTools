#ifndef __nrrdCommon_h__
#define __nrrdCommon_h__
#include "mex.h"
#include "matrix.h"
#include "string.h"
// #include <teem/nrrd.h>
// #include <teem/ten.h>
#include <NrrdIO.h> //FROM ITK
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <stdio.h>

#include "Debug.h"

#define NRRD_MAX_ERROR_MSG_SIZE 8095

#define FIELDNAME_INDEX_data 0
#define FIELDNAME_INDEX_space 1
#define FIELDNAME_INDEX_spacedirections 2
#define FIELDNAME_INDEX_centerings 3
#define FIELDNAME_INDEX_kinds 4
#define FIELDNAME_INDEX_spaceunits 5
#define FIELDNAME_INDEX_spacedefinition 6
#define FIELDNAME_INDEX_spaceorigin 7
#define FIELDNAME_INDEX_measurementframe 8
#define FIELDNAME_INDEX_modality 9
#define FIELDNAME_INDEX_bvalue 10
#define FIELDNAME_INDEX_gradientdirections 11
#define FIELDNAME_INDEX_MANDATORYFEILDS 9
#define FIELDNAME_INDEX_MAXFEILDS 12

// field_num = mxGetFieldNumber(pa, "field_name");
// mxGetFieldByNumber(pa, index, field_num);

class MatlabStructManager
{
public:
  typedef std::map<std::string, int> FieldMapType;
  MatlabStructManager(const mxArray * const structMx)
    : m_structMx(structMx)
  {
    std::fill( m_errBuff, m_errBuff + NRRD_MAX_ERROR_MSG_SIZE, '\0' );

    m_StructMap[std::string("data")] = 0;
    m_StructMap[std::string("space")] = 1;
    m_StructMap[std::string("spacedirections")] = 2;
    m_StructMap[std::string("centerings")] = 3;
    m_StructMap[std::string("kinds")] = 4;
    m_StructMap[std::string("spaceunits")] = 5;
    m_StructMap[std::string("spacedefinition")] = 6;
    m_StructMap[std::string("spaceorigin")] = 7;
    m_StructMap[std::string("measurementframe")] = 8;
    m_StructMap[std::string("modality")] = 9;
    m_StructMap[std::string("bvalue")] = 10;
    m_StructMap[std::string("gradientdirections")] = 11;
    this->isValidStruct();
  }

  ~MatlabStructManager()
  {
  }

  int GetNumberOfFields() const
  {
    return mxGetNumberOfFields(m_structMx);
  }

  bool isDataComplex() const
  {
    return mxIsComplex( this->GetField("data") );
  }

  /**
    * Verify that given field name is a valid for
    * the matlab itk/nrrd interface interactions.
    */
  bool isValidField( const std::string & name) const
  {
    const int match_count = m_StructMap.count(name);

    if( match_count == 0 )
      {
      snprintf(m_errBuff, NRRD_MAX_ERROR_MSG_SIZE, "ERROR: Invalid structure element requested %s\n", name.c_str() );
      return false;
      }
    return true;
  }

  bool structHasField( const std::string & name) const
  {
    const mxArray * const temp_field = mxGetField(this->m_structMx, 0, name.c_str() );

    if( isValidField(name) && temp_field )
      {
      return true;
      }
    return false;
  }

  /**
   * Verify that this matlab struct only contains valid structure names.
   */
  bool isValidStruct() const
  {
    for( int field_index = 0; field_index < this->GetNumberOfFields(); ++field_index )
      {
      const char * c_name = mxGetFieldNameByNumber(m_structMx, field_index);
      if( c_name == 0 )
        {
        continue;
        }
      const std::string name(c_name);
      if( !isValidField(name) )
        {
        return false;
        }
      }
    /* Error checking on the data */
    if( this->isDataComplex() )
      {
      snprintf(m_errBuff, NRRD_MAX_ERROR_MSG_SIZE, "ERROR: sorry, array must be real\n");
      mexErrMsgTxt(m_errBuff);
      return false;
      }
    if( nrrdTypeUnknown == this->GetDataType() )
      {
      snprintf(m_errBuff, NRRD_MAX_ERROR_MSG_SIZE, "ERROR: Sorry, can't handle type %s", mxGetClassName(GetField(
                                                                                                          "data") ) );
      mexErrMsgTxt(m_errBuff);
      return false;
      }
    return true;
  }

  /**
   * Determine if this is dwi data.
   */
  bool isDWIdata() const
  {
    if( this->structHasField("modality") )
      {
      const mxArray * const modalityField = this->GetField("modality");
      if( modalityField )
        {
        char modalityCharArray[NRRD_MAX_ERROR_MSG_SIZE];
        mxGetString(modalityField,  modalityCharArray, NRRD_MAX_ERROR_MSG_SIZE - 1);
        if( std::string("DWMRI") == std::string(modalityCharArray) )
          {
          return true;
          }
        }
      if( this->GetNumberOfDimensions("data") == 4 ) // Assume DWI if dim is 4
        {
        return true;
        }
      }
    return false;
  }

#if 0
  /*TODO*/ snprintf(m_errBuff, NRRD_MAX_ERROR_MSG_SIZE, "HERE %d: %s ", __LINE__, __FILE__); mexWarnMsgTxt(m_errBuff);
#endif
  const mwSize * GetDimensions( const std::string & name ) const
  {
    return mxGetDimensions( this->GetField(name) );
  }

  unsigned int GetNumberOfDimensions(const std::string & name) const
  {
    return static_cast<unsigned int>( mxGetNumberOfDimensions( this->GetField( name ) ) );
  }

  int GetDataType() const
  {
    const mxClassID mtype = mxGetClassID( this->GetField("data") );

    return typeMtoN(mtype);
  }

  mxArray * GetField(const std::string & name) const
  {
    const int index = 0; // Always 0 in our case

    if( !isValidField(name) )
      {
      snprintf(m_errBuff, NRRD_MAX_ERROR_MSG_SIZE, "WARNING: sorry, can not access %s\n", name.c_str() );
      mexErrMsgTxt(m_errBuff);
      return NULL;
      }
    const int field_num = mxGetFieldNumber(m_structMx, name.c_str() );
    return mxGetFieldByNumber(m_structMx, index, field_num);
  }

private:

  int typeMtoN(const mxClassID mtype) const
  {
    int ntype;

    switch( mtype )
      {
      case mxINT8_CLASS:
        ntype = nrrdTypeChar;
        break;
      case mxUINT8_CLASS:
        ntype = nrrdTypeUChar;
        break;
      case mxINT16_CLASS:
        ntype = nrrdTypeShort;
        break;
      case mxUINT16_CLASS:
        ntype = nrrdTypeUShort;
        break;
      case mxINT32_CLASS:
        ntype = nrrdTypeInt;
        break;
      case mxUINT32_CLASS:
        ntype = nrrdTypeUInt;
        break;
      case mxINT64_CLASS:
        ntype = nrrdTypeLLong;
        break;
      case mxUINT64_CLASS:
        ntype = nrrdTypeULLong;
        break;
      case mxSINGLE_CLASS:
        ntype = nrrdTypeFloat;
        break;
      case mxDOUBLE_CLASS:
        ntype = nrrdTypeDouble;
        break;
      default:
        ntype = nrrdTypeUnknown;
        break;
      }
    return ntype;
  }

  FieldMapType          m_StructMap;
  mutable char          m_errBuff[NRRD_MAX_ERROR_MSG_SIZE];
  const mxArray * const m_structMx;
};

static const char GRADIENT_PREFIX[] = "DWMRI_gradient_";

extern mxClassID typeNtoM(const int ntype);

extern bool isGradientAxis(const unsigned int kindFlag);

#endif /* __nrrdCommon_h__ */
