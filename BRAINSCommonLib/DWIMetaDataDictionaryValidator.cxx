/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing
 *                     and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
*=========================================================================*/
/*
 *  \authors Hans J. Johnson, David M. Welch, Ali Ghayoor
 *
 *
 * This class is used to build a meta-data dictionary and validate that the
 * meta-data dictionary is suitable to writing valid DWI data to disk in nrrd
 * file format.
 *
 */
#include "DWIMetaDataDictionaryValidator.h"

#include "itkNumberToString.h"

DWIMetaDataDictionaryValidator::DWIMetaDataDictionaryValidator()
{
  this->m_dict.Clear();
}

DWIMetaDataDictionaryValidator::~DWIMetaDataDictionaryValidator() {}

std::string DWIMetaDataDictionaryValidator::GetGradientKeyString(int index) const
{
  char               tmpStr[64];
  sprintf(tmpStr, "DWMRI_gradient_%04u", index);
  std::string key(tmpStr);
  return key;
}


void DWIMetaDataDictionaryValidator::SetStringDictObject(const std::string key, const std::string value)
{
  itk::EncapsulateMetaData<std::string>(m_dict, key, value);
}

DWIMetaDataDictionaryValidator::MetaDataDictionaryType DWIMetaDataDictionaryValidator::GetMetaDataDictionary()
{
  return m_dict;
}

void DWIMetaDataDictionaryValidator::SetMetaDataDictionary(DWIMetaDataDictionaryValidator::ConstMetaDataDictionaryType rhs)
{
  m_dict = rhs;
}

DWIMetaDataDictionaryValidator::RotationMatrixType DWIMetaDataDictionaryValidator::GetMeasurementFrame() const
{
  RotationMatrixType myRM;
  std::vector<std::vector<double> > retval;
  itk::ExposeMetaData<std::vector<std::vector<double> > >(m_dict, "NRRD_measurement frame", retval);
  if (retval.size() != 0)
    {
      for (size_t r = 0; r < myRM.RowDimensions; ++r)
      {
        for(size_t c=0; c < myRM.ColumnDimensions; ++c)
        {
          myRM(r, c) = retval[r][c];
        }
      }
    }
  else
    {
    itkGenericExceptionMacro("Measurement frame not found in metadata");
    }
  return myRM;
}

void DWIMetaDataDictionaryValidator::SetMeasurementFrame(const RotationMatrixType & input)
{
  std::string       key = "NRRD_measurement frame";
  typedef std::vector<std::vector<double> > MFTYPE;
  MFTYPE dictInput;
  dictInput.resize(input.RowDimensions);
  for(size_t r=0; r < input.RowDimensions; ++r)
  {
    dictInput[r].resize(input.ColumnDimensions);
    for(size_t c=0; c < input.ColumnDimensions; ++c)
    {
      dictInput[r][c] = input(r,c);
    }
  }
itk::EncapsulateMetaData<MFTYPE>(m_dict, key, dictInput);
}

DWIMetaDataDictionaryValidator::Double3x1ArrayType DWIMetaDataDictionaryValidator::GetGradient(int index) const
{
  DWIMetaDataDictionaryValidator::Double3x1ArrayType   currentGradient;
  std::string key = DWIMetaDataDictionaryValidator::GetGradientKeyString(index);
  if (m_dict.HasKey(key))
   {
   std::string   NrrdValue;
   itk::ExposeMetaData<std::string>(m_dict, key, NrrdValue);
   sscanf(NrrdValue.c_str(), "%lf %lf %lf", &currentGradient[0], &currentGradient[1], &currentGradient[2]);
   }
  else
   {
   itkGenericExceptionMacro("Gradient not found in metadata");
   }
  return currentGradient;
}

void DWIMetaDataDictionaryValidator::SetGradient(int index, DWIMetaDataDictionaryValidator::Double3x1ArrayType & gradient)
{
  std::string       key = DWIMetaDataDictionaryValidator::GetGradientKeyString(index);
  char              tmp[64];
  sprintf(tmp, "%lf %lf %lf", gradient[0], gradient[1], gradient[2]);
  std::string       value(tmp);
  DWIMetaDataDictionaryValidator::SetStringDictObject(key, value);
}

int DWIMetaDataDictionaryValidator::GetGradientCount()
{
  DWIMetaDataDictionaryValidator::StringVectorType keys = m_dict.GetKeys();
  int count = 0;
  for( DWIMetaDataDictionaryValidator::StringVectorType::iterator it = keys.begin(); it != keys.end(); ++it )
     {
     if((*it).find("DWMRI_gradient_") != std::string::npos)
       {
       count++;
       }
     // Check for repeated components
     else if((*it).find("DWMRI_NEX") != std::string::npos)
       {
       int repeats=0;
       if( itk::ExposeMetaData<int>(m_dict, "DWMRI_NEX", repeats) )
         {
         count += repeats;
         }
       }
     }
  return count;
}

DWIMetaDataDictionaryValidator::GradientTableType DWIMetaDataDictionaryValidator::GetGradientTable() const
{
  DWIMetaDataDictionaryValidator::GradientTableType myGradientTable;
  DWIMetaDataDictionaryValidator::StringVectorType allKeys = m_dict.GetKeys();
  int count = 0;
  for (DWIMetaDataDictionaryValidator::StringVectorType::iterator it = allKeys.begin();
      it != allKeys.end();
      ++it)
   {
   if((*it).find("DWMRI_gradient_") != std::string::npos)
     {
     std::string   gradientString = (*it).substr(15,4);
     int index = std::atoi(gradientString.c_str());
     DWIMetaDataDictionaryValidator::Double3x1ArrayType curGradientDirection = GetGradient(index);
     myGradientTable.push_back(curGradientDirection);
     }
   count++;
   }
  return myGradientTable;
}

void DWIMetaDataDictionaryValidator::SetGradientTable(GradientTableType & myGradientTable)
{
  int count = 0;
  for( DWIMetaDataDictionaryValidator::GradientTableType::iterator it = myGradientTable.begin();
       it != myGradientTable.end(); ++it )
     {
     DWIMetaDataDictionaryValidator::SetGradient(count++, *it);
     }
  // Remove additional gradients
  std::string nextKey = DWIMetaDataDictionaryValidator::GetGradientKeyString(count++);
  while(m_dict.HasKey(nextKey))
     {
     m_dict.Erase(nextKey);
     nextKey = DWIMetaDataDictionaryValidator::GetGradientKeyString(count++);
     }
}

void DWIMetaDataDictionaryValidator::DeleteGradientTable()
{
  for (itk::MetaDataDictionary::MetaDataDictionaryMapType::iterator it = m_dict.Begin();
       it != m_dict.End(); ++it)
      {
      //http://stackoverflow.com/questions/8234779/how-to-remove-from-a-map-while-iterating-it
      if( (*it).first.find("DWMRI_gradient_") != std::string::npos )
        {
        m_dict.Erase((*it).first);
        }
      }
}

double DWIMetaDataDictionaryValidator::GetBValue() const
{
  double      retval = 0.0;
  std::string valstr;

  itk::ExposeMetaData<std::string>(m_dict, "DWMRI_b-value", valstr);
  std::stringstream ss(valstr);
  ss >> retval;
  return retval;
}

void DWIMetaDataDictionaryValidator::SetBValue(const double bvalue)
{
  const std::string        key = "DWMRI_b-value";
  itk::NumberToString<double> doubleConvert;
  const std::string        valstr(doubleConvert(bvalue));
  DWIMetaDataDictionaryValidator::SetStringDictObject(key, valstr);
}

std::string DWIMetaDataDictionaryValidator::GetIndexedKeyString(const std::string base_key_name, const size_t index) const
{
  std::ostringstream itos;
  itos << index;
  return base_key_name+"["+itos.str()+"]";
}

void DWIMetaDataDictionaryValidator::GenericSetStringVector(const std::vector<std::string> & values, const std::string & KeyBaseName)
{
  for(size_t index=0; index< values.size(); ++index)
    {
    const std::string currKey = this->GetIndexedKeyString(KeyBaseName,index);
    /*
     * Only those axises that have information associated with it need to be set here.
     * itkNrrdIO uses "???" for the list/vector axis.
     * Also, it handles the permutation properly based on the "kinds" field in image.
     */
    if( values[index] != "???" )
      {
      itk::EncapsulateMetaData< std::string >(this->m_dict,
                                              currKey,
                                              values[index]);
      }
    }
}

void DWIMetaDataDictionaryValidator::GenericSetDoubleVector(const std::vector<double> & values, const std::string & KeyBaseName)
{
  for(size_t index=0; index< values.size(); ++index)
    {
    const std::string currKey = this->GetIndexedKeyString(KeyBaseName,index);
    /*
     * thickness is a per axis value, and the only axis that has information associated with it needs to be set here.
     * itkNrrdIO uses "nan" for the other axis.
     * Also, itkNrrdIO handles the correct permutation based on the "kinds" field in image. It sets the thickness value
     * for the 3rd space/domain.
     */
    if( !vnl_math_isnan(values[index]) )
      {
      itk::EncapsulateMetaData< double >(this->m_dict,
                                         currKey,
                                         values[index]);
      }
    }
}


std::vector<std::string> DWIMetaDataDictionaryValidator::GenericGetStringVector(const std::string & KeyBaseName,
                                                                                const size_t  numElements,
                                                                                const std::string defaultValue) const
{
  std::vector<std::string> values(numElements);
  for(size_t index=0; index< values.size() ; ++index)
    {
    const std::string currKey = this->GetIndexedKeyString(KeyBaseName,index);
    std::string temp("0.0");
    if (itk::ExposeMetaData(this->m_dict, currKey, temp ) )
      {
      values[index]=temp;
      }
    else
      {
      values[index]=defaultValue;
      }
    }
  return values;
}

std::vector<double> DWIMetaDataDictionaryValidator::GenericGetDoubleVector(const std::string & KeyBaseName,
                                                                           const size_t  numElements,
                                                                           const double defaultValue) const
{
  std::vector<double> values(numElements);
  for(size_t index=0; index< values.size() ; ++index)
    {
    const std::string currKey = this->GetIndexedKeyString(KeyBaseName,index);
    double temp=0.0;
    if (itk::ExposeMetaData(this->m_dict, currKey, temp ) )
      {
      values[index]=temp;
      }
    else
      {
      values[index]=defaultValue;
      }
    }
  return values;
}

std::vector<std::string>  DWIMetaDataDictionaryValidator::GetCenterings() const
{
  // centerings units are always 4  for DWI
  const std::string KeyBaseName("NRRD_centerings");
  return this->GenericGetStringVector(KeyBaseName,4,"???");
}


void DWIMetaDataDictionaryValidator::SetCenterings(const std::vector<std::string> & values)
{
  const std::string _cellString("cell");
  const std::string _dkString("???");

  std::vector<std::string> outValues=values;
  while(outValues.size() < 4)  //PADD to 4D for diffusion last element is list
    {
    outValues.push_back(_dkString);
    }

  for(size_t i=0; i< outValues.size(); ++i)
    {
    const std::string x = outValues[i];
    if ( (x != _cellString ) && (x != _dkString) )
      {
      std::cout << "ERROR: " << i << " " << x << " Not a valid NRRD_centerings" << std::endl;
      }
    }
  const std::string KeyBaseName("NRRD_centerings");
  this->GenericSetStringVector(outValues,KeyBaseName);
}

std::vector<double> DWIMetaDataDictionaryValidator::GetThicknesses() const
{
  // Thickness units are always 4 for DWI
  const std::string KeyBaseName("NRRD_thicknesses");
  return this->GenericGetDoubleVector(KeyBaseName,4,std::numeric_limits<double>::quiet_NaN());
}

void DWIMetaDataDictionaryValidator::SetThicknesses(const std::vector<double> & values)
{
  const std::string KeyBaseName("NRRD_thicknesses");
  this->GenericSetDoubleVector(values,KeyBaseName);
}


std::string DWIMetaDataDictionaryValidator::GetModality() const
{
  const std::string KeyBaseName("modality");
  std::string valstr;
  itk::ExposeMetaData<std::string>(m_dict, KeyBaseName, valstr);
  return valstr;
}

void DWIMetaDataDictionaryValidator::SetModality(const std::string & value)
{
  const std::string KeyBaseName("modality");
  if ( value != "DWMRI" )
    {
    std::cout << "ERROR: " << value << " Not a valid modality" << std::endl;
    }
  this->SetStringDictObject(KeyBaseName, value);
}
