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
 * TODO: Determine if extensions may be made to facilitate writing
 * nifti files with associated bvec and bval files
 *
 */
#ifndef DWIMetaDataDictionaryValidator_h_
#define DWIMetaDataDictionaryValidator_h_

#include <string>
#include <vector>
#include "NrrdIO.h"
#include "itkMetaDataObject.h"
#include "itkMetaDataDictionary.h"

/*
 * Note that the following fields do not need to be set by Metadata validator,
   since itkNrrdImageIO extract them from image information:
   - type e.g. short
   - dimension i.e. 4
   - space i.e. left-posterior-superior
   - sizes
   - space directions
   - kinds
   - endian
   - encoding
   - space origin

 * Also, "space units" is an optional field, and ITK does not set that
   in the output MetaDataDictionary even if it exists in the header of
   read Nrrd file. Therefore, this field does not need to be set by the
   validator.

 * The following fields should be managed by validator:
    - thickness
    - centerings
    - measurement frame
    - modality
    - b-value
    - gradient_xxxx
 */


// This class is wrapper for std::vector to make sure it is defined with a fixed size.
// An easy alternative was using std::array<Type, Size>
// However, std::array is only available in c++ 11. We need to build with c++03, since Slicer does not complie with c++11.
template<class Type, int MySize>
class MyArrayWrapper: public std::vector<Type>
{
public:
  MyArrayWrapper(): std::vector<Type>(MySize) {}

private:
  void resize(size_t) ITK_DELETED_FUNCTION;
  void resize(size_t, Type) ITK_DELETED_FUNCTION;
};

class DWIMetaDataDictionaryValidator
{
 private:
  itk::MetaDataDictionary m_dict;
  void SetStringDictObject(const std::string, const std::string);

 protected:
  std::string GetGradientKeyString(int) const;
  std::string GetIndexedKeyString(const std::string base_key_name, const size_t index) const;

  void GenericSetStringVector(const std::vector<std::string> & values, const std::string & KeyBaseName);
  void GenericSetDoubleVector(const std::vector<double> & values, const std::string & KeyBaseName);

  std::vector<std::string> GenericGetStringVector(const std::string & KeyBaseName,
                                                  const size_t numElemnents,
                                                  const std::string defaultValue) const;

  std::vector<double> GenericGetDoubleVector(const std::string & KeyBaseName,
                                             const size_t  numElements,
                                             const double defaultValue) const;
 public:
  // 3D
  typedef MyArrayWrapper<int, 3>                   Integer3x1ArrayType;
  typedef MyArrayWrapper<double, 3>                Double3x1ArrayType;
  typedef MyArrayWrapper<std::string, 3>           String3x1ArrayType;
  // 4D
  typedef MyArrayWrapper<int, 4>                   Integer4x1ArrayType;
  typedef MyArrayWrapper<double, 4>                Double4ArrayType;
  typedef MyArrayWrapper<std::string, 4>           String4x1ArrayType;

  typedef std::vector<double>                      DoubleVectorType;
  typedef std::vector<std::string>                 StringVectorType;
  typedef std::vector<MyArrayWrapper<double, 3> >  GradientTableType;
  typedef std::vector<std::vector<double> >        MeasurementFrameType;
  typedef MyArrayWrapper<DoubleVectorType, 3>      SpaceDirectionType;

  typedef itk::MetaDataDictionary &            MetaDataDictionaryType;
  typedef const itk::MetaDataDictionary &      ConstMetaDataDictionaryType;

  DWIMetaDataDictionaryValidator();
  ~DWIMetaDataDictionaryValidator();
  // metadata dictionary methods
  MetaDataDictionaryType GetMetaDataDictionary();
  void SetMetaDataDictionary( ConstMetaDataDictionaryType );

  // measurement frame
  std::vector<std::vector<double> > GetMeasurementFrame() const;
  void SetMeasurementFrame(const std::vector<std::vector<double> > & );

  // gradients
  Double3x1ArrayType GetGradient(int) const;
  int GetGradientCount();
  void SetGradient(int, Double3x1ArrayType & );
  GradientTableType GetGradientTable() const;
  void SetGradientTable(GradientTableType & );
  void DeleteGradientTable();

// GetBValue is defined as a macro on windows. Need to undef it if present.
// The macro is usually defined in WinGDI.h itself included by Windows.h.
#ifdef GetBValue
#undef GetBValue
#endif

  // b-value
  double GetBValue() const;
  void SetBValue(const double);

  // centering
  std::vector<std::string> GetCenterings() const;
  void SetCenterings(const std::vector<std::string> & );

  //thickness
  std::vector<double> GetThicknesses() const;
  void SetThicknesses(const std::vector<double> & values);

  //modality
  std::string GetModality() const;
  void SetModality(const std::string & value);
};

#endif // DWIMetaDataDictionaryValidator_h_
