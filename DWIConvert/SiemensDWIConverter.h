/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
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
 *
 *=========================================================================*/
#ifndef __SiemensDWIConverter_h
#define __SiemensDWIConverter_h
#include "DWIConverter.h"
#include "StringContains.h"
#include "DWIDICOMConverterBase.h"

/** specific converter for Siemens scanners*/
class SiemensDWIConverter : public DWIDICOMConverterBase
{
public:
  SiemensDWIConverter(DCMTKFileVector &                  allHeaders,
                      DWIConverter::FileNamesContainer & inputFileNames,
                      const bool                         useBMatrixGradientDirections,
                      const double                       smallGradientThreshold);
  ~SiemensDWIConverter() override;

  template <typename T>
  T
  CSAExtractFromString(const char * ptr);

  class CSAItem : public std::vector<std::string>
  {
  public:
    using SuperClass = std::vector<std::string>;

    itk::uint32_t vm;
    std::string   vr;

    CSAItem(unsigned int length)
      : SuperClass(length)
      , vm(0)
    {}
    CSAItem()
      : SuperClass()
    {}
    CSAItem(const CSAItem & other)
      : SuperClass(other.size())
    {
      *this = other;
    }
    CSAItem &
    operator=(const CSAItem & other)
    {
      this->resize(0);
      for (auto it = other.begin(); it != other.end(); ++it)
      {
        this->push_back(*it);
      }
      this->vm = other.vm;
      this->vr = other.vr;
      return *this;
    }

    template <typename T>
    std::vector<T>
    AsVector() const
    {
      std::vector<T> rval;
      for (unsigned i = 0; i < this->size(); ++i)
      {
        if (!(*this)[i].empty())
        {
          T                 val = 0;
          std::stringstream convert((*this)[i]);
          convert >> val;
          rval.push_back(val);
        }
      }
      return rval;
    }
    void
    DebugPrint() const
    {
      std::cerr << "  VM = " << this->vm << " VR = " << this->vr << std::endl << "    ";
      bool firstTime(false);
      for (auto it = this->begin(); it != this->end(); ++it)
      {
        if (firstTime)
        {
          firstTime = false;
        }
        else
        {
          std::cerr << " ";
        }
        std::cerr << *it;
      }
      std::cerr << std::endl;
    }
  };

  class CSAHeader : public std::map<std::string, CSAItem>
  {
  public:
    void
    DebugPrint() const
    {
      for (auto it = this->begin(); it != this->end(); ++it)
      {
        std::cerr << it->first << std::endl;
        it->second.DebugPrint();
      }
    }
  };


  void
  DecodeCSAHeader(CSAHeader & header, const std::string & infoString);


  /** Siemens datasets are either in the
   *  normal sequential volume arrangement or
   *  in mosaic datasets, where each slice contains
   *  a collection of 2D slices arranged in a single
   *  mosaic slice.
   */
  void
  LoadDicomDirectory() override;

  double
  ExtractBValue(CSAHeader * csaHeader, unsigned int strideVolume);

  bool
  ExtractGradientDirection(CSAHeader * csaHeader, unsigned int strideVolume, vnl_vector_fixed<double, 3> & gradient);

  bool
  ExtractBMatrix(CSAHeader * csaHeader, unsigned int strideVolume, vnl_matrix_fixed<double, 3, 3> & bMatrix);
  /**
   * @brief  find the bvalues and gradient vectors
   */
  void
  ExtractDWIData() override;

private:
  static bool
  IsZeroMag(DWIMetaDataDictionaryValidator::GradientDirectionType vec);

protected:
  /** turn a mosaic image back into a sequential volume image */
  void
  DeMosaic();

  unsigned int
  ConvertFromCharPtr(const char * s);
  /** pull data out of Siemens scans.
   *
   *  Siemens sticks most of the DTI information into a single
   *  OB-format entry.  This is actually rigidly structured, but
   *  this function depends on the needed data living at fixed offset
   *  from the beginning of the name of each tag, and ignores the
   *  internal structure documented in the Siemens Dicom Compliance
   *  document.
   */
  unsigned int
  ExtractSiemensDiffusionInformation(const std::string &   tagString,
                                     const std::string &   nameString,
                                     std::vector<double> & valueArray);

  void
  CheckCSAHeaderAvailable();

  void
  AddFlagsToDictionary() override;

private:
  double       m_SmallGradientThreshold;
  unsigned int m_MMosaic;
  unsigned int m_NMosaic;
  unsigned int m_Stride;
  bool         m_HasCSAHeader;
};
#include "SiemensDWIConverter.hxx"

#endif // __SiemensDWIConverter_h
