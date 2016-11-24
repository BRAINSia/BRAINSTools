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

#ifndef __HitachiDWIConverter_h
#define __HitachiDWIConverter_h
#include "DWIConverter.h"

/** specific converter for Hitachi scanners */
class HitachiDWIConverter : public DWIDICOMConverterBase
{
public:
  HitachiDWIConverter(DWIDICOMConverterBase::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections) : DWIDICOMConverterBase(allHeaders,inputFileNames,
                                                                        useBMatrixGradientDirections)
    {
    }

  virtual ~HitachiDWIConverter() {}
  /* load dicom directory -- no postprocessing necessary after letting
   * superclass do its thing.
   */
  virtual void LoadDicomDirectory() ITK_OVERRIDE
    {
      this->m_SliceOrderIS = false;
      this->SetDirectionsFromSliceOrder();
      this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
    }
  /** extract gradient vectors.
   *  Hitachi apparently supports the Supplement 49 definition
   *  for Diffusion data.-- see page 94 of the Supplement 49 document:
   *  ftp://medical.nema.org/medical/dicom/final/sup49_ft.pdf
   */
  void ExtractDWIData() ITK_OVERRIDE
    {
      for(unsigned int k = 0; k < this->m_NSlice; k += this->m_SlicesPerVolume)
        {
        itk::DCMTKSequence SharedFunctionalGroupsSequence;
        this->m_Headers[k]->GetElementSQ(0x5200,0x9229,SharedFunctionalGroupsSequence);
        double b = 0.0;
        SharedFunctionalGroupsSequence.GetElementFD(0x0018,0x9087,b);
        this->m_BValues.push_back(b);
        double doubleArray[3];
        SharedFunctionalGroupsSequence.GetElementFD(0x0018,0x9089,3,doubleArray);
        vnl_vector_fixed<double, 3> vect3d;
        for(unsigned g = 0; g < 3; ++g)
          {
          vect3d[g] = doubleArray[g];
          }
        this->m_DiffusionVectors.push_back(vect3d);
        }
    }
protected:
  virtual void AddFlagsToDictionary() ITK_OVERRIDE
    {
    }
};

#endif // __HitachiDWIConverter_h
