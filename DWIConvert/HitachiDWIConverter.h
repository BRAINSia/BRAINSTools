#ifndef __HitachiDWIConverter_h
#define __HitachiDWIConverter_h
#include "DWIConverter.h"

/** specific converter for Hitachi scanners */
class HitachiDWIConverter : public DWIConverter
{
public:
  HitachiDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                      DWIConverter::FileNamesContainer &inputFileNames,
                      bool useBMatrixGradientDirections) : DWIConverter(allHeaders,inputFileNames,
                                                                        useBMatrixGradientDirections)
    {
    }

  virtual ~HitachiDWIConverter() {}
  /* load dicom directory -- no postprocessing necessary after letting
   * superclass do its thing.
   */
  virtual void LoadDicomDirectory()
    {
      this->DWIConverter::LoadDicomDirectory();
      this->m_SliceOrderIS = false;
      this->SetDirectionsFromSliceOrder();
      this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
    }
  /** extract gradient vectors.
   *  Hitachi apparently supports the Supplement 49 definition
   *  for Diffusion data.-- see page 94 of the Supplement 49 document:
   *  ftp://medical.nema.org/medical/dicom/final/sup49_ft.pdf
   */
  void ExtractDWIData()
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
  virtual void AddFlagsToDictionary()
    {
    }
};

#endif // __HitachiDWIConverter_h
