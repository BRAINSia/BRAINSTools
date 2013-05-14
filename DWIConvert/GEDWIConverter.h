#ifndef __GEDWIConverter_h
#define __GEDWIConverter_h
#include "DWIConverter.h"

/** Specific converter for GE Scanners */
class GEDWIConverter : public DWIConverter
{
public:
  GEDWIConverter(DWIConverter::DCMTKFileVector &allHeaders,
                 DWIConverter::FileNamesContainer &inputFileNames,
                 bool useBMatrixGradientDirections) : DWIConverter(allHeaders,inputFileNames,
                                                                   useBMatrixGradientDirections)
    {
    }
  virtual ~GEDWIConverter() {}
  virtual void LoadDicomDirectory()
    {
      this->DWIConverter::LoadDicomDirectory();
      this->m_MeasurementFrame = this->m_LPSDirCos;
      this->DetermineSliceOrderIS();
      this->SetDirectionsFromSliceOrder();
      this->m_NVolume = this->m_NSlice / this->m_SlicesPerVolume;
    }

  void ExtractDWIData()
    {
      std::string ModelName;
      // OK, so there is an accomdation made on the basis of one site
      // having garbage BVal/GVectors.  It has to do with variations
      // of behavior of the Signa HDxt scanner.
      // In all cases, the data is thus:
      // BVal = [0043,1039]
      // GVec[0] = [0019.10bb] GVec[1] = [0019,10bc] GVec[2] = [0019,10bd]
      // there are 3 possible encodings of this data for GE scanners:
      // 1. As IS/DS -- integer an decimal strings -- the normal
      // behavior
      // 2. As OB, but the byte data is binary and may need byte
      // swapping.
      // 3. As OB, but it's actuall as case 1 -- numbers represented
      // as strings.
      // I'm accounting for these cases by looking specifically for
      // the Signa HDxt scanner, and if it doesn't find IS/DS data,
      // look for char strings in the OB data.
      // Honestly this is not an optimal way to handle this
      // situation. In an ideal world we'd have accurate knowledge of
      // what each Scanner/Software Version is doing in these tags,
      // and handle them accordingly. But we don't live in that world.
      bool isSignaHDxt(false);
      if( this->m_Headers[0]->GetElementLO(0x0008, 0x001090, ModelName, false) == EXIT_SUCCESS &&
          (ModelName == "Signa HDxt" || ModelName == "SIGNA HDx"
           || ModelName == "GENESIS_SIGNA") )
        {
        isSignaHDxt = true;
        }

      // assume volume interleaving
      std::cout << "Number of Slices: " << this->m_NSlice << std::endl;
      std::cout << "Number of Volume: " << this->m_NVolume << std::endl;
      std::cout << "Number of Slices in each volume: " << this->m_SlicesPerVolume << std::endl;
      for( unsigned int k = 0; k < this->m_NSlice; k += this->m_SlicesPerVolume )
        {
        // parsing bvalue and gradient directions
        DWIConverter::DiffusionVectorType vect3d;
        vect3d.fill( 0 );
        // for some weird reason this item in the GE dicom
        // header is stored as an IS (Integer String) element.
        ::itk::int32_t intb;
        if( !isSignaHDxt )
          {
          this->m_Headers[k]->GetElementISorOB(0x0043, 0x1039, intb);
          }
        else
          {
          if( this->m_Headers[k]->GetElementIS(0x0043, 0x1039, intb, false) != EXIT_SUCCESS )
            {
            std::string val;
            this->m_Headers[k]->GetElementOB(0x0043, 0x1039, val);
            size_t slashpos = val.find('\\');
            val = val.substr(0, slashpos);
            std::stringstream s(val);
            s >> intb;
            }
          }
        float b = static_cast<float>(intb);
        for( unsigned elementNum = 0x10bb; elementNum <= 0x10bd; ++elementNum )
          {
          int vecI(elementNum - 0x10bb);
          if( !isSignaHDxt )
            {
            this->m_Headers[k]->GetElementDSorOB(0x0019, elementNum, vect3d[vecI]);
            }
          else
            {
            if( this->m_Headers[k]->GetElementDS(0x0019, elementNum, 1, &vect3d[vecI], false) != EXIT_SUCCESS )
              {
              std::string val;
              this->m_Headers[k]->GetElementOB(0x0019, elementNum, val);
              std::stringstream s(val);
              s >> vect3d[vecI];
              }
            }
          }

        vect3d[0] = -vect3d[0];
        vect3d[1] = -vect3d[1];

        this->m_BValues.push_back( b );
        if( b == 0 )
          {
          vect3d.fill( 0 );
          this->m_DiffusionVectors.push_back(vect3d);
          }
        else
          {
          // vect3d.normalize();
          this->m_DiffusionVectors.push_back(vect3d);
          }

        std::cout << "B-value: " << b
                  << "; diffusion direction: "
                  << this->m_DoubleConvert(vect3d[0])
                  << ", "
                  << this->m_DoubleConvert(vect3d[1])
                  << ", "
                  << this->m_DoubleConvert(vect3d[2]) << std::endl;
        }
    }
protected:
  virtual void AddFlagsToDictionary()
    {
      // these have to be dynamically allocated because otherwise there's
      // a malloc error after main exits.

      // relevant GE tags
      DcmDictEntry *GEDictBValue = new DcmDictEntry(0x0043, 0x1039, DcmVR(EVR_IS),
                                                    "B Value of diffusion weighting", 1, 1, 0, true,
                                                    "dicomtonrrd");
      DcmDictEntry *GEDictXGradient = new DcmDictEntry(0x0019, 0x10bb, DcmVR(EVR_DS),
                                                       "X component of gradient direction", 1, 1, 0, true,
                                                       "dicomtonrrd");
      DcmDictEntry *GEDictYGradient = new DcmDictEntry(0x0019, 0x10bc, DcmVR(EVR_DS),
                                                       "Y component of gradient direction", 1, 1, 0, true,
                                                       "dicomtonrrd");
      DcmDictEntry *GEDictZGradient = new DcmDictEntry(0x0019, 0x10bd, DcmVR(EVR_DS),
                                                       "Z component of gradient direction", 1, 1, 0, true,
                                                       "dicomtonrrd");

      itk::DCMTKFileReader::AddDictEntry(GEDictBValue);
      itk::DCMTKFileReader::AddDictEntry(GEDictXGradient);
      itk::DCMTKFileReader::AddDictEntry(GEDictYGradient);
      itk::DCMTKFileReader::AddDictEntry(GEDictZGradient);

    }
};

#endif // __GEDWIConverter_h
