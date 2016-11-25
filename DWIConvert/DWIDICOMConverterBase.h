//
// Created by Johnson, Hans J on 11/24/16.
//
#ifndef BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
#define BRAINSTOOLS_DWIDICOMCONVERTERBASE_H

#include <vector>
#include <iostream>
#include "DWIConverter.h"
#include "itkDCMTKSeriesFileNames.h"
#include "itkMacro.h"
#include "itkDCMTKImageIO.h"
#include "itkImage.h"
#include "itkDCMTKFileReader.h"
#include "itkNumberToString.h"
#include "DWIConvertUtils.h"

class DWIDICOMConverterBase : public DWIConverter {
 public:

  typedef itk::DCMTKSeriesFileNames           InputNamesGeneratorType;
  typedef std::vector<itk::DCMTKFileReader *> DCMTKFileVector;

  DWIDICOMConverterBase(const DCMTKFileVector &allHeaders,
               const FileNamesContainer &inputFileNames,
               const bool useBMatrixGradientDirections,
               const bool FSLFileFormatHorizontalBy3Rows) :
                                                    DWIConverter(inputFileNames, FSLFileFormatHorizontalBy3Rows),
                                                    m_UseBMatrixGradientDirections(useBMatrixGradientDirections),
                                                    m_Headers(allHeaders)

    {

    }

  virtual void LoadFromDisk(){
    this->LoadDicomDirectory();
  }

  virtual void LoadDicomDirectory()
    {
      //
      // add vendor-specific flags to dictionary
      this->AddFlagsToDictionary();
      //
      // load the volume, either single or multivolume.
      m_NSlice = this->m_InputFileNames.size();
      itk::DCMTKImageIO::Pointer dcmtkIO = itk::DCMTKImageIO::New();
      if( this->m_InputFileNames.size() > 1 )
        {
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetImageIO( dcmtkIO );
        reader->SetFileNames( this->m_InputFileNames );
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::__1::cerr << "Exception thrown while reading DICOM volume"
                         << std::__1::endl;
          std::__1::cerr << excp << std::__1::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = false;
        }
      else
        {
        itk::ImageFileReader<VolumeType>::Pointer reader =
          itk::ImageFileReader<VolumeType>::New();
        reader->SetImageIO( dcmtkIO );
        reader->SetFileName( this->m_InputFileNames[0] );
        m_NSlice = this->m_InputFileNames.size();
        try
          {
          reader->Update();
          }
        catch( itk::ExceptionObject & excp )
          {
          std::__1::cerr << "Exception thrown while reading the series" << std::__1::endl;
          std::__1::cerr << excp << std::__1::endl;
          throw;
          }
        m_Volume = reader->GetOutput();
        m_MultiSliceVolume = true;
        }

      // figure out image dimensions
      {
        //TODO:  Remove this redundant code.
        unsigned short rows, cols;
        m_Headers[0]->GetElementUS(0x0028, 0x0010, rows);
        m_Headers[0]->GetElementUS(0x0028, 0x0011, cols);


        if(cols != this->GetCols() )
        {
          itkGenericExceptionMacro(<< "ERROR:  Cols do not match what was read by image " << cols <<  " != " <<
            this->m_Volume->GetLargestPossibleRegion().GetSize()[1] << std::endl
          )
        }

        if(rows != this->GetRows() )
        {
          itkGenericExceptionMacro(<< "ERROR:  Rows do not match what was read by image " << rows << " != " <<
            this->m_Volume->GetLargestPossibleRegion().GetSize()[0] << std::endl
          )
        }
      }

      {
      // origin
      double origin[3];
      m_Headers[0]->GetOrigin(origin);
      VolumeType::PointType imOrigin;
      imOrigin[0] = origin[0];
      imOrigin[1] = origin[1];
      imOrigin[2] = origin[2];
      this->m_Volume->SetOrigin(imOrigin);
      }
      // spacing
      {

      double spacing[3];
      m_Headers[0]->GetSpacing(spacing);
      SpacingType imSpacing;
      imSpacing[0] = spacing[0];
      imSpacing[1] = spacing[1];
      imSpacing[2] = spacing[2];
      m_Volume->SetSpacing(imSpacing);
      }

      // a map of ints keyed by the slice location string
      // reported in the dicom file.  The number of slices per
      // volume is the same as the number of unique slice locations
      std::__1::map<std::__1::string, int> sliceLocations;
      //
      // check for interleave
      if( !this->m_MultiSliceVolume )
        {
        // Make a hash of the sliceLocations in order to get the correct
        // count.  This is more reliable since SliceLocation may not be available.
        std::__1::vector<int>         sliceLocationIndicator;
        std::__1::vector<std::__1::string> sliceLocationStrings;

        sliceLocationIndicator.resize( this->m_NSlice );
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          std::__1::string originString;
          this->m_Headers[k]->GetElementDS(0x0020, 0x0032, originString );
          sliceLocationStrings.push_back( originString );
          sliceLocations[originString]++;
          }

        // this seems like a crazy way to figure out if slices are
        // interleaved, but it works. Perhaps replace with comparing
        // the reported location between the first two slices?
        // Would be less clever-looking and devious, but would require
        // less computation.
        for( unsigned int k = 0; k < this->m_NSlice; ++k )
          {
          std::map<std::string, int>::iterator it = sliceLocations.find( sliceLocationStrings[k] );
          sliceLocationIndicator[k] = distance( sliceLocations.begin(), it );
          }

        // sanity check on # of volumes versus # of dicom files
        if(this->m_Headers.size() % sliceLocations.size() != 0)
          {
          itkGenericExceptionMacro(<< "Missing DICOM Slice files: Number of slice files ("
                            << this->m_Headers.size() << ") not evenly divisible by"
                            << " the number of slice locations ");
          }

        this->m_SlicesPerVolume = sliceLocations.size();
        std::__1::cout << "=================== this->m_SlicesPerVolume:" << this->m_SlicesPerVolume << std::__1::endl;


        // if the this->m_SlicesPerVolume == 1, de-interleaving won't do
        // anything so there's no point in doing it.
        if( this->m_NSlice >= 2 && this->m_SlicesPerVolume > 1 )
          {
          if( sliceLocationIndicator[0] != sliceLocationIndicator[1] )
            {
            std::__1::cout << "Dicom images are ordered in a volume interleaving way." << std::__1::endl;
            }
          else
            {
            std::__1::cout << "Dicom images are ordered in a slice interleaving way." << std::__1::endl;
            this->m_IsInterleaved = true;
            // reorder slices into a volume interleaving manner
            DeInterleaveVolume();
            }
          }
        }

    {
    VolumeType::DirectionType LPSDirCos;
    LPSDirCos.SetIdentity();

    // check ImageOrientationPatient and figure out slice direction in
    // L-P-I (right-handed) system.
    // In Dicom, the coordinate frame is L-P by default. Look at
    // http://medical.nema.org/dicom/2007/07_03pu.pdf ,  page 301
    double dirCosArray[6];
    // 0020,0037 -- Image Orientation (Patient)
    this->m_Headers[0]->GetDirCosArray(dirCosArray);
    double *dirCosArrayP = dirCosArray;
    for( unsigned i = 0; i < 2; ++i )
      {
      for( unsigned j = 0; j < 3; ++j, ++dirCosArrayP )
        {
        LPSDirCos[j][i] = *dirCosArrayP;
        }
      }

    // Cross product, this gives I-axis direction
    LPSDirCos[0][2] = LPSDirCos[1][0] * LPSDirCos[2][1] - LPSDirCos[2][0] * LPSDirCos[1][1];
    LPSDirCos[1][2] = LPSDirCos[2][0] * LPSDirCos[0][1] - LPSDirCos[0][0] * LPSDirCos[2][1];
    LPSDirCos[2][2] = LPSDirCos[0][0] * LPSDirCos[1][1] - LPSDirCos[1][0] * LPSDirCos[0][1];

    this->m_Volume->SetDirection(LPSDirCos);
    }
    std::__1::cout << "ImageOrientationPatient (0020:0037): ";
    std::__1::cout << "LPS Orientation Matrix" << std::__1::endl;
    std::__1::cout << this->m_Volume->GetDirection() << std::__1::endl;


    std::__1::cout << "this->m_SpacingMatrix" << std::__1::endl;
    std::__1::cout << this->GetSpacingMatrix() << std::__1::endl;

    std::__1::cout << "NRRDSpaceDirection" << std::__1::endl;
    std::__1::cout << this->GetNRRDSpaceDirection() << std::__1::endl;

    }

protected:
  /* determine if slice order is inferior to superior */
  void DetermineSliceOrderIS()
    {
      double image0Origin[3];
      image0Origin[0]=this->m_Volume->GetOrigin()[0];
      image0Origin[1]=this->m_Volume->GetOrigin()[1];
      image0Origin[2]=this->m_Volume->GetOrigin()[2];
      std::__1::cout << "Slice 0: " << image0Origin[0] << " "
                     << image0Origin[1] << " " << image0Origin[2] << std::__1::endl;

      // assume volume interleaving, i.e. the second dicom file stores
      // the second slice in the same volume as the first dicom file
      double image1Origin[3];

      unsigned long nextSlice = 0;
      if (this->m_Headers.size() > 1)
        {
        // assuming multiple files is invalid for single-file volume: http://www.na-mic.org/Bug/view.php?id=4105
        nextSlice = this->m_IsInterleaved ? this->m_NVolume : 1;
        }

      this->m_Headers[nextSlice]->GetOrigin(image1Origin);
      std::__1::cout << "Slice " << nextSlice << ": " << image1Origin[0] << " " << image1Origin[1]
                     << " " << image1Origin[2] << std::__1::endl;

      image1Origin[0] -= image0Origin[0];
      image1Origin[1] -= image0Origin[1];
      image1Origin[2] -= image0Origin[2];
      const RotationMatrixType & NRRDSpaceDirection = this->GetNRRDSpaceDirection();
      double x1 = image1Origin[0] * (NRRDSpaceDirection[0][2])
        + image1Origin[1] * (NRRDSpaceDirection[1][2])
        + image1Origin[2] * (NRRDSpaceDirection[2][2]);
      if( x1 < 0 )
        {
        this->m_SliceOrderIS = false;
        }
    }

  /** force use of the BMatrix to compute gradients in Siemens data instead of
   *  the reported gradients. which are in many cases bogus.
   */
  const bool m_UseBMatrixGradientDirections;
  /** one file reader per DICOM file in dataset */
  const DCMTKFileVector     m_Headers;
};

#endif //BRAINSTOOLS_DWIDICOMCONVERTERBASE_H
