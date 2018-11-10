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
/*
 * Author: Hans J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef landmarksConstellationModelIO_h
#define landmarksConstellationModelIO_h

#include "landmarksConstellationCommon.h"
#include "landmarksConstellationTrainingDefinitionIO.h"
#include "landmarksConstellationModelBase.h"

#include "itkByteSwapper.h"
#include "itkIO.h"
#include "itkNumberToString.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <map>

#include "BRAINSConstellationDetectorVersion.h"

////
inline void
defineTemplateIndexLocations
  (const float r,
  const float h,
  std::vector<SImageType::PointType::VectorType> & indexLocations)
{
  using CoordType = SImageType::PointType::VectorType::ComponentType;

  // Reserve space that will be needed
  indexLocations.reserve( static_cast<unsigned int>( 4 * r * r * h ) );

  const CoordType h_2 = h / 2;
  const CoordType r2 = r * r;
  for( CoordType SI = -r; SI <= r; SI += 1.0F )
    {
    for( CoordType PA = -r; PA <= r; PA += 1.0F )
      {
      for( CoordType LR = -h_2; LR <= h_2; LR += 1.0F )
        {
        if( ( SI * SI + PA * PA )  <= r2 )  // a suspicious place
          {
          SImageType::PointType::VectorType temp;
          temp[0] = LR;
          temp[1] = PA;
          temp[2] = SI;
          indexLocations.push_back(temp);
          }
        }
      }
    }
}

class landmarksConstellationModelIO : public landmarksConstellationModelBase
{
private:
  using Self = landmarksConstellationModelIO;
  using ioErr = enum { readFail, writeFail };
  typedef enum { file_signature = 0x12345678,
                 swapped_file_signature = 0x78563412 } fileSig;
public:
  using IndexLocationVectorType = std::vector<SImageType::PointType::VectorType>;
  using FloatVectorType = std::vector<float>;
  using Float2DVectorType = std::vector<FloatVectorType>;
  using Float3DVectorType = std::vector<Float2DVectorType>;

  using FloatVectorIterator = FloatVectorType::iterator;
  using ConstFloatVectorIterator = FloatVectorType::const_iterator;

  using Float2DVectorIterator = Float2DVectorType::iterator;
  using ConstFloat2DVectorIterator = Float2DVectorType::const_iterator;

  using Float3DVectorIterator = Float3DVectorType::iterator;
  using ConstFloat3DVectorIterator = Float3DVectorType::const_iterator;
public:

  landmarksConstellationModelIO()
  {
    m_Swapped = false;
    m_RPPC_to_RPAC_angleMean = 0;
    m_RPAC_over_RPPCMean = 0;
  }

  ~landmarksConstellationModelIO() override
  {
  }

  void SetCMtoRPMean(const SImageType::PointType::VectorType & CMtoRPMean)
  {
    this->m_CMtoRPMean = CMtoRPMean;
  }

  void SetRPtoXMean(std::string name, const SImageType::PointType::VectorType & RPtoXMean)
  {
    this->m_RPtoXMean[name] = RPtoXMean;
  }

  void SetRPtoCECMean(const SImageType::PointType::VectorType & RPtoCECMean)
  {
    this->m_RPtoCECMean = RPtoCECMean;
  }

  void SetRPPC_to_RPAC_angleMean(float RPPC_to_RPAC_angleMean)
  {
    this->m_RPPC_to_RPAC_angleMean = RPPC_to_RPAC_angleMean;
  }

  void SetRPAC_over_RPPCMean(float RPAC_over_RPPCMean)
  {
    this->m_RPAC_over_RPPCMean = RPAC_over_RPPCMean;
  }

  const SImageType::PointType::VectorType & GetCMtoRPMean() const
  {
    return this->m_CMtoRPMean;
  }

  const SImageType::PointType::VectorType & GetRPtoXMean(std::string name)
  {
    return this->m_RPtoXMean[name];
  }

  const SImageType::PointType::VectorType & GetRPtoCECMean() const
  {
    return this->m_RPtoCECMean;
  }

  float GetRPPC_to_RPAC_angleMean() const
  {
    return this->m_RPPC_to_RPAC_angleMean;
  }

  float GetRPAC_over_RPPCMean() const
  {
    return this->m_RPAC_over_RPPCMean;
  }

  /**
   * Creates a mean for each of the angles by collapsing the image
   * dimensions.
   *
   * @author hjohnson (9/6/2008)
   *
   * @param output
   * @param input
   */
  void ComputeAllMeans(Float2DVectorType & output, const Float3DVectorType & input)
  {
    // First allocate the output mememory
    output.resize( input[0].size() );
    for( unsigned int q = 0; q < output.size(); q++ )
      {
      output[q].resize( input[0][0].size() );
      for( FloatVectorIterator oit = output[q].begin(); oit != output[q].end(); ++oit )
        {
        *oit = 0;
        }
      }
    for( ConstFloat3DVectorIterator curr_dataset = input.begin(); curr_dataset != input.end(); ++curr_dataset )
      {
      ConstFloat2DVectorIterator input_angleit = curr_dataset->begin();
      Float2DVectorIterator      output_angleit = output.begin();
      while( input_angleit != curr_dataset->end() && output_angleit != output.end() )
        {
        ConstFloatVectorIterator init = input_angleit->begin();
        FloatVectorIterator      outit = output_angleit->begin();
        while( init != input_angleit->end() && outit != output_angleit->end() )
          {
          *outit += *init;
          ++outit;
          ++init;
          }

        ++input_angleit;
        ++output_angleit;
        }
      }
    // Now divide by number of data sets
    const float inv_size = 1.0 / input.size();
    for( unsigned int q = 0; q < output.size(); q++ )
      {
      for( FloatVectorIterator oit = output[q].begin(); oit != output[q].end(); ++oit )
        {
        *oit *= inv_size;
        }
      }
  }

  // Access the internal memory locations for modification.
  FloatVectorType & AccessTemplate(std::string name, const unsigned int indexDataSet, const unsigned int indexAngle)
  {
    if( this->m_Templates.find(name) == this->m_Templates.end() )
      {
      itkGenericExceptionMacro(<< "Attempt to access an undifined landmark template for "
                               << name);
      }
    return this->m_Templates[name][indexDataSet][indexAngle];
  }

  // Access the mean vectors of templates
  const FloatVectorType & AccessTemplateMean(std::string name, const unsigned int indexAngle)
  {
    if( this->m_TemplateMeansComputed.find(name) == this->m_TemplateMeansComputed.end()
        || this->m_Templates.find(name) == this->m_Templates.end() )
      {
      itkGenericExceptionMacro(<< "Attempt to access an undifined landmark template (mean)!");
      }
    ComputeAllMeans(this->m_TemplateMeans[name], this->m_Templates[name]);
    return this->m_TemplateMeans[name][indexAngle];
  }

  const Float2DVectorType & AccessAllTemplateMeans(std::string name)
  {
    if( this->m_TemplateMeansComputed.find(name) == this->m_TemplateMeansComputed.end()
        && this->m_Templates.find(name) == this->m_Templates.end() )
      {
      itkGenericExceptionMacro(<< "Attempt to access an undefined landmark template (mean)!");
      }
    ComputeAllMeans(this->m_TemplateMeans[name], this->m_Templates[name]);
    return this->m_TemplateMeans[name];
  }

  const Float2DVectorType & GetTemplateMeans(const std::string& name)
  {
    return this->m_TemplateMeans[name];
  }

  void WriteModelFile(const std::string & filename)
  {
    //
    //
    // //////////////////////////////////////////////////////////////////////////
    itk::NumberToString<double> doubleToString;
    std::ofstream  output( filename.c_str() ); // open setup file for reading

    if( !output.is_open() )
      {
      std::cerr << "Can't write " << filename << std::endl;
      std::cerr.flush();
      }
    try
      {
      this->Write<unsigned int>(output, file_signature);   // Write out the
                                                           // signature first
      /*
       * WEI: As the new model deals with arbitrary landmarks, we need to
       * write the landmark names to the model file in order to
       * differentiate for example the search radius of each landmark.
       * That means ofstream::write() alone is not sufficient for reading
       * landmark names with variable length.
       */
      this->Write<char>(output, '\n');

      output << BCDVersionString << std::endl;

      std::map<std::string, bool>::const_iterator it2;
      for( it2 = this->m_TemplateMeansComputed.begin();
           it2 != this->m_TemplateMeansComputed.end(); ++it2 )
        {
        output << it2->first.c_str() << std::endl;
        output << doubleToString(this->GetRadius(it2->first) ) << std::endl;
        output << doubleToString(this->GetHeight(it2->first) ) << std::endl;
        }
      output << "END" << std::endl;

      this->Write<unsigned int>( output, this->GetSearchboxDims() );
      this->Write<float>( output, this->GetResolutionUnits() );
      this->Write<unsigned int>( output, this->GetNumDataSets() );
      this->Write<unsigned int>( output, this->GetNumRotationSteps() );

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&output, "before vectors", -1);
#endif
      for( it2 = this->m_TemplateMeansComputed.begin();
           it2 != this->m_TemplateMeansComputed.end(); ++it2 )
        {
        ComputeAllMeans(this->m_TemplateMeans[it2->first],
                        this->m_Templates[it2->first]);
        this->Write(output, this->m_TemplateMeans[it2->first]);
        this->WritedebugMeanImages(it2->first);
        }

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&output, "after vectors", -1);
#endif

      // We only need those RPtoXMeans
      this->Write<double>(output, this->m_RPtoXMean["PC"][0]);
      this->Write<double>(output, this->m_RPtoXMean["PC"][1]);
      this->Write<double>(output, this->m_RPtoXMean["PC"][2]);
      this->Write<double>(output, this->m_CMtoRPMean[0]);
      this->Write<double>(output, this->m_CMtoRPMean[1]);
      this->Write<double>(output, this->m_CMtoRPMean[2]);
      this->Write<double>(output, this->m_RPtoXMean["VN4"][0]);
      this->Write<double>(output, this->m_RPtoXMean["VN4"][1]);
      this->Write<double>(output, this->m_RPtoXMean["VN4"][2]);
      this->Write<double>(output, this->m_RPtoCECMean[0]);
      this->Write<double>(output, this->m_RPtoCECMean[1]);
      this->Write<double>(output, this->m_RPtoCECMean[2]);
      this->Write<double>(output, this->m_RPtoXMean["AC"][0]);
      this->Write<double>(output, this->m_RPtoXMean["AC"][1]);
      this->Write<double>(output, this->m_RPtoXMean["AC"][2]);
      this->Write<float>(output, this->m_RPPC_to_RPAC_angleMean);
      this->Write<float>(output, this->m_RPAC_over_RPPCMean);

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&output, "end of file", -1);
#endif
      }
    catch( ioErr e )
      {
      std::cerr << "Write failed for " << filename << std::endl;
      std::cerr << e << std::endl;
      }
    output.close();
  }

  void PrintHeaderInfo(void)
  {
    //
    //
    // //////////////////////////////////////////////////////////////////////////
    std::cout << "SearchboxDims"    << ": " << this->GetSearchboxDims() << std::endl;
    std::cout << "ResolutionUnits"  << ": " << this->GetResolutionUnits() << std::endl;
    std::cout << "NumDataSets"      << ": " << this->GetNumDataSets() << std::endl;
    std::cout << "NumRotationSteps" << ": " << this->GetNumRotationSteps() << std::endl;

    std::map<std::string, bool>::const_iterator it2;

    for( it2 = this->m_TemplateMeansComputed.begin();
         it2 != this->m_TemplateMeansComputed.end(); ++it2 )
      {
      std::cout << it2->first << "Radius: " << this->GetRadius(it2->first) << std::endl;
      std::cout << it2->first << "Height: " << this->GetHeight(it2->first) << std::endl;
      }

    std::cout << "CMtoRPMean: " << this->m_CMtoRPMean << std::endl;
    std::cout << "RPtoECEMean: " << this->m_RPtoCECMean << std::endl;
    std::cout << "RPtoACMean: " << this->m_RPtoXMean["AC"] << std::endl;
    std::cout << "RPtoPCMean: " << this->m_RPtoXMean["PC"] << std::endl;
    std::cout << "RPtoVN4Mean: " << this->m_RPtoXMean["VN4"] << std::endl;
    std::cout << "RPPC_to_RPAC_angleMean: " << this->m_RPPC_to_RPAC_angleMean << std::endl;
    std::cout << "RPAC_over_RPPCMean: " << this->m_RPAC_over_RPPCMean << std::endl;
  }

  void ReadModelFile(const std::string & filename)
  {
    //
    //
    // //////////////////////////////////////////////////////////////////////////

    std::ifstream input( filename.c_str() ); // open setup file for reading

    if( !input.is_open() )
      {
      itkGenericExceptionMacro(<< "Can't read " << filename);
      }
    try
      {
      unsigned int sig = 0;
      this->Read<unsigned int>(input, sig);
      if( sig != file_signature && sig != swapped_file_signature )
        {
        this->m_Swapped = false;
        }
      else if( sig == swapped_file_signature )
        {
        this->m_Swapped = true;
        }

        {
        std::string Version("INVALID");

        input >> Version;

        std::cout << "Input model file version: " << Version << std::endl;
        if( Version.compare( BCDVersionString ) != 0 )
          {
          itkGenericExceptionMacro(<<"Input model file is outdated.\n"
            << "Input model file version: " << Version
            << ", Required version: " << BCDVersionString << std::endl);
          }
         }
      /*
       * WEI: As the new model deals with arbitrary landmarks, we need to
       * write the landmark names to the model file in order to
       * differentiate for example the search radius of each landmark.
       * That means ofstream::write() alone is not sufficient for reading
       * landmark names with variable length.
       */
      std::string name;
      input >> name;
      while( name.compare("END") != 0 )
        {
        input >> this->m_Radius[name];
        input >> this->m_Height[name];
        this->m_TemplateMeansComputed[name] = false;
        input >> name;
        }

        {
        char tmp = '\0';
        this->Read<char>(input, tmp);
        }


      this->Read<unsigned int>(input, this->m_SearchboxDims);
      this->Read<float>(input, this->m_ResolutionUnits);
      this->Read<unsigned int>(input, this->m_NumDataSets);
      this->Read<unsigned int>(input, this->m_NumRotationSteps);

      std::cout << "NumberOfDataSets: " << this->m_NumDataSets << std::endl;
      std::cout << "SearchBoxDims: " << this->m_SearchboxDims << std::endl;
      std::cout << "ResolutionUnits: " << this->m_ResolutionUnits << std::endl;
      std::cout << "NumberOfRotationSteps: " << this->m_NumRotationSteps << std::endl;

      // initalize the size of m_VectorIndexLocations and m_TemplateMeans
      InitializeModel(false);

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&input, "before vectors", -1);
#endif

        {
        std::map<std::string, bool>::const_iterator it2;
        for( it2 = this->m_TemplateMeansComputed.begin();
             it2 != this->m_TemplateMeansComputed.end(); ++it2 )
          {
          this->Read(input, this->m_TemplateMeans[it2->first]);
          this->m_TemplateMeansComputed[it2->first] = true;
          }
        }

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&input, "after vectors", -1);
#endif

      this->Read<double>(input, this->m_RPtoXMean["PC"][0]);
      this->Read<double>(input, this->m_RPtoXMean["PC"][1]);
      this->Read<double>(input, this->m_RPtoXMean["PC"][2]);
      this->Read<double>(input, this->m_CMtoRPMean[0]);
      this->Read<double>(input, this->m_CMtoRPMean[1]);
      this->Read<double>(input, this->m_CMtoRPMean[2]);
      this->Read<double>(input, this->m_RPtoXMean["VN4"][0]);
      this->Read<double>(input, this->m_RPtoXMean["VN4"][1]);
      this->Read<double>(input, this->m_RPtoXMean["VN4"][2]);
      this->Read<double>(input, this->m_RPtoCECMean[0]);
      this->Read<double>(input, this->m_RPtoCECMean[1]);
      this->Read<double>(input, this->m_RPtoCECMean[2]);
      this->Read<double>(input, this->m_RPtoXMean["AC"][0]);
      this->Read<double>(input, this->m_RPtoXMean["AC"][1]);
      this->Read<double>(input, this->m_RPtoXMean["AC"][2]);
      this->Read<float>(input, this->m_RPPC_to_RPAC_angleMean);
      this->Read<float>(input, this->m_RPAC_over_RPPCMean);

#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&input, "end of file", -1);
#endif
      }
    catch( ioErr e )
      {
#ifdef __USE_OFFSET_DEBUGGING_CODE__
      this->debugOffset(&input, "end of file", -1);
#endif
      std::cerr << "ioErr " << e << std::endl;
      throw;
      }
    input.close();
  }

  class debugImageDescriptor
  {
public:
    Float2DVectorType &     mean;
    IndexLocationVectorType locations;
    float                   r;
    float                   h;
    const char *            debugImageName;

    debugImageDescriptor(Float2DVectorType & _mean, float _r, float _h, const char *_debugImageName) :
      mean(_mean), r(_r), h(_h), debugImageName(_debugImageName)
    {
      defineTemplateIndexLocations(_r, _h, this->locations);
    }
  };

  void WritedebugMeanImages(std::string name)
  {
    debugImageDescriptor DID(this->m_TemplateMeans[name],
                             this->GetRadius(name),
                             this->GetHeight(name),
                             "");

    Float2DVectorType::iterator meanIt = DID.mean.begin();

    for( int j = 0; meanIt != DID.mean.end(); ++meanIt, j++ )
      {
      using FloatImageType = itk::Image<float, 3>;
      using FloatImageRegion = FloatImageType::RegionType;
      using FloatImageSpacing = FloatImageType::SpacingType;
      FloatImageSpacing spacing;
      spacing[0] = spacing[1] = spacing[2] = 1.0;
      FloatImageRegion region;
      region.SetSize
        ( 0, static_cast<FloatImageRegion::SizeValueType>( DID.h + 1 ) );
      region.SetSize
        (1, static_cast<FloatImageRegion::SizeValueType>( 2.0 * DID.r ) + 1);
      region.SetSize
        (2, static_cast<FloatImageRegion::SizeValueType>( 2.0 * DID.r ) + 1);
      FloatImageRegion::IndexType _index;
      _index[0] = _index[1] = _index[2] = 0;
      region.SetIndex(_index);
      FloatImageType::Pointer debugImage = FloatImageType::New();
      debugImage->SetSpacing(spacing);
      // Just use defaults debugImage->SetOrigin(in->GetOrigin());
      // Just use defaults debugImage->SetDirection(in->GetDirection());
      debugImage->SetRegions(region);
      debugImage->Allocate();
      debugImage->FillBuffer(0.0);
      float center[3];

      center[0] = ( DID.h + 1.0 ) * 0.5;
      center[1] = center[2] = ( DID.r + 0.5 );

      IndexLocationVectorType::const_iterator locIt = DID.locations.begin();
      IndexLocationVectorType::const_iterator locItEnd = DID.locations.end();
      for( FloatVectorType::iterator floatIt = ( *meanIt ).begin();
           locIt != locItEnd;
           ++floatIt, ++locIt )
        {
        FloatImageType::IndexType ind;
        for( unsigned k = 0; k < 3; k++ )
          {
          ind[k] = static_cast<FloatImageType::IndexValueType>( center[k] + ( *locIt )[k] );
          }
        debugImage->SetPixel( ind, ( *floatIt ) );
        }
      char buf[512];
      sprintf( buf, "%d%s", j, name.c_str() );
      std::string fname(buf);
      fname += "Mean.nii.gz";
      itkUtil::WriteImage<FloatImageType>(debugImage, fname);
      }
  }

  void InitializeModel(bool CreatingModel)
  {
    //
    //
    // ///////////////////////////////////////////////////////////////////////////////
    // NOTE:  THIS IS NOW AT FULL RESOLUTION, but there was an optimization that
    // had reduced
    //       the resolution in the IS/RP direction to gain some speed
    // efficiency.
    std::map<std::string, bool>::const_iterator it2;

    for( it2 = this->m_TemplateMeansComputed.begin(); it2 != this->m_TemplateMeansComputed.end(); ++it2 )
      {
      defineTemplateIndexLocations(this->GetRadius(it2->first), this->GetHeight(it2->first),
                                   this->m_VectorIndexLocations[it2->first]);
      printf( "%s template size = %u voxels\n", it2->first.c_str(),
              static_cast<unsigned int>( this->m_VectorIndexLocations[it2->first].size() ) );
      if( CreatingModel )
        {
        // Allocate the outter dim for all datasets
        this->m_Templates[it2->first].resize( this->GetNumDataSets() );
        for( unsigned int currdataset = 0; currdataset < this->GetNumDataSets(); ++currdataset )
          {
          // Allocate for number of angles
          this->m_Templates[it2->first][currdataset].resize( this->GetNumRotationSteps() );
          for( unsigned int currrotstep = 0; currrotstep < this->GetNumRotationSteps(); ++currrotstep )
            {
            // Allocate for number of intensity values
            this->m_Templates[it2->first][currdataset][currrotstep].resize
              ( this->m_VectorIndexLocations[it2->first].size() );
            }
          }
        }
      else // make room for reading in template means result
        {
        this->m_TemplateMeans[it2->first].resize( this->GetNumRotationSteps() );
        for( unsigned i = 0; i < this->GetNumRotationSteps(); ++i )
          {
          this->m_TemplateMeans[it2->first][i].resize
            ( this->m_VectorIndexLocations[it2->first].size() );
          }
        }
      }
  }

  void CopyFromModelDefinition(const landmarksConstellationTrainingDefinitionIO & mDef)
  {
    this->SetNumDataSets( mDef.GetNumDataSets() );
    this->SetSearchboxDims( mDef.GetSearchboxDims() );
    this->SetResolutionUnits( mDef.GetResolutionUnits() );
    this->SetInitialRotationAngle( mDef.GetInitialRotationAngle() );
    this->SetInitialRotationStep( mDef.GetInitialRotationStep() );
    this->SetNumRotationSteps( mDef.GetNumRotationSteps() );

      {
      ValMapConstIterator it2;
      for( it2 = mDef.GetRadii().begin(); it2 != mDef.GetRadii().end(); ++it2 )
        {
        this->SetRadius(it2->first, it2->second);
        this->SetHeight( it2->first, mDef.GetHeight(it2->first) );
        this->m_TemplateMeansComputed[it2->first] = false;
        }
      }
  }

  bool NE(double const a, double const b)
  {
    if( ( a < 0.0 && b >= 0.0 )
        || ( b < 0.0 && a >= 0.0 ) )
      {
      return true;
      }
    double absa( fabs(a) ),
    absb( fabs(b) );
    double absdiff( fabs(absa - absb) );
    double avg( ( absa + absb ) / 2.0 );
    if( absdiff > ( avg / 1000.0 ) )
      {
      return true;
      }
    return false;
  }

  bool NE(const std::string & label,
          const Float2DVectorType & a,
          const Float2DVectorType & b)
  {
    itk::NumberToString<double> doubleToString;

    for( unsigned i = 0; i < a.size(); i++ )
      {
      for( unsigned j = 0; j < a[i].size(); j++ )
        {
        if( NE(a[i][j], b[i][j]) )
          {
          std::cerr << label << " a[" << i
                    << "] = " << doubleToString(a[i][j])
                    << "b[" << j << "] = "
                    << doubleToString(b[i][j]) << std::endl;
          return true;
          }
        }
      }
    return false;
  }

  bool operator==(Self & other)
  {
    if( ( NE(this->m_SearchboxDims, other.m_SearchboxDims) )
        || ( NE(this->m_ResolutionUnits, other.m_ResolutionUnits) )
        || ( NE(this->m_NumDataSets, other.m_NumDataSets) )
        || ( NE(this->m_NumRotationSteps, other.m_NumRotationSteps) )
        || ( NE(this->m_RPPC_to_RPAC_angleMean, other.m_RPPC_to_RPAC_angleMean) )
        || ( NE(this->m_RPAC_over_RPPCMean, other.m_RPAC_over_RPPCMean) )
        || ( NE(this->m_CMtoRPMean[0], other.m_CMtoRPMean[0]) )
        || ( NE(this->m_CMtoRPMean[1], other.m_CMtoRPMean[1]) )
        || ( NE(this->m_CMtoRPMean[2], other.m_CMtoRPMean[2]) ) )
      {
      return false;
      }

    std::map<std::string, bool>::const_iterator it2;
    for( it2 = this->m_TemplateMeansComputed.begin(); it2 != this->m_TemplateMeansComputed.end(); ++it2 )
      {
      if( ( NE( this->GetRadius(it2->first), other.GetRadius(it2->first) ) )
          || ( NE( this->GetHeight(it2->first), other.GetHeight(it2->first) ) )
          || ( NE(it2->first + " template mean",
                  this->m_TemplateMeans[it2->first], other.m_TemplateMeans[it2->first]) ) )
        {
        return false;
        }
      }

    if( ( NE(this->m_RPtoXMean["AC"][0], other.m_RPtoXMean["AC"][0]) )
        || ( NE(this->m_RPtoXMean["AC"][1], other.m_RPtoXMean["AC"][1]) )
        || ( NE(this->m_RPtoXMean["PC"][0], other.m_RPtoXMean["PC"][0]) )
        || ( NE(this->m_RPtoXMean["PC"][1], other.m_RPtoXMean["PC"][1]) )
        || ( NE(this->m_RPtoXMean["VN4"][0], other.m_RPtoXMean["VN4"][0]) )
        || ( NE(this->m_RPtoXMean["VN4"][1], other.m_RPtoXMean["VN4"][1]) ) )
      {
      return false;
      }

    return true;
  }

  // HACK:  Kent please wrap these in member variables
  std::map<std::string, landmarksConstellationModelIO::IndexLocationVectorType> m_VectorIndexLocations;
private:
  bool m_Swapped;

  template <typename T>
  void Write(std::ofstream & f, T var)
  {
    if( f.bad() || f.eof() )
      {
      throw landmarksConstellationModelIO::writeFail;
      }
    f.write( reinterpret_cast<char *>( &var ), sizeof( T ) );
  }

  template <typename T>
  void Read(std::ifstream & f, T & var)
  {
    if( f.bad() || f.eof() )
      {
      throw landmarksConstellationModelIO::readFail;
      }

    f.read( reinterpret_cast<char *>( &var ), sizeof( T ) );
    if( this->m_Swapped )
      {
      if( itk::ByteSwapper<T>::SystemIsBigEndian() )
        {
        itk::ByteSwapper<T>::SwapFromSystemToLittleEndian(&var);
        }
      else
        {
        itk::ByteSwapper<T>::SwapFromSystemToBigEndian(&var);
        }
      }
  }

  void Write(std::ofstream & f, const Float2DVectorType & vec)
  {
    for( ConstFloat2DVectorIterator it1 = vec.begin();
         it1 != vec.end(); ++it1 )
      {
      for( ConstFloatVectorIterator it2 = it1->begin();
           it2 != it1->end(); ++it2 )
        {
        this->Write<float>(f, *it2);
        }
      }
  }

  void Read(std::ifstream & f, Float2DVectorType & vec)
  {
    for( Float2DVectorIterator it1 = vec.begin();
         it1 != vec.end(); ++it1 )
      {
      for( FloatVectorIterator it2 = it1->begin();
           it2 != it1->end(); ++it2 )
        {
        this->Read<float>(f, *it2);
        }
      }
  }

#ifdef __USE_OFFSET_DEBUGGING_CODE__
  void debugOffset(std::ifstream *f, const std::string & msg, long int x)
  {
    std::ostream::pos_type offset(-1);

    offset = f->tellg();
    std::cerr << msg << " " << x << " offset " << offset << std::endl;
    std::cerr.flush();
  }

  void debugOffset(std::ofstream *f, const std::string & msg, long int x)
  {
    std::ostream::pos_type offset(-1);

    offset = f->tellp();
    std::cerr << msg << " " << x << " offset " << offset << std::endl;
    std::cerr.flush();
  }

#endif
  // The templates are arrays of data sets, arrayed by number of angles, and a
  // vector of model intensity values.
  std::map<std::string, Float3DVectorType> m_Templates;
  std::map<std::string, Float2DVectorType> m_TemplateMeans;
  std::map<std::string, bool>              m_TemplateMeansComputed;
  std::map<std::string,
           SImageType::PointType::VectorType> m_RPtoXMean;

  SImageType::PointType::VectorType m_RPtoCECMean;           // CEC is a little
                                                             // different
  SImageType::PointType::VectorType m_CMtoRPMean;
  float                             m_RPPC_to_RPAC_angleMean;
  float                             m_RPAC_over_RPPCMean;
};

#endif // landmarksConstellationModelIO_h
