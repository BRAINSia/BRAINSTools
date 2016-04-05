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

#ifndef landmarksConstellationTrainingDefinitionIO_h
#define landmarksConstellationTrainingDefinitionIO_h
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <itkPoint.h>
#include <itkContinuousIndex.h>
#include "landmarksDataSet.h"
#include "landmarksConstellationModelBase.h"
#include "itkNumberToString.h"


// HACK:  Remove the multiple inheritance here.
class landmarksConstellationTrainingDefinitionIO : public landmarksConstellationModelBase,
  public std::vector<landmarksDataSet>   // This should be a private member
                                         // variable
{
private:
  typedef enum { eof, bad } err_flags;
public:
  landmarksConstellationTrainingDefinitionIO()
  {
  }

  explicit landmarksConstellationTrainingDefinitionIO(const std::string & filename)
  {
    this->ReadFile(filename);
  }

  const ValMapType & GetRadii() const ITK_OVERRIDE
  {
    return landmarksConstellationModelBase::GetRadii();
  }

  template <class type>
  void Read(std::ifstream & s, type & var)
  {
    if( s.bad() )
      {
      throw landmarksConstellationTrainingDefinitionIO::bad;
      }
    if( s.eof() )
      {
      throw landmarksConstellationTrainingDefinitionIO::eof;
      }
    s >> var;
  }

  int ReadFile(const std::string & filename)
  {
    std::ifstream input( filename.c_str() );

    if( input.bad() )
      {
      return -1;
      }
    try
      {
      this->Read<unsigned int>(input, this->m_NumDataSets);
      // Note, more than 2 datasets are needed in order to get valid sample
      // means and variances.
      if( m_NumDataSets < 2 )
        {
        itkGenericExceptionMacro(<< "NumberOfDatasets must be greater than 2");
        }
      this->Read<unsigned int>(input, this->m_SearchboxDims);
      if( m_SearchboxDims % 2 == 0 )
        {
        itkGenericExceptionMacro(<< "SearchBoxDims must be odd");
        }
      this->Read<float>(input, this->m_ResolutionUnits);
      if( m_ResolutionUnits <= 0 )
        {
        itkGenericExceptionMacro(<< "m_ResolutionUnits must be greater than zero");
        }
      this->Read<float>(input, this->m_InitialRotationAngle);
      this->Read<float>(input, this->m_InitialRotationStep);
      this->Read<unsigned int>(input, this->m_NumRotationSteps);

      // Read in template size for each landmark
      std::string name;
      float       val = 0;
      this->Read<std::string>(input, name);

      while( name.compare("END") != 0 )
        {
        this->Read<float>(input, val);
        this->m_Radius[name] = val;
        this->Read<float>(input, val);
        this->m_Height[name] = val;
        this->Read<std::string>(input, name);
        }
      // Read in landmarks
      for( unsigned int i = 0; i < this->m_NumDataSets; i++ )
        {
        std::string TrainingImageFilename;
        this->Read<std::string>(input, TrainingImageFilename);

        std::string LandmarksFilename;
        this->Read<std::string>(input, LandmarksFilename);

        landmarksDataSet DataSet(TrainingImageFilename, LandmarksFilename);
        this->push_back(DataSet);

        std::string Finished;
        this->Read<std::string>(input, Finished);
        // This separator line says 'END'
        }
      input.close();
      }
    catch( err_flags f )
      {
      std::cerr << "File read error "
                << ( f == landmarksConstellationTrainingDefinitionIO::eof ?
           "unexpected end of file" :
           "file read error" )
                << std::endl;
      std::cerr.flush();
      input.close();
      return -1;
      }
    return 0;
  }

private:
};

inline
std::ostream & operator<<(std::ostream & os, const landmarksConstellationTrainingDefinitionIO & def)
{
  itk::NumberToString<double> doubleToString;

  os << def.GetNumDataSets() << std::endl;
  os << def.GetSearchboxDims() << " "
     << def.GetResolutionUnits() << std::endl;
  os << doubleToString(def.GetRadius("RP") ) << " "
     << doubleToString(def.GetHeight("RP") ) << std::endl;
  os << doubleToString(def.GetRadius("AC") ) << " "
     << doubleToString(def.GetHeight("AC") ) << std::endl;
  os << doubleToString(def.GetRadius("PC") ) << " "
     << doubleToString(def.GetHeight("PC") ) << std::endl;
  os << doubleToString(def.GetRadius("VN4") ) << " "
     << doubleToString(def.GetHeight("VN4") ) << std::endl;
  os << doubleToString(def.GetInitialRotationAngle() ) << " "
     << def.GetInitialRotationStep() << " "
     << def.GetNumRotationSteps() << std::endl << std::endl;
  for( unsigned int i = 0; i < def.GetNumDataSets(); i++ )
    {
    os << def[i] << std::endl;
    }
  return os;
}

#endif // landmarksConstellationTrainingDefinitionIO_h
